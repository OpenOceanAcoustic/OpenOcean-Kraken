#ifndef OPENOCEAN_KRAKENC_THREAD_POOL_H
#define OPENOCEAN_KRAKENC_THREAD_POOL_H

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <functional>
#include <future>
#include <iostream>
#include <memory>
#include <mutex>
#include <queue>
#include <stdexcept>
#include <string>
#include <thread>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace OpenOceanKrakenc {

/**
 * A fixed-size executor.  ThreadPool itself is a handle; ExecutorLease keeps the
 * internal executor alive while an operation is scheduling or awaiting work.
 */
class ThreadPool
{
private:
    struct ExecutorState;
    struct ExecutorControl;

public:
    class ExecutorLease
    {
    public:
        ExecutorLease() = default;

        explicit operator bool() const noexcept { return static_cast<bool>(control_); }
        bool ownsCurrentThread() const noexcept;

        template <class F, class... Args>
        auto enqueue(F &&f, Args &&...args)
            -> std::future<std::invoke_result_t<F, Args...>>;

        template <class F, class... Args>
        auto enqueue_with_id(const std::string &id, F &&f, Args &&...args)
            -> std::future<std::invoke_result_t<F, Args...>>;

        void wait_completion() const;
        bool wait_completion_for(const std::chrono::milliseconds &timeout) const;
        void wait_id(const std::string &id) const;
        bool wait_id_for(const std::string &id,
                         const std::chrono::milliseconds &timeout) const;

    private:
        friend class ThreadPool;
        explicit ExecutorLease(std::shared_ptr<ExecutorControl> control)
            : control_(std::move(control)) {}

        std::shared_ptr<ExecutorControl> control_;
    };

    class BorrowedExecutor
    {
    public:
        BorrowedExecutor() = default;
        bool expired() const noexcept { return control_.expired(); }
        ExecutorLease lock() const noexcept
        {
            return ExecutorLease(control_.lock());
        }

    private:
        friend class ThreadPool;
        explicit BorrowedExecutor(const std::shared_ptr<ExecutorControl> &control)
            : control_(control) {}

        std::weak_ptr<ExecutorControl> control_;
    };

    explicit ThreadPool(std::size_t threads);
    ~ThreadPool() = default;

    ThreadPool(const ThreadPool &) = delete;
    ThreadPool &operator=(const ThreadPool &) = delete;
    ThreadPool(ThreadPool &&) = delete;
    ThreadPool &operator=(ThreadPool &&) = delete;

    ExecutorLease acquireExecutor() const noexcept { return ExecutorLease(control_); }
    BorrowedExecutor borrowExecutor() const noexcept { return BorrowedExecutor(control_); }

    template <class F, class... Args>
    auto enqueue(F &&f, Args &&...args)
        -> std::future<std::invoke_result_t<F, Args...>>
    {
        return acquireExecutor().enqueue(std::forward<F>(f),
                                         std::forward<Args>(args)...);
    }

    template <class F, class... Args>
    auto enqueue_with_id(const std::string &id, F &&f, Args &&...args)
        -> std::future<std::invoke_result_t<F, Args...>>
    {
        return acquireExecutor().enqueue_with_id(
            id, std::forward<F>(f), std::forward<Args>(args)...);
    }

    void wait_completion() { acquireExecutor().wait_completion(); }
    bool wait_completion_for(const std::chrono::milliseconds &timeout)
    {
        return acquireExecutor().wait_completion_for(timeout);
    }
    bool ownsCurrentThread() const noexcept
    {
        return acquireExecutor().ownsCurrentThread();
    }
    void wait_id(const std::string &id) { acquireExecutor().wait_id(id); }
    bool wait_id_for(const std::string &id,
                     const std::chrono::milliseconds &timeout)
    {
        return acquireExecutor().wait_id_for(id, timeout);
    }
    void cleanup_finished_ids();

private:
    template <typename Func, typename Tuple, std::size_t... I>
    static auto invoke_impl(Func &&func, Tuple &&tuple, std::index_sequence<I...>)
        -> decltype(std::invoke(std::forward<Func>(func),
                                std::get<I>(std::forward<Tuple>(tuple))...))
    {
        return std::invoke(std::forward<Func>(func),
                           std::get<I>(std::forward<Tuple>(tuple))...);
    }

    template <typename Func, typename Tuple>
    static auto invoke_tuple(Func &&func, Tuple &&tuple)
        -> decltype(invoke_impl(
            std::forward<Func>(func), std::forward<Tuple>(tuple),
            std::make_index_sequence<std::tuple_size_v<std::decay_t<Tuple>>>{}))
    {
        constexpr auto count = std::tuple_size_v<std::decay_t<Tuple>>;
        return invoke_impl(std::forward<Func>(func), std::forward<Tuple>(tuple),
                           std::make_index_sequence<count>{});
    }

    struct ExecutorState : std::enable_shared_from_this<ExecutorState>
    {
        std::vector<std::thread> workers;
        std::queue<std::function<void()>> tasks;
        std::atomic<std::size_t> activeTasks{0};
        std::unordered_map<std::string, std::size_t> idTaskCounts;
        mutable std::mutex queueMutex;
        mutable std::mutex idMapMutex;
        std::condition_variable condition;
        std::condition_variable idCondition;
        std::atomic<bool> stop{false};
        std::atomic<bool> shutdownStarted{false};

        void start(std::size_t threads);
        void shutdown() noexcept;
        void push(std::function<void()> task);
        void pushIdentified(const std::string &id, std::function<void()> task);
    };

    struct ExecutorControl
    {
        explicit ExecutorControl(std::shared_ptr<ExecutorState> value)
            : state(std::move(value)) {}
        ~ExecutorControl() { state->shutdown(); }
        std::shared_ptr<ExecutorState> state;
    };

    inline static thread_local const ExecutorState *currentExecutor_ = nullptr;
    inline static thread_local const std::string *currentTaskId_ = nullptr;
    std::shared_ptr<ExecutorControl> control_;
};

inline ThreadPool::ThreadPool(std::size_t threads)
{
    if (threads == 0)
    {
        throw std::invalid_argument("ThreadPool constructor: threads must be > 0");
    }
    auto state = std::make_shared<ExecutorState>();
    state->start(threads);
    control_ = std::make_shared<ExecutorControl>(std::move(state));
}

inline void ThreadPool::ExecutorState::start(std::size_t threadCount)
{
    try
    {
        const std::shared_ptr<ExecutorState> self = shared_from_this();
        for (std::size_t i = 0; i < threadCount; ++i)
        {
            workers.emplace_back([self] {
                currentExecutor_ = self.get();
                while (true)
                {
                    std::function<void()> task;
                    {
                        std::unique_lock<std::mutex> lock(self->queueMutex);
                        self->condition.wait(lock, [&] {
                            return self->stop.load(std::memory_order_acquire) ||
                                   !self->tasks.empty();
                        });
                        if (self->stop.load(std::memory_order_acquire) &&
                            self->tasks.empty())
                        {
                            currentExecutor_ = nullptr;
                            return;
                        }
                        task = std::move(self->tasks.front());
                        self->tasks.pop();
                    }
                    task();
                    const std::size_t previous =
                        self->activeTasks.fetch_sub(1, std::memory_order_acq_rel);
                    if (previous == 1)
                    {
                        self->condition.notify_all();
                    }
                }
            });
        }
    }
    catch (...)
    {
        shutdown();
        throw;
    }
}

inline void ThreadPool::ExecutorState::shutdown() noexcept
{
    bool expected = false;
    if (!shutdownStarted.compare_exchange_strong(expected, true,
                                                 std::memory_order_acq_rel))
    {
        return;
    }
    {
        std::lock_guard<std::mutex> queueLock(queueMutex);
        stop.store(true, std::memory_order_release);
    }
    condition.notify_all();
    idCondition.notify_all();
    for (std::thread &worker : workers)
    {
        if (!worker.joinable())
        {
            continue;
        }
        if (worker.get_id() == std::this_thread::get_id())
        {
            worker.detach();
        }
        else
        {
            worker.join();
        }
    }
}

inline void ThreadPool::ExecutorState::push(std::function<void()> task)
{
    {
        std::lock_guard<std::mutex> lock(queueMutex);
        if (stop.load(std::memory_order_acquire))
        {
            throw std::runtime_error("enqueue on stopped ThreadPool");
        }
        activeTasks.fetch_add(1, std::memory_order_relaxed);
        tasks.push(std::move(task));
    }
    condition.notify_one();
}

inline void ThreadPool::ExecutorState::pushIdentified(
    const std::string &id, std::function<void()> task)
{
    {
        std::lock_guard<std::mutex> lock(queueMutex);
        if (stop.load(std::memory_order_acquire))
        {
            throw std::runtime_error("enqueue with ID on stopped ThreadPool");
        }
        {
            std::lock_guard<std::mutex> idLock(idMapMutex);
            ++idTaskCounts[id];
        }
        activeTasks.fetch_add(1, std::memory_order_relaxed);
        tasks.push([self = shared_from_this(), id, task = std::move(task)]() mutable {
            const std::string *previousTaskId = currentTaskId_;
            currentTaskId_ = &id;
            try
            {
                task();
            }
            catch (...)
            {
                currentTaskId_ = previousTaskId;
                throw;
            }
            currentTaskId_ = previousTaskId;
            bool notify = false;
            {
                std::lock_guard<std::mutex> idLock(self->idMapMutex);
                const auto it = self->idTaskCounts.find(id);
                if (it != self->idTaskCounts.end() && --it->second == 0)
                {
                    self->idTaskCounts.erase(it);
                    notify = true;
                }
            }
            if (notify)
            {
                self->idCondition.notify_all();
            }
        });
    }
    condition.notify_one();
}

inline bool ThreadPool::ExecutorLease::ownsCurrentThread() const noexcept
{
    return control_ && currentExecutor_ == control_->state.get();
}

template <class F, class... Args>
auto ThreadPool::ExecutorLease::enqueue(F &&f, Args &&...args)
    -> std::future<std::invoke_result_t<F, Args...>>
{
    if (!control_)
    {
        throw std::runtime_error("external ThreadPool expired");
    }
    using Return = std::invoke_result_t<F, Args...>;
    auto packaged = std::make_shared<std::packaged_task<Return()>>(
        [function = std::forward<F>(f),
         tuple = std::make_tuple(std::forward<Args>(args)...)]() mutable -> Return {
            return ThreadPool::invoke_tuple(function, tuple);
        });
    std::future<Return> future = packaged->get_future();
    control_->state->push([packaged] { (*packaged)(); });
    return future;
}

template <class F, class... Args>
auto ThreadPool::ExecutorLease::enqueue_with_id(
    const std::string &id, F &&f, Args &&...args)
    -> std::future<std::invoke_result_t<F, Args...>>
{
    if (!control_)
    {
        throw std::runtime_error("external ThreadPool expired");
    }
    using Return = std::invoke_result_t<F, Args...>;
    auto packaged = std::make_shared<std::packaged_task<Return()>>(
        [function = std::forward<F>(f),
         tuple = std::make_tuple(std::forward<Args>(args)...)]() mutable -> Return {
            return ThreadPool::invoke_tuple(function, tuple);
        });
    std::future<Return> future = packaged->get_future();
    control_->state->pushIdentified(id, [packaged] { (*packaged)(); });
    return future;
}

inline void ThreadPool::ExecutorLease::wait_completion() const
{
    if (!control_)
    {
        throw std::runtime_error("external ThreadPool expired");
    }
    if (ownsCurrentThread())
    {
        throw std::logic_error(
            "ThreadPool current task cannot wait for executor completion");
    }
    std::unique_lock<std::mutex> lock(control_->state->queueMutex);
    control_->state->condition.wait(lock, [&] {
        return control_->state->tasks.empty() &&
               control_->state->activeTasks.load(std::memory_order_acquire) == 0;
    });
}

inline bool ThreadPool::ExecutorLease::wait_completion_for(
    const std::chrono::milliseconds &timeout) const
{
    if (!control_)
    {
        throw std::runtime_error("external ThreadPool expired");
    }
    if (ownsCurrentThread())
    {
        return false;
    }
    std::unique_lock<std::mutex> lock(control_->state->queueMutex);
    return control_->state->condition.wait_for(lock, timeout, [&] {
        return control_->state->tasks.empty() &&
               control_->state->activeTasks.load(std::memory_order_acquire) == 0;
    });
}

inline void ThreadPool::ExecutorLease::wait_id(const std::string &id) const
{
    if (!control_)
    {
        throw std::runtime_error("external ThreadPool expired");
    }
    if (ownsCurrentThread() && currentTaskId_ &&
        *currentTaskId_ == id)
    {
        throw std::logic_error(
            "ThreadPool current task ID cannot wait for itself");
    }
    std::unique_lock<std::mutex> lock(control_->state->idMapMutex);
    control_->state->idCondition.wait(lock, [&] {
        return control_->state->idTaskCounts.find(id) ==
               control_->state->idTaskCounts.end();
    });
}

inline bool ThreadPool::ExecutorLease::wait_id_for(
    const std::string &id, const std::chrono::milliseconds &timeout) const
{
    if (!control_)
    {
        throw std::runtime_error("external ThreadPool expired");
    }
    if (ownsCurrentThread() && currentTaskId_ &&
        *currentTaskId_ == id)
    {
        return false;
    }
    std::unique_lock<std::mutex> lock(control_->state->idMapMutex);
    const bool complete = control_->state->idCondition.wait_for(lock, timeout, [&] {
        return control_->state->idTaskCounts.find(id) ==
               control_->state->idTaskCounts.end();
    });
    if (!complete)
    {
        std::cerr << "ThreadPool: Timeout waiting for ID '" << id << "'\n";
    }
    return complete;
}

inline void ThreadPool::cleanup_finished_ids()
{
    const ExecutorLease lease = acquireExecutor();
    std::lock_guard<std::mutex> lock(lease.control_->state->idMapMutex);
    for (auto it = lease.control_->state->idTaskCounts.begin();
         it != lease.control_->state->idTaskCounts.end();)
    {
        if (it->second == 0)
        {
            it = lease.control_->state->idTaskCounts.erase(it);
        }
        else
        {
            ++it;
        }
    }
}

} // namespace OpenOceanKrakenc

#endif // OPENOCEAN_KRAKENC_THREAD_POOL_H
