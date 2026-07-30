#ifndef THREAD_POOL_H
#define THREAD_POOL_H

#include <vector>
#include <queue>
#include <memory>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <future>
#include <functional>
#include <stdexcept>
#include <unordered_map>
#include <atomic>
#include <chrono>
#include <iostream>

#if defined(OPENOCEAN_KRAKEN_FIELD_RUNTIME)
#include <openocean/field/runtime.hpp>
#endif

/**
 * @brief 高性能线程池类，支持普通任务、带ID任务、等待完成、超时等功能
 */
class ThreadPool
{
public:
    explicit ThreadPool(size_t threads);
#if defined(OPENOCEAN_KRAKEN_FIELD_RUNTIME)
    explicit ThreadPool(OpenOcean::Field::FieldRuntime &runtime);
#endif

    /**
     * @brief 提交一个任务到线程池
     * @param f 函数对象
     * @param args 参数
     * @return std::future<T> 可用于获取返回值或异常
     */
    template <class F, class... Args>
    auto enqueue(F &&f, Args &&...args)
        -> std::future<std::invoke_result_t<F, Args...>>;

    void wait_completion();

    bool wait_completion_for(const std::chrono::milliseconds &timeout);

    template <class F, class... Args>
    auto enqueue_with_id(const std::string &id, F &&f, Args &&...args)
        -> std::future<std::invoke_result_t<F, Args...>>;

    void wait_id(const std::string &id);

    bool wait_id_for(const std::string &id, const std::chrono::milliseconds &timeout);

    void cleanup_finished_ids();

    ~ThreadPool();

private:
    // 辅助：解包 tuple 并调用函数（C++17 兼容）
    template <typename Func, typename Tuple, std::size_t... I>
    static auto invoke_impl(Func &&func, Tuple &&t, std::index_sequence<I...>)
        -> decltype(std::invoke(std::forward<Func>(func), std::get<I>(std::forward<Tuple>(t))...))
    {
        return std::invoke(std::forward<Func>(func), std::get<I>(std::forward<Tuple>(t))...);
    }

    template <typename Func, typename Tuple>
    static auto invoke_tuple(Func &&func, Tuple &&t)
        -> decltype(invoke_impl(std::forward<Func>(func), std::forward<Tuple>(t),
                                std::make_index_sequence<std::tuple_size_v<std::decay_t<Tuple>>>{}))
    {
        constexpr auto N = std::tuple_size_v<std::decay_t<Tuple>>;
        return invoke_impl(std::forward<Func>(func), std::forward<Tuple>(t),
                           std::make_index_sequence<N>{});
    }

    std::vector<std::thread> workers;
    std::queue<std::function<void()>> tasks;
    std::atomic<size_t> active_tasks{0};
    std::unordered_map<std::string, size_t> id_task_counts; // 改为非 atomic，由 mutex 保护
    mutable std::mutex queue_mutex;
    mutable std::mutex id_map_mutex;
    std::condition_variable condition;
    std::condition_variable id_condition;
    std::atomic<bool> stop{false};
#if defined(OPENOCEAN_KRAKEN_FIELD_RUNTIME)
    OpenOcean::Field::FieldRuntime *fieldRuntime = nullptr;
#endif
};

// ==================== 实现 ====================

inline ThreadPool::ThreadPool(size_t threads)
    : stop(false)
{
    if (threads == 0)
    {
        throw std::invalid_argument("ThreadPool constructor: threads must be > 0");
    }

    for (size_t i = 0; i < threads; ++i)
    {
        workers.emplace_back([this]
                             {
            while (true) {
                std::function<void()> task;
                {
                    std::unique_lock<std::mutex> lock(queue_mutex);
                    condition.wait(lock, [this] { 
                        return stop.load(std::memory_order_relaxed) || !tasks.empty(); 
                    });
                    if (stop.load(std::memory_order_relaxed) && tasks.empty()) {
                        return;
                    }
                    task = std::move(tasks.front());
                    tasks.pop();
                }
                task(); // 执行任务

                // ✅ 修复：使用 fetch_sub 安全 decrement
                size_t prev = active_tasks.fetch_sub(1, std::memory_order_acq_rel);
                if (prev == 1) {
                    condition.notify_all();
                }
            } });
    }
}

#if defined(OPENOCEAN_KRAKEN_FIELD_RUNTIME)
inline ThreadPool::ThreadPool(
    OpenOcean::Field::FieldRuntime &runtime)
    : stop(false),
      fieldRuntime(&runtime)
{
}
#endif

template <class F, class... Args>
auto ThreadPool::enqueue(F &&f, Args &&...args)
    -> std::future<std::invoke_result_t<F, Args...>>
{
    using return_type = std::invoke_result_t<F, Args...>;

    auto bound_task = std::make_shared<std::packaged_task<return_type()>>(
        [func = std::forward<F>(f), args_tuple = std::make_tuple(std::forward<Args>(args)...)]() mutable -> return_type
        {
            return invoke_tuple(func, args_tuple);
        });

    std::future<return_type> res = bound_task->get_future();
    std::function<void()> queued = [bound_task, this]()
    {
        (*bound_task)();
        size_t prev =
            active_tasks.fetch_sub(1, std::memory_order_acq_rel);
        if (prev == 1)
        {
            condition.notify_all();
        }
    };
    bool useFieldRuntime = false;
    {
        std::unique_lock<std::mutex> lock(queue_mutex);
        if (stop.load(std::memory_order_relaxed))
        {
            throw std::runtime_error("enqueue on stopped ThreadPool");
        }

        active_tasks.fetch_add(1, std::memory_order_relaxed);
#if defined(OPENOCEAN_KRAKEN_FIELD_RUNTIME)
        if (fieldRuntime != nullptr)
        {
            useFieldRuntime = true;
        }
        else
#endif
        {
            tasks.emplace([bound_task]()
                          { (*bound_task)(); });
        }
    }
#if defined(OPENOCEAN_KRAKEN_FIELD_RUNTIME)
    if (useFieldRuntime)
    {
        try
        {
            fieldRuntime->submit(std::move(queued));
        }
        catch (...)
        {
            active_tasks.fetch_sub(1, std::memory_order_acq_rel);
            condition.notify_all();
            throw;
        }
    }
#endif
#if defined(OPENOCEAN_KRAKEN_FIELD_RUNTIME)
    if (fieldRuntime == nullptr)
#endif
    {
        condition.notify_one();
    }
    return res;
}

inline void ThreadPool::wait_completion()
{
    std::unique_lock<std::mutex> lock(queue_mutex);
    condition.wait(lock, [this]
                   { return tasks.empty() && active_tasks.load(std::memory_order_acquire) == 0; });
}

inline bool ThreadPool::wait_completion_for(const std::chrono::milliseconds &timeout)
{
    std::unique_lock<std::mutex> lock(queue_mutex);
    return condition.wait_for(lock, timeout, [this]
                              { return tasks.empty() && active_tasks.load(std::memory_order_acquire) == 0; });
}

template <class F, class... Args>
auto ThreadPool::enqueue_with_id(const std::string &id, F &&f, Args &&...args)
    -> std::future<std::invoke_result_t<F, Args...>>
{
    using return_type = std::invoke_result_t<F, Args...>;

    auto bound_task = std::make_shared<std::packaged_task<return_type()>>(
        [func = std::forward<F>(f), args_tuple = std::make_tuple(std::forward<Args>(args)...)]() mutable -> return_type
        {
            return invoke_tuple(func, args_tuple);
        });

    std::future<return_type> res = bound_task->get_future();
    bool useFieldRuntime = false;
    std::function<void()> identifiedTask;
    {
        std::unique_lock<std::mutex> lock(queue_mutex);
        if (stop.load(std::memory_order_relaxed))
        {
            throw std::runtime_error("enqueue with ID on stopped ThreadPool");
        }

        active_tasks.fetch_add(1, std::memory_order_relaxed);

        {
            std::unique_lock<std::mutex> id_lock(id_map_mutex);
            id_task_counts[id]++; // ✅ 现在是普通 size_t，由 mutex 保护
        }

        identifiedTask = [bound_task, this, id]()
        {
            (*bound_task)();

            // ✅ 安全 decrement global counter
            size_t prev = active_tasks.fetch_sub(1, std::memory_order_acq_rel);
            if (prev == 1) {
                condition.notify_all();
            }

            // ✅ 更新 ID 计数（加锁）
            bool should_notify = false;
            {
                std::unique_lock<std::mutex> lock(id_map_mutex);
                auto it = id_task_counts.find(id);
                if (it != id_task_counts.end()) {
                    if (--(it->second) == 0) {
                        id_task_counts.erase(it);
                        should_notify = true;
                    }
                }
            }
            if (should_notify) {
                id_condition.notify_all();
            }
        };
#if defined(OPENOCEAN_KRAKEN_FIELD_RUNTIME)
        if (fieldRuntime != nullptr)
        {
            useFieldRuntime = true;
        }
        else
#endif
        {
            tasks.emplace(std::move(identifiedTask));
        }
    }
#if defined(OPENOCEAN_KRAKEN_FIELD_RUNTIME)
    if (useFieldRuntime)
    {
        try
        {
            fieldRuntime->submit(std::move(identifiedTask));
        }
        catch (...)
        {
            active_tasks.fetch_sub(1, std::memory_order_acq_rel);
            {
                std::unique_lock<std::mutex> id_lock(id_map_mutex);
                auto match = id_task_counts.find(id);
                if (match != id_task_counts.end() &&
                    --match->second == 0)
                {
                    id_task_counts.erase(match);
                }
            }
            condition.notify_all();
            id_condition.notify_all();
            throw;
        }
    }
#endif
#if defined(OPENOCEAN_KRAKEN_FIELD_RUNTIME)
    if (fieldRuntime == nullptr)
#endif
    {
        condition.notify_one();
    }
    return res;
}

inline void ThreadPool::wait_id(const std::string &id)
{
    std::unique_lock<std::mutex> lock(id_map_mutex);
    id_condition.wait(lock, [this, &id]
                      {
                          auto it = id_task_counts.find(id);
                          return it == id_task_counts.end(); // 只要不存在，说明已完成
                      });
}

inline bool ThreadPool::wait_id_for(const std::string &id, const std::chrono::milliseconds &timeout)
{
    std::unique_lock<std::mutex> lock(id_map_mutex);
    bool result = id_condition.wait_for(lock, timeout, [this, &id]
                                        {
        auto it = id_task_counts.find(id);
        return it == id_task_counts.end(); });

    if (!result)
    {
        std::cerr << "ThreadPool: Timeout waiting for ID '" << id << "'\n";
    }
    return result;
}

inline void ThreadPool::cleanup_finished_ids()
{
    std::unique_lock<std::mutex> lock(id_map_mutex);
    for (auto it = id_task_counts.begin(); it != id_task_counts.end();)
    {
        if (it->second == 0)
        {
            it = id_task_counts.erase(it);
        }
        else
        {
            ++it;
        }
    }
}

inline ThreadPool::~ThreadPool()
{
#if defined(OPENOCEAN_KRAKEN_FIELD_RUNTIME)
    if (fieldRuntime != nullptr)
    {
        wait_completion();
        return;
    }
#endif
    {
        std::lock_guard<std::mutex> lock(queue_mutex);
        stop.store(true, std::memory_order_relaxed);
        condition.notify_all();
    }
    {
        std::lock_guard<std::mutex> lock(id_map_mutex);
        id_condition.notify_all();
    }

    for (std::thread &worker : workers)
    {
        if (worker.joinable())
        {
            worker.join();
        }
    }

    cleanup_finished_ids();
}

#endif // THREAD_POOL_H
