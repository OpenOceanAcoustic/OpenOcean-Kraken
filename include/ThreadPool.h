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

/**
 * @brief 高性能线程池类，支持普通任务、带ID任务、等待完成、超时等功能
 */
class ThreadPool {
public:
    explicit ThreadPool(size_t threads);
    
    /**
     * @brief 提交一个任务到线程池
     * @param f 函数对象
     * @param args 参数
     * @return std::future<T> 可用于获取返回值或异常
     */
    template<class F, class... Args>
    auto enqueue(F&& f, Args&&... args) 
        -> std::future<std::invoke_result_t<F, Args...>>;

    void wait_completion();
    
    bool wait_completion_for(const std::chrono::milliseconds& timeout);
    
    template<class F, class... Args>
    auto enqueue_with_id(const std::string& id, F&& f, Args&&... args) 
        -> std::future<std::invoke_result_t<F, Args...>>;

    void wait_id(const std::string& id);
    
    bool wait_id_for(const std::string& id, const std::chrono::milliseconds& timeout);
    
    void cleanup_finished_ids();
    
    ~ThreadPool();

private:
    std::vector<std::thread> workers;
    std::queue<std::function<void()>> tasks;
    std::atomic<size_t> active_tasks{0};
    std::unordered_map<std::string, std::atomic<size_t>> id_task_counts;
    mutable std::mutex queue_mutex;
    mutable std::mutex id_map_mutex;
    std::condition_variable condition;
    std::condition_variable id_condition;
    std::atomic<bool> stop{false};

    // 辅助函数：解包 tuple 并调用函数
    template<typename Func, typename Tuple, std::size_t... I>
    static auto invoke_impl(Func&& func, Tuple&& t, std::index_sequence<I...>) 
        -> decltype(std::invoke(std::forward<Func>(func), std::get<I>(std::forward<Tuple>(t))...))
    {
        return std::invoke(std::forward<Func>(func), std::get<I>(std::forward<Tuple>(t))...);
    }

    template<typename Func, typename Tuple>
    static auto invoke_tuple(Func&& func, Tuple&& t)
        -> decltype(invoke_impl(std::forward<Func>(func), std::forward<Tuple>(t),
                                std::make_index_sequence<std::tuple_size_v<std::decay_t<Tuple>>>{}))
    {
        constexpr auto N = std::tuple_size_v<std::decay_t<Tuple>>;
        return invoke_impl(std::forward<Func>(func), std::forward<Tuple>(t),
                          std::make_index_sequence<N>{});
    }
};

// 构造函数（未变）
inline ThreadPool::ThreadPool(size_t threads)
    : stop(false)
{
    if (threads == 0) {
        throw std::invalid_argument("ThreadPool constructor: threads must be > 0");
    }

    for (size_t i = 0; i < threads; ++i) {
        workers.emplace_back([this] {
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
            }
        });
    }
}

// 提交普通任务（✅ 已改为 C++17 兼容写法）
template<class F, class... Args>
auto ThreadPool::enqueue(F&& f, Args&&... args) 
    -> std::future<std::invoke_result_t<F, Args...>>
{
    using return_type = std::invoke_result_t<F, Args...>;

    // 使用 tuple 存储参数，避免 C++20 的 pack capture
    auto bound_task = std::make_shared<std::packaged_task<return_type()>>(
        [func = std::forward<F>(f), args_tuple = std::make_tuple(std::forward<Args>(args)...)]() mutable -> return_type {
            return invoke_tuple(func, args_tuple);
        }
    );

    std::future<return_type> res = bound_task->get_future();
    {
        std::unique_lock<std::mutex> lock(queue_mutex);
        if (stop.load(std::memory_order_relaxed)) {
            throw std::runtime_error("enqueue on stopped ThreadPool");
        }

        active_tasks.fetch_add(1, std::memory_order_relaxed);

        tasks.emplace([bound_task, this]() {
            (*bound_task)();

            size_t old_active;
            do {
                old_active = active_tasks.load(std::memory_order_relaxed);
                if (old_active == 0) {
                    std::cerr << "Error: active_tasks was 0 when trying to decrement!\n";
                    break;
                }
            } while (!active_tasks.compare_exchange_weak(
                old_active, old_active - 1,
                std::memory_order_release,
                std::memory_order_relaxed));

            if (old_active == 1) {
                condition.notify_all();
            }
        });
    }
    condition.notify_one();
    return res;
}

// 等待所有任务完成（未变）
inline void ThreadPool::wait_completion() {
    std::unique_lock<std::mutex> lock(queue_mutex);
    condition.wait(lock, [this] {
        return tasks.empty() && active_tasks.load(std::memory_order_acquire) == 0;
    });
}

// 带超时等待（未变）
inline bool ThreadPool::wait_completion_for(const std::chrono::milliseconds& timeout) {
    std::unique_lock<std::mutex> lock(queue_mutex);
    return condition.wait_for(lock, timeout, [this] {
        return tasks.empty() && active_tasks.load(std::memory_order_acquire) == 0;
    });
}

// 提交带 ID 的任务（✅ 已改为 C++17 兼容写法）
template<class F, class... Args>
auto ThreadPool::enqueue_with_id(const std::string& id, F&& f, Args&&... args) 
    -> std::future<std::invoke_result_t<F, Args...>>
{
    using return_type = std::invoke_result_t<F, Args...>;

    auto bound_task = std::make_shared<std::packaged_task<return_type()>>(
        [func = std::forward<F>(f), args_tuple = std::make_tuple(std::forward<Args>(args)...)]() mutable -> return_type {
            return invoke_tuple(func, args_tuple);
        }
    );

    std::future<return_type> res = bound_task->get_future();
    {
        std::unique_lock<std::mutex> lock(queue_mutex);
        if (stop.load(std::memory_order_relaxed)) {
            throw std::runtime_error("enqueue with ID on stopped ThreadPool");
        }

        active_tasks.fetch_add(1, std::memory_order_relaxed);

        {
            std::unique_lock<std::mutex> id_lock(id_map_mutex);
            ++id_task_counts[id];
        }

        tasks.emplace([bound_task, this, id]() {
            (*bound_task)();

            size_t old_active;
            do {
                old_active = active_tasks.load(std::memory_order_relaxed);
                if (old_active == 0) {
                    std::cerr << "Error: active_tasks was 0 when trying to decrement in ID task!\n";
                    break;
                }
            } while (!active_tasks.compare_exchange_weak(
                old_active, old_active - 1,
                std::memory_order_release,
                std::memory_order_relaxed));

            if (old_active == 1) {
                condition.notify_all();
            }

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
        });
    }
    condition.notify_one();
    return res;
}

// 等待指定 ID 的任务（未变）
inline void ThreadPool::wait_id(const std::string& id) {
    std::unique_lock<std::mutex> lock(id_map_mutex);
    id_condition.wait(lock, [this, &id] {
        auto it = id_task_counts.find(id);
        return it == id_task_counts.end() || it->second.load(std::memory_order_acquire) == 0;
    });
}

// 带超时等待 ID（未变）
inline bool ThreadPool::wait_id_for(const std::string& id, const std::chrono::milliseconds& timeout) {
    std::unique_lock<std::mutex> lock(id_map_mutex);
    bool result = id_condition.wait_for(lock, timeout, [this, &id] {
        auto it = id_task_counts.find(id);
        return it == id_task_counts.end() || it->second.load(std::memory_order_acquire) == 0;
    });

    if (!result) {
        std::cerr << "ThreadPool: Timeout waiting for ID '" << id << "'\n";
    }
    return result;
}

// 清理已完成的任务 ID（未变）
inline void ThreadPool::cleanup_finished_ids() {
    std::unique_lock<std::mutex> lock(id_map_mutex);
    for (auto it = id_task_counts.begin(); it != id_task_counts.end();) {
        if (it->second.load(std::memory_order_acquire) == 0) {
            it = id_task_counts.erase(it);
        } else {
            ++it;
        }
    }
}

// 析构函数（未变）
inline ThreadPool::~ThreadPool() {
    {
        std::lock_guard<std::mutex> lock(queue_mutex);
        stop.store(true, std::memory_order_relaxed);
        condition.notify_all();
    }
    {
        std::lock_guard<std::mutex> lock(id_map_mutex);
        id_condition.notify_all();
    }

    for (std::thread& worker : workers) {
        if (worker.joinable()) {
            worker.join();
        }
    }

    cleanup_finished_ids();
}

#endif // THREAD_POOL_H