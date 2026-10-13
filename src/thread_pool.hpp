#pragma once

#include <atomic>
#include <condition_variable>
#include <cstdint>
#include <deque>
#include <exception>
#include <functional>
#include <memory>
#include <mutex>
#include <thread>
#include <vector>

namespace wqfm {

// Fork-join loops and tasks on a fixed set of threads, nestable and callable from any thread.
//
// parallel_for(n, fn) calls fn(i, slot) for i in [0, n): the caller takes part, idle threads join in,
// and items are claimed dynamically (slow cores and uneven items do not stall it). `slot` identifies
// the executing thread (0 = the thread that created the pool), for per-thread scratch space; items
// must not block. Several parallel_for calls may be active at once.
//
// submit(fn) queues a task; wait(task) runs queued tasks and loop items while it waits, so a thread
// never idles while there is work (fork-join recursion: run one child inline, submit the other).
class ThreadPool {
public:
    struct Task {
        std::function<void()> fn;
        std::atomic<bool> done{false};
        std::exception_ptr error;
    };

    explicit ThreadPool(int threads);
    ~ThreadPool();
    ThreadPool(const ThreadPool&) = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;

    int threads() const { return static_cast<int>(workers_.size()) + 1; }

    void parallel_for(int n, const std::function<void(int, int)>& fn);
    std::shared_ptr<Task> submit(std::function<void()> fn);
    void wait(const std::shared_ptr<Task>& task);

private:
    struct Job {
        const std::function<void(int, int)>* fn = nullptr;
        int n = 0;
        std::atomic<int> next{0}, done{0}, helpers{0};
        std::atomic<bool> failed{false};
        std::exception_ptr error;
    };

    void worker(int slot);
    void help(Job& job);
    Job* claimable_job();                // under the lock
    void run(const std::shared_ptr<Task>& task);
    bool help_once(bool tasks_too);      // one loop item or one task, if any

    std::vector<std::thread> workers_;
    std::mutex m_;
    std::condition_variable work_cv_;    // idle threads: new items or tasks, or stop
    std::condition_variable done_cv_;    // loop callers and task waiters
    std::vector<Job*> jobs_;
    std::deque<std::shared_ptr<Task>> tasks_;
    bool stop_ = false;
};

}  // namespace wqfm
