#include "thread_pool.hpp"

#include <algorithm>

namespace wqfm {
namespace {
thread_local int tls_slot = 0;   // 0 for the creating (and any outside) thread, 1.. for workers
}

ThreadPool::ThreadPool(int threads) {
    for (int slot = 1; slot < threads; ++slot) workers_.emplace_back([this, slot] { worker(slot); });
}

ThreadPool::~ThreadPool() {
    {
        std::lock_guard<std::mutex> lk(m_);
        stop_ = true;
    }
    work_cv_.notify_all();
    for (auto& t : workers_) t.join();
}

ThreadPool::Job* ThreadPool::claimable_job() {
    for (Job* j : jobs_)
        if (j->next.load() < j->n) return j;
    return nullptr;
}

void ThreadPool::worker(int slot) {
    tls_slot = slot;
    std::unique_lock<std::mutex> lk(m_);
    for (;;) {
        work_cv_.wait(lk, [&] { return stop_ || claimable_job() || !tasks_.empty(); });
        if (Job* job = claimable_job()) {
            ++job->helpers;                       // under the lock: the caller cannot retire the job meanwhile
            lk.unlock();
            help(*job);
            lk.lock();
            if (--job->helpers == 0) done_cv_.notify_all();
            continue;
        }
        if (!tasks_.empty()) {
            std::shared_ptr<Task> task = std::move(tasks_.front());
            tasks_.pop_front();
            lk.unlock();
            run(task);
            lk.lock();
            continue;
        }
        if (stop_) return;
    }
}

void ThreadPool::help(Job& job) {
    for (int i; (i = job.next.fetch_add(1)) < job.n;) {
        if (!job.failed.load()) {
            try {
                (*job.fn)(i, tls_slot);
            } catch (...) {
                if (!job.failed.exchange(true)) job.error = std::current_exception();
            }
        }
        if (job.done.fetch_add(1) + 1 == job.n) {
            std::lock_guard<std::mutex> lk(m_);
            done_cv_.notify_all();
        }
    }
}

void ThreadPool::parallel_for(int n, const std::function<void(int, int)>& fn) {
    if (workers_.empty() || n <= 1) {
        for (int i = 0; i < n; ++i) fn(i, tls_slot);
        return;
    }
    Job job;
    job.fn = &fn;
    job.n = n;
    {
        std::lock_guard<std::mutex> lk(m_);
        jobs_.push_back(&job);
    }
    work_cv_.notify_all();
    done_cv_.notify_all();                // threads waiting for a task help too
    help(job);
    {
        // in-flight items of other threads are short: just wait for them
        std::unique_lock<std::mutex> lk(m_);
        done_cv_.wait(lk, [&] { return job.done.load() == n && job.helpers.load() == 0; });
        jobs_.erase(std::find(jobs_.begin(), jobs_.end(), &job));
    }
    if (job.error) std::rethrow_exception(job.error);
}

std::shared_ptr<ThreadPool::Task> ThreadPool::submit(std::function<void()> fn) {
    auto task = std::make_shared<Task>();
    task->fn = std::move(fn);
    if (workers_.empty()) {               // no workers: run it now
        run(task);
        return task;
    }
    {
        std::lock_guard<std::mutex> lk(m_);
        tasks_.push_back(task);
    }
    work_cv_.notify_one();
    done_cv_.notify_all();                // a waiter may pick it up
    return task;
}

void ThreadPool::run(const std::shared_ptr<Task>& task) {
    try {
        task->fn();
    } catch (...) {
        task->error = std::current_exception();
    }
    {
        std::lock_guard<std::mutex> lk(m_);
        task->done = true;
    }
    done_cv_.notify_all();
}

bool ThreadPool::help_once(bool tasks_too) {
    std::unique_lock<std::mutex> lk(m_);
    if (Job* job = claimable_job()) {
        ++job->helpers;
        lk.unlock();
        help(*job);
        lk.lock();
        if (--job->helpers == 0) done_cv_.notify_all();
        return true;
    }
    if (tasks_too && !tasks_.empty()) {
        std::shared_ptr<Task> task = std::move(tasks_.front());
        tasks_.pop_front();
        lk.unlock();
        run(task);
        return true;
    }
    return false;
}

void ThreadPool::wait(const std::shared_ptr<Task>& task) {
    while (!task->done.load()) {
        if (help_once(true)) continue;
        std::unique_lock<std::mutex> lk(m_);
        done_cv_.wait(lk, [&] { return task->done.load() || claimable_job() || !tasks_.empty(); });
    }
    if (task->error) std::rethrow_exception(task->error);
}

}  // namespace wqfm
