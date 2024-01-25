/////////////////////////////////////////////////////////////////////////////
/// @file ThreadPool.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2017 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef __CAPD_THREADING_THREADPOOL_H__
#define __CAPD_THREADING_THREADPOOL_H__

#include <vector>
#include <list>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <exception>

namespace capd{
namespace threading{
/// @addtogroup threading 
/// @{

template<class F>
struct LambdaTask;

/**
 * This class is an abstract Task to be executed in a thread pool.
 * An inheriting class should override the method run which will be exectued in the task.
 *
 * An instance of processed task should exists until the thread pool to which the task is submitted
 * completes the executon of the task or until the pool has been interrupted.
 * Otherwise the behaviour is undefined.
 */
class Task {
public:
  typedef std::unique_lock<std::mutex> Lock;

  Task(const Task&) = delete;
  Task() : ready(false) {}
  virtual ~Task(){}
  virtual void run(unsigned threadId) = 0;

  /// wait for task completion
  void join(){
    Lock lock(m);
    while(!this->ready)
      completed.wait(lock);
  }

  /// This method can be invoked by the user but it is recomended to run tasks through a ThreadPool instance.
  void execute(unsigned threadId) {
    Lock lock(m);
    if(!this->ready)
    {
      this->run(threadId);
      this->ready = true;
    }
    completed.notify_all();
  }

  template<class F>
  static LambdaTask<F>* createFromLambda(F& f);
private:

  std::mutex m;
  std::condition_variable completed;
  bool ready;
};

template<class F>
struct LambdaTask : public capd::threading::Task{
  LambdaTask(F& _f) : f(_f){}
  void run(unsigned threadId) { this->f(threadId); }
  F& f;
};

template<class F>
LambdaTask<F>* Task::createFromLambda(F& f){
  return new LambdaTask<F>(f);
}

// #####################################################################

/**
 * This class realizes a simple thread pool of fixed size. All the threads are keeped alive until
 * the instance of ThreadPool exists or it has been marked as interrupted.
 */
class ThreadPool {
public:
  typedef std::unique_lock<std::mutex> UniqueLock;
  typedef std::lock_guard<std::mutex> LockGuard;

  /// Contructor takes number of threads in the pool.
  ThreadPool() = delete;
  ThreadPool(const ThreadPool&) = delete;
  ThreadPool(ThreadPool&) = delete;
  ThreadPool(ThreadPool&&) = delete;

  explicit ThreadPool(unsigned noThreads) : activeWorkers(noThreads)
  {
    for(unsigned i=0;i<noThreads;++i)
      workers.emplace_back(run,this,i);
  }

  /// Task submission. It should be a valid pointer to an existing instance od Task.
  void process(Task* task){
    LockGuard lock(m);
    if(!this->interrupted){
      tasks.push_back(task);
      nonempty.notify_one();
    }
    else
      throw std::runtime_error("A Task sent for processing to an interrupted ThreadPool.");
  }

  /// This method marks a ThreadPool for interruption.
  /// The tasks already running may continue for a while (even they can complete).
  /// No new tasks will be started after call to interrupt.
  /// Submissions of new tasks to an interrupted ThreadPool leads in throwing an exception.
  void interrupt(){
    {
      LockGuard lock(m);
      interrupted = true;
      tasks.clear();
    }
    nonempty.notify_all();
  }

  ~ThreadPool(){
    this->interrupt();
    for(auto& t : workers)
      t.join();
  }

  void waitUntilIdle(){
    UniqueLock lock(m);
    while(activeWorkers!=0 or !tasks.empty())
      idle.wait(lock);
  }

  unsigned poolSize() const { return workers.size(); }
private:

  Task* next(){
    UniqueLock lock(m);
    while(tasks.empty() and !this->interrupted){
      idle.notify_one();
      nonempty.wait(lock);
    }

    if(this->interrupted) {
      return nullptr;
    }
    Task* task = tasks.front();
    tasks.pop_front();
    return task;
  }

  static void run(ThreadPool* pool, unsigned threadId){
    while(true){
      {
        LockGuard lock(pool->m);
        pool->activeWorkers--;
      }
      Task* task = pool->next();
      if(task!=nullptr){
        {
          LockGuard lock(pool->m);
          pool->activeWorkers++;
        }
        task->execute(threadId);
      } else {
        pool->nonempty.notify_one();
        pool->idle.notify_one();
        return;
      }
    }
  }

  std::mutex m;
  std::vector<std::thread> workers;
  std::list<Task*> tasks;
  volatile bool interrupted = false;
  std::condition_variable nonempty;
  std::condition_variable idle;
  int activeWorkers;
};

/// @}
}} // namespace capd::threading

#endif
