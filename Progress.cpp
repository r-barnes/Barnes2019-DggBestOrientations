#include "Progress.hpp"
#include <iostream>
#include <iomanip>

#ifdef _OPENMP
  ///Multi-threading - yay!
  #include <omp.h>
#else
  ///Macros used to disguise the fact that we do not have multithreading enabled.
  #define omp_get_thread_num()  0
  #define omp_get_num_threads() 1
#endif



///@brief Start/reset the progress bar.
///@param total_work  The amount of work to be completed, usually specified in cells.
ProgressBar::ProgressBar(uint32_t total_work0){
  total_work = total_work0;
  clearConsoleLine();
}



///Clear current line on console so a new progress bar can be written
void ProgressBar::clearConsoleLine() const {
  #ifdef NOPROGRESS
    return;
  #endif

  std::cerr<<"\r\033[2K"<<std::flush;
}



///@brief Update the visible progress bar, but only if enough work has been done.
///
///Define the global `NOPROGRESS` flag to prevent this from having an
///effect. Doing so may speed up the program's execution.
void ProgressBar::update(uint32_t work_done0){
  //Provide simple way of optimizing out progress updates
  #ifdef NOPROGRESS
    return;
  #endif

  work_done = work_done0; //Update the amount of work done

  //Quick return if this isn't the main thread
  if(omp_get_thread_num()!=0)
    return;

  if(work_done<next_update)
    return;

  const auto     num_threads = omp_get_num_threads();
  const double   elapsed     = timer.elapsed();       //Time since start
  const double   lap_time    = elapsed-prev_time;     //Time since last update
  const uint32_t lap_work    = work_done - prev_work; //Work completed in lap

  // std::cerr<<std::endl;
  // std::cerr<<"work_done = "<<work_done<<std::endl;
  // std::cerr<<"prev_work_done = "<<prev_work<<std::endl;
  // std::cerr<<"elapsed = "<<elapsed<<std::endl;
  // std::cerr<<"lap_time = "<<lap_time<<std::endl;
  // std::cerr<<"lap_work = "<<lap_work<<std::endl;

  prev_time = elapsed;    //Update time when last update occurred
  prev_work = work_done;

  const double works_per_sec = lap_work/lap_time;
  avg_time += 1/works_per_sec;
  avg_count++;

  // std::cerr<<"works_per_sec = "<<works_per_sec<<std::endl;

  //Update the next time at which we'll do the expensive update stuff
  next_update += works_per_sec*secs_between_updates;

  // std::cerr<<"next_update = "<<next_update<<std::endl;
  // std::cerr<<"works_per_sec = "<<works_per_sec<<std::endl;
  // std::cerr<<"avg_time  = "<<avg_count<<std::endl;
  // std::cerr<<"avg_count = "<<avg_count<<std::endl;
  const double time_remaining = avg_time/avg_count*(total_work/(double)num_threads-work_done);

  //Use a uint16_t because using a uint8_t will cause the result to print as a
  //character instead of a number
  double percent = 100*work_done*num_threads/total_work;

  //Handle overflows
  if(percent>100)
    percent=100;

  //Print an update string which looks like this:
  //  [================================================  ] (96% - 1.0s - 4 threads)
  std::cerr<<"\r\033[2K["
           <<std::string(percent/2, '=')<<std::string(50-percent/2, ' ')
           <<"] ("
           <<std::fixed<<std::setprecision(1)<<percent<<"% - "
           <<std::fixed<<std::setprecision(1)<<time_remaining
           <<"s - "
           <<omp_get_num_threads()<< " threads)"<<std::flush;
}

///Increment by one the work done and update the progress bar
ProgressBar& ProgressBar::operator++(){
  #ifdef NOPROGRESS
    return;
  #endif

  //Quick return if this isn't the main thread
  if(omp_get_thread_num()!=0)
    return *this;

  update(work_done+1);
  return *this;
}

///Stop the progress bar. Throws an exception if it wasn't started.
///@return The number of seconds the progress bar was running.
double ProgressBar::stop(){
  clearConsoleLine();

  return timer.elapsed();
}

///@return Return the time the progress bar ran for.
double ProgressBar::time_it_took(){
  return timer.elapsed();
}

uint32_t ProgressBar::cellsProcessed() const {
  return work_done;
}
