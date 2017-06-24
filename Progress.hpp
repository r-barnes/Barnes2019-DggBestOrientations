#ifndef _progress_hpp_
#define _progress_hpp_

#include <cstdint>
#include "Timer.hpp"

///@brief Manages a console-based progress bar to keep the user entertained.
///
///Defining the global `NOPROGRESS` will
///disable all progress operations, potentially speeding up a program. The look
///of the progress bar is shown in ProgressBar.hpp.
class ProgressBar{
 private:
  int secs_between_updates = 5;
  uint32_t total_work;      ///< Total work to be accomplished
  uint32_t prev_work   = 0; ///< Work done at last update
  double   prev_time   = 0; ///< Time of last update
  double   avg_time    = 0; ///< Time taken to complete a given amount of work
  uint32_t avg_count   = 0; ///< Number of times stored for averaging
  uint32_t next_update = 1; ///< Next point to update the visible progress bar
  uint32_t work_done   = 0;
  Timer    timer;         ///< Used for generating ETA

  ///Clear current line on console so a new progress bar can be written
  void clearConsoleLine() const;

 public:
  ///@brief Start the progress bar.
  ///@param total_work  The amount of work to be completed, usually specified in cells.
  ProgressBar(uint32_t total_work0);

  ///@brief Update the visible progress bar, but only if enough work has been done.
  ///
  ///Define the global `NOPROGRESS` flag to prevent this from having an
  ///effect. Doing so may speed up the program's execution.
  void update(uint32_t work_done0);

  ///Increment by one the work done and update the progress bar
  ProgressBar& operator++();

  ///Stop the progress bar. Throws an exception if it wasn't started.
  ///@return The number of seconds the progress bar was running.
  double stop();

  ///@return Return the time the progress bar ran for.
  double time_it_took();

  uint32_t cellsProcessed() const;
};

#endif