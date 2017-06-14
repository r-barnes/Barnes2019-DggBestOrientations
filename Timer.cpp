#include "Timer.hpp"

Timer::Timer(){
  reset();
}

void Timer::reset(){
  start_time = clock::now();
}

double Timer::elapsed() const { 
  return std::chrono::duration_cast<second> (clock::now() - start_time).count(); 
}