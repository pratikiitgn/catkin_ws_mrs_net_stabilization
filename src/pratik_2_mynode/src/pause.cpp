//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// pause.cpp
//
// Code generation for function 'pause'
//

// Include files
#include "pause.h"
#include "coder_posix_time.h"

// Variable Definitions
static unsigned char pauseState;

// Function Definitions
namespace coder {
void pause()
{
  coderTimespec b_timespec;
  if (pauseState == 0) {
    b_timespec.tv_sec = 3.0;
    b_timespec.tv_nsec = 0.0;
    coderTimeSleep(&b_timespec);
  }
}

} // namespace coder
void cpause_init()
{
  pauseState = 0U;
}

// End of code generation (pause.cpp)
