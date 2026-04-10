#ifndef POLYTOPE_TIMINGUTILITIES_HH
#define POLYTOPE_TIMINGUTILITIES_HH
//------------------------------------------------------------------------------
// timingUtilities
//
// A set of inline helper methods to encapsulate how we do timing.
//
// JMO:  Tue Dec  9 10:31:14 PST 2008
//------------------------------------------------------------------------------

#include <chrono>

namespace polytope {
//------------------------------------------------------------------------------
// Get the current clock time.
//------------------------------------------------------------------------------
struct Timing {
  typedef std::chrono::time_point<std::chrono::high_resolution_clock> Time;
  typedef std::chrono::microseconds duration;
  static Time currentTime() { return std::chrono::high_resolution_clock::now(); }
  static double difference(const Time& t1, const Time& t2) { return (std::chrono::duration_cast<duration>(t2 - t1)).count(); }
};
}

#endif
