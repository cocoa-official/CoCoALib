//   Copyright (c)  2005,2013  John Abbott and Anna M. Bigatti

//   This file is part of the source of CoCoALib, the CoCoA Library.
//
//   CoCoALib is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   CoCoALib is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with CoCoALib.  If not, see <http://www.gnu.org/licenses/>.


#include "CoCoA/time.H"

#include <chrono>

namespace CoCoA
{

  static const auto StartTime = std::chrono::steady_clock::now();

  double ElapsedTime() noexcept
  {
    auto CurrTime = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed = CurrTime-StartTime;
    return elapsed.count();
  }
}

#if defined(__MINGW32__) || defined(_WINDOWS)
#include <ctime>
// this works for linux too
double SecondsFromClock() noexcept
{
  const double floatSpan = clock();
  const double seconds = (floatSpan / CLOCKS_PER_SEC);
  return seconds;
}


// Non Unix-like environment, so use "lobotomized" time fns.
namespace CoCoA
{

  double CpuTime() noexcept
  {
    return SecondsFromClock();
  }

  double RealTime() noexcept
  {
    return 0.0;
  }


  void DateTime(long& date, long& time) noexcept
  {
    date = 0;
    time = 0;
  }
  
} // end of namespace CoCoA

#else

// Good news: we're in a unix-like environment.
// These are the normal defns.

#include <ctime>
using std::time;
#include <sys/time.h>
#include <sys/resource.h>
//not needed???? #include <unistd.h>

namespace CoCoA
{

  double CpuTime() noexcept
  {
    static struct rusage used;

    getrusage(RUSAGE_SELF, &used);
    // Use volatile to force truncation to IEEE doubles in x86 processors;
    // otherwise it is possible to get negative time differences!
    volatile double seconds = used.ru_utime.tv_sec + used.ru_stime.tv_sec +
      (used.ru_utime.tv_usec + used.ru_stime.tv_usec) / 1.0e6;
    return seconds;
  }

  double RealTime() noexcept
  {
    static struct timeval tv;

    gettimeofday(&tv, nullptr);
    // Use volatile to force truncation to IEEE doubles in x86 processors;
    // otherwise it is possible to get negative time differences!
    volatile double seconds = tv.tv_sec + tv.tv_usec/1.0e6;
    return seconds;
  }

  void DateTime(long& date, long& time) noexcept
  {
    std::time_t SecondsSinceEpoch;
    std::time(&SecondsSinceEpoch); // gets current time - Async-signal-safe
//   const string date = std::ctime_r(&SecondsSinceEpoch);
//   cout << "Time & date: " << date << endl;

    std::tm CurrentTime;
    localtime_r(&SecondsSinceEpoch, &CurrentTime); // thread-safe
    const int year = CurrentTime.tm_year + 1900;
    const int month = CurrentTime.tm_mon + 1;
    const int day = CurrentTime.tm_mday;
    const int hour = CurrentTime.tm_hour;
    const int minute = CurrentTime.tm_min;
    const int second = CurrentTime.tm_sec;

    date = 10000*year + 100*month + day;
    time = 10000*hour + 100*minute + second;
  }
  

} // end of namespace CoCoA


#endif
