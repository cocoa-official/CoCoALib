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


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/time.C,v 1.14 2022/02/18 14:12:02 abbott Exp $
// $Log: time.C,v $
// Revision 1.14  2022/02/18 14:12:02  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.13  2021/02/10 19:40:01  abbott
// Summary: Added noexcept (sometimes instead of throw()) -- redmine 1572
//
// Revision 1.12  2021/01/07 15:23:38  abbott
// Summary: Corrected copyright
//
// Revision 1.11  2020/12/22 15:04:50  abbott
// Summary: Added ElapsedTime (redmine 1553)
//
// Revision 1.10  2019/03/19 11:07:08  abbott
// Summary: Replaced 0 by nullptr where appropriate
//
// Revision 1.9  2013/05/20 15:57:54  abbott
// Corrected copyright date.
//
// Revision 1.8  2013/05/15 09:33:11  bigatti
// -- added tentative CpuTime for windows
//
// Revision 1.7  2012/04/02 15:29:09  abbott
// Added extra CPP check, so that lobotomized versions are used on Windoze.
//
// Revision 1.6  2011/10/04 13:03:02  bigatti
// -- new logo for gui
//
// Revision 1.5  2011/07/15 16:54:24  abbott
// Added some missing consts.
//
// Revision 1.4  2011/07/06 15:49:40  bigatti
// -- added DateTime
//
// Revision 1.3  2011/06/04 08:43:43  abbott
// Added MINGW workaround.
//
// Revision 1.2  2007/10/30 17:14:06  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.1  2005/11/15 12:24:35  cocoa
// New timing functions.
//
