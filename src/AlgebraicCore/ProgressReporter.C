//   Copyright (c)  2014  John Abbott and Anna M. Bigatti

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


#include "CoCoA/ProgressReporter.H"
#include "CoCoA/time.H"

#include <cmath>
using std::floor;
#include <iostream>
using std::cout;
using std::endl;

namespace CoCoA
{

  ProgressReporter::ProgressReporter(double interval):
      myIntervalCount(0),
      myCheckCount(1),
      myTotalCount(0),
      myLastCheckTime(CpuTime()),
      myTargetInterval(interval),
      myNextPrintTime(myLastCheckTime+interval),
      myLastPrintCount(0),
      myLastPrintTime(0)
  {}


  void ProgressReporter::myReset() noexcept
  {
    myTotalCount += myIntervalCount;
    myIntervalCount = 0;
    myCheckCount = 1;
    myLastCheckTime = CpuTime();
    myNextPrintTime = myLastCheckTime+myTargetInterval;
  }


  bool ProgressReporter::myIsTimeToPrint()
  {
    myTotalCount += myIntervalCount;
    myIntervalCount = 0;
    const double t = CpuTime();
    double ratio = (t-myLastCheckTime)/myTargetInterval;
    myLastCheckTime = t;
    if (ratio >= 0.125) { do { decrease125(myCheckCount); ratio *= 0.4;} while (ratio >= 0.125); }
    if (ratio < 0.03125 && myCheckCount < ourMaxCheckingInterval)  // UPB for myCheckCount to avoid overflow
    {
      increase125(myCheckCount);
      // Now adjust myIntervalCount, so next print will be at a multiple of myCheckCount
      const long tmp = myTotalCount%myCheckCount;
      if (tmp != 0)
      {
        myTotalCount -= tmp;
        myIntervalCount = tmp;
        return false;
      }
    }
    if (t < myNextPrintTime) return false;
    myNextPrintTime = t+myTargetInterval;
    return true;
  }

  // compute measured rate (iters/sec), and update myLastPrintCount & myLastPrintTime
  double ProgressReporter::myRate()
  {
    const double CurrTime = CpuTime();
    const double rate = (myTotalCount-myLastPrintCount)/(CurrTime-myLastPrintTime);
    myLastPrintCount = myTotalCount;
    myLastPrintTime = CurrTime;
    return rate;
  }

  void ProgressReporter::myPrintReport()
  {
    cout << "--> Progress count=" << myTotalCount << "   time=" << myLastCheckTime << " \trate=" << myRate() << " iter/sec" << endl;
  }

  void ProgressReporter::myPrintReport(long arg1)
  {
    cout << "--> Progress at "<< arg1 << "   count=" << myTotalCount << "   time=" << myLastCheckTime << " \trate=" << myRate() << " iter/sec" << endl;
  }

  void ProgressReporter::myPrintReport(long arg1, long arg2)
  {
    cout << "--> Progress at (" << arg1 << ", " << arg2 << ")   count=" << myTotalCount << "   time=" << myLastCheckTime << " \trate=" << myRate() << " iter/sec" << endl;
  }


  // Quick makeshift impl.
  std::ostream& operator<<(std::ostream& out, const ProgressReporter& PR)
  {
    if (!out) return out;  // short-cut for bad ostreams
    out << "ProgressReporter(intvl=" << PR.myTargetInterval << ")";
    return out;
  }


  void increase125(long& n)
  {
    long pwr10 = 1;
    while (n >= 10) { n /= 10; pwr10 *= 10; }
    if (n == 1) { n = 2*pwr10; return; }
    if (n == 2) { n = 5*pwr10; return; }
    n = 10*pwr10;
  }

  void decrease125(long& n)
  {
    if (n == 1) return;
    long pwr10 = 1;
    while (n >= 10) { n /= 10; pwr10 *= 10; }
    if (n == 1) { n = pwr10/2; return; }
    if (n == 2) { n = pwr10; return; }
    n = 2*pwr10;
  }


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/ProgressReporter.C,v 1.11 2024/06/09 20:10:20 abbott Exp $
// $Log: ProgressReporter.C,v $
// Revision 1.11  2024/06/09 20:10:20  abbott
// Summary: Renamed myPrepareForNewLoop to myReset (ProgressReporter & CpuTimeLimit); remine 1705
//
// Revision 1.10  2022/11/04 15:40:03  abbott
// Summary: Added myPrepareForNewLoop (redmine 1705)
//
// Revision 1.9  2022/02/18 14:11:56  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.8  2021/01/07 15:16:52  abbott
// Summary: Corrected copyright
//
// Revision 1.7  2019/09/16 17:23:48  abbott
// Summary: Added faster decrease for myCheckCount (if ratio is large)
//
// Revision 1.6  2018/06/25 15:12:58  abbott
// Summary: Fixed bug when UPB for interval count was reached
//
// Revision 1.5  2018/06/25 12:30:26  abbott
// Summary: Added overflow protection (when increasing interval)
//
// Revision 1.4  2016/11/11 14:15:33  abbott
// Summary: Added short-cut to operator<< when ostream is in bad state
//
// Revision 1.3  2014/10/13 18:05:43  abbott
// Summary: Added code to print out measured rate in each progress report
// Author: JAA
//
// Revision 1.2  2014/04/22 14:53:43  abbott
// Summary: Changed format of progress reports; now respects the chosen interval (rather than aiming for multiples of chosen interval)
// Author: JAA
//
// Revision 1.1  2014/04/22 13:26:01  abbott
// Summary: New class ProgressReporter
// Author: JAA
//
//
