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
