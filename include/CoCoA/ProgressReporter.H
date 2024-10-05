#ifndef CoCoA_ProgressReporter_H
#define CoCoA_ProgressReporter_H

//   Copyright (c)  2014  John Abbott,  Anna M. Bigatti

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


#include <iosfwd>
//using std::ostream;

namespace CoCoA
{

  class ProgressReporter
  {
  public:
    explicit ProgressReporter(double interval);
    void operator()() { if (++myIntervalCount == myCheckCount && myIsTimeToPrint()) myPrintReport(); }
    void operator()(long arg1) { if (++myIntervalCount == myCheckCount && myIsTimeToPrint()) myPrintReport(arg1); }
    void operator()(long arg1, long arg2) { if (++myIntervalCount == myCheckCount && myIsTimeToPrint()) myPrintReport(arg1, arg2); }
    void myReset() noexcept;

  private: // data members
    long myIntervalCount;
    long myCheckCount;
    long myTotalCount;
    double myLastCheckTime;
    double myTargetInterval;
    double myNextPrintTime;
    long myLastPrintCount;  // used only for computing the rate
    double myLastPrintTime; // (ditto)
    friend std::ostream& operator<<(std::ostream& out, const ProgressReporter& PR);
  private: // impl details
    double myRate(); // called only by the printing fns; auto-updates myLastPrintCount & myLastPrintTime
    void myPrintReport();
    void myPrintReport(long arg1);
    void myPrintReport(long arg1, long arg2);
    bool myIsTimeToPrint();
    static constexpr long ourMaxCheckingInterval = 1000000;
  };

  std::ostream& operator<<(std::ostream& out, const ProgressReporter& PR);

  void increase125(long& n);
  void decrease125(long& n);


} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/ProgressReporter.H,v 1.6 2024/06/09 20:10:20 abbott Exp $
// $Log: ProgressReporter.H,v $
// Revision 1.6  2024/06/09 20:10:20  abbott
// Summary: Renamed myPrepareForNewLoop to myReset (ProgressReporter & CpuTimeLimit); remine 1705
//
// Revision 1.5  2022/11/04 15:40:03  abbott
// Summary: Added myPrepareForNewLoop (redmine 1705)
//
// Revision 1.4  2022/02/18 14:11:47  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.3  2018/06/25 15:32:22  abbott
// Summary: Corrected copyright line
//
// Revision 1.2  2014/10/13 18:05:43  abbott
// Summary: Added code to print out measured rate in each progress report
// Author: JAA
//
// Revision 1.1  2014/04/28 16:48:41  abbott
// Summary: Header file for ProgressReporter
// Author: JAA
//
//

#endif
