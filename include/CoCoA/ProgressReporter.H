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

#endif
