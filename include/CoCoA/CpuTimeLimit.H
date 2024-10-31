#ifndef CoCoA_CpuTimeLimit_H
#define CoCoA_CpuTimeLimit_H

//   Copyright (c)  2017,2018,2021  John Abbott, Anna Bigatti

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


#include "CoCoA/error.H"
#include "CoCoA/MachineInt.H"

#include <iosfwd>
//using std::ostream;

namespace CoCoA
{

  class TimeoutException: public ErrorInfo
  {
  public:
    TimeoutException(const char* const context): ErrorInfo(ERR::TimedOut, context) {}
//    TimeoutException(const std::string& context): ErrorInfo(ERR::TimedOut, context) {}
    virtual ~TimeoutException();
    virtual void myOutputSelf(std::ostream& out) const override;
  private: // same data members as class ErrorInfo
  };

  //------------------------------------------------------------------

  enum class IterationVariability { low, medium, high };

  class CpuTimeLimit
  {
  public:
    explicit CpuTimeLimit(double interval, IterationVariability v = IterationVariability::medium); // NOT noexcept (arg check)
    bool IamTimedOut() const noexcept { if (--myCountdown != 0) return false; /*else*/ return IamTimedOut_ProperCheck(); /*resets myCountdown*/ }
    void operator()(const char* const context) const { if (IamTimedOut()) ThrowException(TimeoutException(context)); }
    void myReset(IterationVariability v = IterationVariability::medium) const;
    bool IamUnlimited() const { return (myCheckingInterval < 0); };

  private: // data members
    mutable int myCountdown;
    mutable int myCheckingInterval; // count between true timeout checks
    mutable long myTotalCount;
    mutable long myRefPt1Count;
    mutable long myRefPt2Count;
    mutable double myRefPt1Time;
    mutable double myRefPt2Time;
    mutable double myTriggerTime;
    double myTriggerTimeCPU;
    double myExtraTime;
    mutable int myVariability;
    mutable int myScaleFactor;
    friend std::ostream& operator<<(std::ostream& out, const CpuTimeLimit& );
  private: // impl details
    bool IamTimedOut_ProperCheck() const noexcept;

    enum class NO_TIME_LIMIT { marker };
    explicit CpuTimeLimit(NO_TIME_LIMIT) noexcept;
    friend const CpuTimeLimit& NoCpuTimeLimit() noexcept;
  };

  std::ostream& operator<<(std::ostream& out, const CpuTimeLimit& TimeLimit);

  const CpuTimeLimit& NoCpuTimeLimit() noexcept;

  inline bool IsUnlimited(const CpuTimeLimit& TimeLimit) noexcept
  { return TimeLimit.IamUnlimited(); }

} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/CpuTimeLimit.H,v 1.22 2024/06/09 20:10:20 abbott Exp $
// $Log: CpuTimeLimit.H,v $
// Revision 1.22  2024/06/09 20:10:20  abbott
// Summary: Renamed myPrepareForNewLoop to myReset (ProgressReporter & CpuTimeLimit); remine 1705
//
// Revision 1.21  2023/02/22 20:16:30  abbott
// Summary: Corrected stupid bug (now the code works again)
//
// Revision 1.20  2023/01/01 11:30:49  abbott
// Summary: IamTimedOut is now public; minor refactoring
//
// Revision 1.19  2022/02/18 14:11:45  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.18  2021/10/30 16:54:44  abbott
// Summary: Used keyword override (redmine 1625)
//
// Revision 1.17  2021/03/03 22:09:32  abbott
// Summary: New enum class (redmine 894)
//
// Revision 1.16  2021/02/18 16:43:27  abbott
// Summary: Revised impl of CpuTimeLimit (redmine 1558)
//
// Revision 1.15  2021/02/10 19:32:00  abbott
// Summary: Added noexcept (redmine 1572)
//
// Revision 1.14  2021/01/12 13:27:25  abbott
// Summary: Added "variability" to ctor; added myScaleFactor to handle when CPU is slower than steady clock
//
// Revision 1.13  2021/01/08 17:41:11  abbott
// Summary: Impl now more flexible (see redmine 1558)
//
// Revision 1.12  2020/06/19 14:58:21  abbott
// Summary: Calls ThrowException instead of directly throwing
//
// Revision 1.11  2020/02/11 17:16:17  abbott
// Summary: All mem fns like myOutput now return void (redmine 969)
//
// Revision 1.10  2020/01/18 21:32:44  abbott
// Summary: Maybe myCountdown should be unsigned?
//
// Revision 1.9  2019/12/21 16:40:16  abbott
// Summary: Added "variability"; revised myPrepareForNewLoop
//
// Revision 1.8  2019/12/20 15:51:38  abbott
// Summary: Major revision to CpuTimeLimit
//
// Revision 1.7  2018/06/27 10:20:16  abbott
// Summary: Updated
//
// Revision 1.6  2018/05/25 09:24:46  abbott
// Summary: Major redesign of CpuTimeLimit (many consequences)
//
// Revision 1.5  2017/09/06 14:02:16  abbott
// Summary: Added new ERR::TimedOut; TimeoutException derives from ErrorInfo
//
// Revision 1.4  2017/07/22 12:56:34  abbott
// Summary: Updated exception type
//
// Revision 1.3  2017/07/21 15:06:09  abbott
// Summary: Major revision -- no longer needs BOOST
//
// Revision 1.2  2017/07/21 13:21:21  abbott
// Summary: Split olf interrupt into two ==> new file SignalWatcher; refactored interrupt and CpuTimeLimit
//
// Revision 1.1  2017/07/15 15:18:43  abbott
// Summary: Added CpuTimeLimit
//
//

#endif
