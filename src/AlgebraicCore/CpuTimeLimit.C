//   Copyright (c)  2017,2018  John Abbott,  Anna M. Bigatti

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


#include "CoCoA/CpuTimeLimit.H"
#include "CoCoA/assert.H"
#include "CoCoA/error.H"
#include "CoCoA/time.H"

#include <iostream>
using std::ostream;

namespace CoCoA
{

  //-------------------------------------------------------
  // Must define this because of virtual mem fns.
  TimeoutException::~TimeoutException()
  {}


  void TimeoutException::myOutputSelf(std::ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams
    out << "CoCoA::TimeoutException(context=\"" << myContext << "\")";
  }
    
  //------------------------------------------------------------------

  namespace // anonymous
  {
    int QuantifyVariability(IterationVariability v)
    {
      if (v == IterationVariability::low) return 1;
      if (v == IterationVariability::medium) return 4;
      return 16; // v == IterationVariability::high
    }
    
  } // end of namespace anonymous
  
  //------------------------------------------------------------------
  
  CpuTimeLimit::CpuTimeLimit(double interval, IterationVariability v):
      myCountdown(1),
      myCheckingInterval(1),
      myTotalCount(0),
      myRefPt1Count(0),
      myRefPt2Count(0),
      myRefPt1Time(ElapsedTime()),
      myRefPt2Time(myRefPt1Time),
      myTriggerTime(myRefPt1Time+1.002*interval),
      myTriggerTimeCPU(CpuTime()+interval),
      myExtraTime(interval/32),
///???      myTriggerTimePlusEpsilon(myTriggerTime+interval/16+0.0625),
      myVariability(QuantifyVariability(v)),
      myScaleFactor(1)
  {
    static const char* const FnName = "CpuTimeLimit ctor";
    if (interval < 0) CoCoA_THROW_ERROR(ERR::ReqNonNegative, FnName);
    if (interval > 1000000) CoCoA_THROW_ERROR(ERR::ArgTooBig, FnName);
///    if (myVariability < 1) CoCoA_THROW_ERROR(ERR::ReqPositive, FnName);
///    if (myVariability > 256) CoCoA_THROW_ERROR(ERR::ArgTooBig, FnName);
  }

  
  void CpuTimeLimit::myReset(IterationVariability v) const
  {
///    static const char* const FnName = "myReset";
    if (IamUnlimited()) return; // do nothing
///    long v = AsSignedLong(variability);
///    if (v < 1) CoCoA_THROW_ERROR(ERR::ReqPositive, FnName);
///    if (v > 256) CoCoA_THROW_ERROR(ERR::ArgTooBig, FnName);
    myCheckingInterval = 1;
    myCountdown = 1;
    myRefPt1Count = myTotalCount;
    myRefPt2Count = myRefPt1Count;
    const double now = ElapsedTime();
    myRefPt1Time = now;
    myRefPt2Time = now;
    myVariability = QuantifyVariability(v);
  }

  
  bool CpuTimeLimit::IamTimedOut_ProperCheck() const noexcept
  {
    if (IamUnlimited()) return false;  // Should (almost) never happen!
//    --myCountdown;
//    if (myCountdown > 0) return false; // Only really check when myCountDown reaches 0
    CoCoA_ASSERT(myCountdown == 0);
    myTotalCount += myCheckingInterval;
    const double now = ElapsedTime();
    if (now >= myTriggerTime)
    {
//??      if (JustUsingElapsedTime) return true;
//      std::clog<<'*';
//      std::clog<<"CPU CHECK  elapsed+"<<now<<"   ElapsedTrigger = "<<myTriggerTime<<std::endl;
      const double CpuNow = CpuTime();
//      std::clog << "CpuTrigger: " << myTriggerTimeCPU << "   actual CPU: "<<CpuNow<<std::endl;
      if (CpuNow >= myTriggerTimeCPU) return true;
      const double TimeToDeadline = myTriggerTimeCPU - CpuNow;
      ++myScaleFactor;
//std::clog<<"ScaleFactor="<<myScaleFactor<<std::endl;
      myTriggerTime = now + TimeToDeadline + myScaleFactor*myExtraTime;
//      std::clog<<"New (elapsed) trigger time: " << myTriggerTime << std::endl;
    }
/// std::clog<<'.';
    if (myTotalCount-myRefPt2Count > 15) { myRefPt1Count = myRefPt2Count; myRefPt1Time = myRefPt2Time; myRefPt2Count = myTotalCount; myRefPt2Time = now; }
    // Compute ave time per count
    const double AveTime = (now-myRefPt1Time)/(myTotalCount-myRefPt1Count);
    const double TimeToDeadline = myTriggerTime - now;
    double EstNextTimeInterval = myVariability*myCheckingInterval*AveTime;
///std::clog<<"myRefPt1Count="<<myRefPt1Count<<"   myRefPt1Time="<<myRefPt1Time<<"   myRefPt2Count="<<myRefPt2Count<<"   myRefPt2Time="<<myRefPt2Time<<std::endl;
///std::clog<<"count diff="<<myTotalCount-myRefPt1Count<<"   AveTime="<<AveTime<<"  rem="<<TimeToDeadline<<"   est="<<EstNextTimeInterval<<"   CheckingInterval="<<myCheckingInterval<<std::endl;
    if (EstNextTimeInterval < TimeToDeadline/4) { if (myCheckingInterval < 32) { myCheckingInterval *= 2; } }
    else if (EstNextTimeInterval > TimeToDeadline)
    {
      while (/*EstNextTimeInterval > 0.1 &&*/ EstNextTimeInterval > TimeToDeadline && myCheckingInterval > 1)
      { EstNextTimeInterval /= 2; myCheckingInterval /= 2; }
    }
///std::clog<<"CHECK AGAIN AFTER "<<myCheckingInterval<<std::endl;
    myCountdown = myCheckingInterval;
    return false;
  }


  // Quick makeshift impl.
  std::ostream& operator<<(std::ostream& out, const CpuTimeLimit& TimeLimit)
  {
    if (!out) return out;  // short-cut for bad ostreams
    
    if (IsUnlimited(TimeLimit)) return out << "CpuTimeLimit(UNLIMITED)";

    out << "CpuTimeLimit(TriggerTime=" << TimeLimit.myTriggerTimeCPU
        << ", CurrTime=" << CpuTime()
        << ",  Countdown=" << TimeLimit.myCountdown
        << ", CheckingInterval=" << TimeLimit.myCheckingInterval << ")";
    return out;
  }


  // Special ctor for "unlimited" CpuTimeLimit object;
  // called only by NoCpuTimeLimit (below).
  CpuTimeLimit::CpuTimeLimit(NO_TIME_LIMIT) noexcept:
      myCountdown(-1),
      myCheckingInterval(-1), // negative myInterval marks out the "unlimited" object
      myTotalCount(-1),
      myRefPt1Count(-1),
      myRefPt2Count(-1),
      myRefPt1Time(-1.0),
      myRefPt2Time(-1.0),
      myTriggerTime(-1.0),
      myTriggerTimeCPU(-1.0),
///???      myTriggerTimePlusEpsilon(-1.0),
      myVariability(-1)
    {}
  

  const CpuTimeLimit& NoCpuTimeLimit() noexcept
  {
    static const CpuTimeLimit SingleCopy(CpuTimeLimit::NO_TIME_LIMIT::marker);
    return SingleCopy;
  }


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/CpuTimeLimit.C,v 1.26 2024/07/02 15:29:22 bigatti Exp $
// $Log: CpuTimeLimit.C,v $
// Revision 1.26  2024/07/02 15:29:22  bigatti
// Summary: Changed error codes: LogZero into ReqNonZero
//
// and Not... into ReqNonNegative, ReqNonNegativeGrading, ReqPositive, ReqPositiveGrading
// Specialized some BadArg with explicit codes
// Tidied up some CoCoA_THROW_ERROR second argument (string)
//
// Revision 1.25  2024/06/09 20:10:20  abbott
// Summary: Renamed myPrepareForNewLoop to myReset (ProgressReporter & CpuTimeLimit); remine 1705
//
// Revision 1.24  2023/01/01 11:30:49  abbott
// Summary: IamTimedOut is now public; minor refactoring
//
// Revision 1.23  2022/02/18 14:11:53  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.22  2021/03/03 22:09:32  abbott
// Summary: New enum class (redmine 894)
//
// Revision 1.21  2021/02/18 16:43:27  abbott
// Summary: Revised impl of CpuTimeLimit (redmine 1558)
//
// Revision 1.20  2021/02/10 19:39:59  abbott
// Summary: Added noexcept (sometimes instead of throw()) -- redmine 1572
//
// Revision 1.19  2021/01/12 13:27:25  abbott
// Summary: Added "varibility" to ctor; added myScaleFactor to handle when CPU is slower than steady clock
//
// Revision 1.18  2021/01/08 17:41:11  abbott
// Summary: Impl now more flexible (see redmine 1558)
//
// Revision 1.17  2020/06/17 15:49:22  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.16  2020/02/11 16:56:40  abbott
// Summary: Corrected last update (see redmine 969)
//
// Revision 1.15  2020/01/18 21:33:18  abbott
// Summary: Added two checks for being unlimited
//
// Revision 1.14  2019/12/21 16:40:16  abbott
// Summary: Added "variability"; revised myPrepareForNewLoop
//
// Revision 1.13  2019/12/20 15:51:38  abbott
// Summary: Major revision to CpuTimeLimit
//
// Revision 1.12  2019/10/29 11:35:46  abbott
// Summary: Replaced using namespace std by th specific using fn-name directives.
//
// Revision 1.11  2018/06/27 10:20:16  abbott
// Summary: Updated
//
// Revision 1.10  2018/06/27 09:37:57  abbott
// Summary: More detailed printout (more helpful for debugging)
//
// Revision 1.9  2018/06/25 12:31:49  abbott
// Summary: Added overflow protection when increasing interval length
//
// Revision 1.8  2018/05/25 09:24:46  abbott
// Summary: Major redesign of CpuTimeLimit (many consequences)
//
// Revision 1.7  2017/09/06 14:08:50  abbott
// Summary: Changed name to TimeoutException
//
// Revision 1.6  2017/07/23 15:32:32  abbott
// Summary: Fixed STUPID bug in myDeactivate
//
// Revision 1.5  2017/07/22 13:03:02  abbott
// Summary: Added new exception InterruptdByTimeout; changed rtn type of myOutputSelf
//
// Revision 1.4  2017/07/21 15:06:10  abbott
// Summary: Major revision -- no longer needs BOOST
//
// Revision 1.3  2017/07/21 13:21:22  abbott
// Summary: Split olf interrupt into two ==> new file SignalWatcher; refactored interrupt and CpuTimeLimit
//
// Revision 1.2  2017/07/15 15:44:44  abbott
// Summary: Corrected error ID name (to ArgTooBig)
//
// Revision 1.1  2017/07/15 15:17:48  abbott
// Summary: Added CpuTimeLimit
//
//
