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
