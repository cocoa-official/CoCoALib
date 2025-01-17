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

#endif
