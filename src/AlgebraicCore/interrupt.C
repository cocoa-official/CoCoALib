//   Copyright (c)  2015,2017  John Abbott,  Anna M. Bigatti

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


#include "CoCoA/interrupt.H"
#include "CoCoA/error.H"
#include "CoCoA/verbose.H"
#include "CoCoA/SignalWatcher.H"

#include <csignal>
//using std::signal;

#include <iostream>
using std::ostream;

namespace CoCoA
{

  const char* const InterruptReceived::ourExceptionMesg = "External interrupt";

  namespace // anonymous
  {
    // This is our signal "handler": just sets a "flag".
    // To be used in conjunction with CheckForInterrupt.
    void AnnounceInterruption(int sig)
    {
      VerboseLog VERBOSE("intr");
      VERBOSE(10) << std::endl;
      VERBOSE(10) << "--------------------------------" << std::endl;
      VERBOSE(10) << "-->>  CoCoALib interrupted  <<--" << std::endl;
      // padding in next line is just to get alignment right (assumes 0 <= sig <= 99)
      const char* const padding = ((sig < 10)?" ":"");
      VERBOSE(11) << "-->>  (by signal " << sig << ")" << padding << "        <<--" << std::endl;
      VERBOSE(10) << "--------------------------------" << std::endl;
    }

  } // end of anonymous namespace


  // ??? Should the first half of this fn be inline???
  void CheckForInterrupt(const char* const context)
  {
    const int sig = GetAndResetSignalReceived(); // resets the signal buffer
    if (sig == 0) return; // no signal received, so just return

    if (VerbosityLevel() >= 10) AnnounceInterruption(sig);
    ThrowException(InterruptedBySignal(sig, context));
  }

  // void CheckForInterrupt(const std::string& context)
  // {
  //   // Check FIRST for any signals, and then afterwards for a timeout.
  //   const int sig = GetAndResetSignalReceived(); // resets the signal buffer
  //   if (sig > 0)
  //   {
  //     if (VerbosityLevel() >= 10) AnnounceInterruption(sig);
  //     throw InterruptedBySignal(sig, context);
  //   }
  //   CheckForTimeout(context); // throws if timeout has occurred
  // }


  //-------------------------------------------------------
  // Must define this because of virtual mem fns.
  InterruptReceived::~InterruptReceived()
  {}


  void InterruptReceived::myOutputSelf(std::ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams
    out << "CoCoA::InterruptReceived(context=\"" << myContext << "\")";
  }

} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/interrupt.C,v 1.21 2023/10/05 17:58:37 abbott Exp $
// $Log: interrupt.C,v $
// Revision 1.21  2023/10/05 17:58:37  abbott
// Summary: Improved layout in printed mesg
//
// Revision 1.20  2023/07/26 11:56:37  abbott
// Summary: Minor tidying
//
// Revision 1.19  2022/10/25 17:19:20  abbott
// Summary: Added comment about making part of CheckForInterrupt inline
//
// Revision 1.18  2022/02/18 14:12:02  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.17  2021/10/07 12:12:30  abbott
// Summary: Cleaned CheckForInterrupt
//
// Revision 1.16  2020/06/19 14:59:00  abbott
// Summary: Calls ThrowException instead of throwing directly
//
// Revision 1.15  2020/02/11 16:56:42  abbott
// Summary: Corrected last update (see redmine 969)
//
// Revision 1.14  2018/05/25 09:24:47  abbott
// Summary: Major redesign of CpuTimeLimit (many consequences)
//
// Revision 1.13  2017/09/06 14:10:58  abbott
// Summary: Improved comment
//
// Revision 1.12  2017/07/22 16:11:50  abbott
// Summary: Added comment
//
// Revision 1.11  2017/07/22 13:05:17  abbott
// Summary: Cleaning after creating new exceptions InterrptedBySignal & InterruptedByTimeout
//
// Revision 1.10  2017/07/21 13:21:23  abbott
// Summary: Split olf interrupt into two ==> new file SignalWatcher; refactored interrupt and CpuTimeLimit
//
// Revision 1.9  2017/07/14 13:52:56  abbott
// Summary: Added new TimerInterruptReceived
//
// Revision 1.8  2017/07/08 19:05:51  abbott
// Summary: major revision to interrupt mechanism
//
// Revision 1.7  2017/03/01 17:11:05  bigatti
// -- modified message (again): now says CoCoALib
//
// Revision 1.6  2017/03/01 17:09:22  bigatti
// -- slightly modified printout (to make it "red")
//
// Revision 1.5  2016/11/18 18:12:09  abbott
// Summary: InterruptReceived now stores signal; added new fns TriggeredBySignal, SignalInterruptsCoCoA
//
// Revision 1.4  2016/11/11 14:15:34  abbott
// Summary: Added short-cut to operator<< when ostream is in bad state
//
// Revision 1.3  2015/06/26 14:58:03  abbott
// Summary: Now InterruptReceived derives from CoCoA::exception; CheckForInterrupt requires context string.
// Author: JAA
//
// Revision 1.2  2015/06/25 16:09:31  abbott
// Summary: InterruptReceived ctor now requires string arg; added printing.
// Author: JAA
//
// Revision 1.1  2015/05/20 14:44:58  abbott
// Summary: New fns for responding to interrupts
// Author: JAA
//
//
