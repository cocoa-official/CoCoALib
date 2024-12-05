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
