//   Copyright (c)  2017  John Abbott,  Anna M. Bigatti

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

#include "CoCoA/SignalWatcher.H"
#include "CoCoA/assert.H"
#include "CoCoA/error.H"

#include <iostream>

namespace CoCoA
{

  namespace // anonymous
  {

    volatile sig_atomic_t SignalReceived = 0; // GLOBAL VARIABLE!!!
    
  } // end of anonymous namespace


  SignalWatcher::SignalWatcher(int sig, void FnPtr(int)):
      mySig(sig)
  {
#if defined (_WIN32) || defined (_WIN64)
    if (FnPtr == nullptr)
      myPrevSigactionPtr = signal(sig, SetSignalReceived); // default CoCoA signal handler
    else
      myPrevSigactionPtr = signal(sig, FnPtr);
#else
    struct sigaction sa;
    if (FnPtr == nullptr)
      sa.sa_handler = &SetSignalReceived; // default CoCoA signal handler
    else
      sa.sa_handler = FnPtr;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = SA_RESTART; // Restart functions if interrupted by handler

    myPrevSigactionPtr = new struct sigaction;
    if (sigaction(sig, &sa, myPrevSigactionPtr) != 0)
    {
      delete myPrevSigactionPtr;
      CoCoA_THROW_ERROR1("Unable to set signal handler");
    }
#endif
  }


  void SignalWatcher::myOutputSelf(std::ostream& out) const
  {
    if (!out)  return;  // short-cut for bad ostreams
    out << "SignalWatcher(sig=" << mySig;
    if (IamActive())
      out << ')';
    out << ", DEACTIVATED)";
  }


  void SignalWatcher::myDeactivate() noexcept
  {
    if (!IamActive())  return;
#if defined (_WIN32) || defined (_WIN64)
    signal(mySig, *myPrevSigactionPtr);
    myPrevSigactionPtr = nullptr;
#else
    sigaction(mySig, myPrevSigactionPtr, nullptr);
    delete myPrevSigactionPtr;
    myPrevSigactionPtr = nullptr;
#endif
  }


  SignalWatcher::~SignalWatcher()
  {
    myDeactivate();
  }


  std::ostream& operator<<(std::ostream& out, const SignalWatcher& SW)
  {
    if (!out) return out;  // short-cut for bad ostreams
    SW.myOutputSelf(out);
    return out;
  }


  int GetAndResetSignalReceived() noexcept
  {
    const int sig = SignalReceived;
    SignalReceived = 0;
    return sig;
  }

  void SetSignalReceived(int sig) noexcept
  {
    // SANITY CHECK on value of sig???
#if defined (_WIN32) || defined (_WIN64)
    signal(sig, SetSignalReceived);
#endif
    CoCoA_ASSERT(sig != 0);
    SignalReceived = sig;
  }


  //-------------------------------------------------------
  // Must define this because of virtual mem fns.
  InterruptedBySignal::~InterruptedBySignal()
  {}

  void InterruptedBySignal::myOutputSelf(std::ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams
    out << "CoCoA::InterruptedBySignal(signal=" << mySignal
        << ", context=\"" << myContext << "\")";
  }


} // end of namespace CoCoA
