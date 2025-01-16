#ifndef CoCoA_SignalWatcher_H
#define CoCoA_SignalWatcher_H

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


#include "CoCoA/interrupt.H"

#include <iosfwd>
//using std::ostream;
#include <csignal>
//using sigaction; (perhaps not in namespace std?)

namespace CoCoA
{

  class InterruptedBySignal: public InterruptReceived
  {
  public:
    InterruptedBySignal(int sig, const char* const context): InterruptReceived(context), mySignal(sig) {}
    virtual ~InterruptedBySignal();
    void myOutputSelf(std::ostream& out) const override;
  private: // data members (in addition to those of CoCoA::InterruptReceived)
    const int mySignal; // signal which caused the interrupt
    friend int TriggeringSignal(const InterruptedBySignal& intr) noexcept;
  };

  inline int TriggeringSignal(const InterruptedBySignal& intr) noexcept
  { return intr.mySignal; }


  // RAII class: ctor sets new sigaction, dtor restores previous one.
  // NB myDeactivate restores previous sigaction.
  class SignalWatcher
  {
  public:
    explicit SignalWatcher(int sig, void FnPtr(int) = nullptr);
    ~SignalWatcher();
    void myOutputSelf(std::ostream& out) const;
    friend void deactivate(SignalWatcher& SW) noexcept { SW.myDeactivate(); }
  private:
    bool IamActive() const noexcept { return (myPrevSigactionPtr != nullptr); }
    void myDeactivate() noexcept;
  private: // data members
    int mySig;
#if defined (_WIN32) || defined (_WIN64)
    void (*myPrevSigactionPtr)(int);
#else
    struct sigaction* myPrevSigactionPtr;
#endif
  };

  void deactivate(SignalWatcher& SW) noexcept;
  std::ostream& operator<<(std::ostream& out, const SignalWatcher& w);

  int GetAndResetSignalReceived() noexcept;
  void SetSignalReceived(int sig) noexcept; // used by interpreter & C5 (gui)

} // end of namespace CoCoA

#endif
