#ifndef CoCoA_interrupt_H
#define CoCoA_interrupt_H

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

#include "CoCoA/exception.H"

//#include <string>
//using std::string;
#include <iosfwd>
//using std::ostream;

namespace CoCoA
{

  class InterruptReceived: public exception
  {
  private:
    static const char* const ourExceptionMesg;
  public:
    explicit InterruptReceived(const char* const context): exception(ourExceptionMesg, context) {};
//    explicit InterruptReceived(const std::string& context): exception(ourExceptionMesg, context) {};
    virtual ~InterruptReceived();
    void myOutputSelf(std::ostream& out) const override;
  };


  void CheckForInterrupt(const char* const context);  // User should call this fn "not too rarely"
//  void CheckForInterrupt(const std::string& context);

  
} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/interrupt.H,v 1.14 2022/02/18 14:11:52 abbott Exp $
// $Log: interrupt.H,v $
// Revision 1.14  2022/02/18 14:11:52  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.13  2021/10/30 17:27:36  abbott
// Summary: Used keyword override (redmine 1625)
//
// Revision 1.12  2021/10/07 12:11:43  abbott
// Summary: Added comment
//
// Revision 1.11  2021/01/07 14:48:33  abbott
// Summary: Corrected copyright
//
// Revision 1.10  2020/02/11 17:16:17  abbott
// Summary: All mem fns like myOutput now return void (redmine 969)
//
// Revision 1.9  2017/07/22 13:05:17  abbott
// Summary: Cleaning after creating new exceptions InterrptedBySignal & InterruptedByTimeout
//
// Revision 1.8  2017/07/21 13:21:22  abbott
// Summary: Split olf interrupt into two ==> new file SignalWatcher; refactored interrupt and CpuTimeLimit
//
// Revision 1.7  2017/07/14 19:22:18  abbott
// Summary: Minor improvements: better comments, removed cruft
//
// Revision 1.6  2017/07/14 13:54:40  abbott
// Summary: Added new class TimerInterruptReceived
//
// Revision 1.5  2017/07/08 19:05:51  abbott
// Summary: major revision to interrupt mechanism
//
// Revision 1.4  2016/11/18 18:10:08  abbott
// Summary: Exception type InterrupReceived now memorizes which signal triggered it; new fns TriggeredBySignal, SignalInterruptsCoCoA
//
// Revision 1.3  2015/06/26 14:58:03  abbott
// Summary: Now InterruptReceived derives from CoCoA::exception; CheckForInterrupt requires context string.
// Author: JAA
//
// Revision 1.2  2015/06/25 16:09:30  abbott
// Summary: InterruptReceived ctor now requires string arg; added printing.
// Author: JAA
//
// Revision 1.1  2015/05/20 14:44:58  abbott
// Summary: New fns for responding to interrupts
// Author: JAA
//
//

#endif
