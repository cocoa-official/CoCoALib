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

#endif
