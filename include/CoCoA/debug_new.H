#ifndef debug_new_H
#define debug_new_H

//   Copyright (c)  2005  John Abbott,  Anna M. Bigatti

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


namespace debug_new
{

  /*-----------------------------------------------------------------*/
  /** \include debug_new.txt  */
  /*-----------------------------------------------------------------*/
  class PrintTrace
  {
  public:
    PrintTrace(bool activate = true);  // turns on printing of messages (if arg is true)
    ~PrintTrace();                     // restores printing to the previous state
  private:
    bool PreviousState;
  };

  void intercepted();                  // a good place to set debugger breakpoints

  void InterceptNew(unsigned long nth);     // call debug_new::intercepted function after nth new
  void InterceptDelete(unsigned long nth);  // call debug_new::intercepted function after nth delete

} // end of namespace debug_new

#endif
