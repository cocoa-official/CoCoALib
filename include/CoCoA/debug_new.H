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



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/debug_new.H,v 1.5 2022/02/18 14:11:51 abbott Exp $
// $Log: debug_new.H,v $
// Revision 1.5  2022/02/18 14:11:51  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.4  2021/01/07 14:48:33  abbott
// Summary: Corrected copyright
//
// Revision 1.3  2011/03/16 15:33:19  abbott
// Changed arg type "unsigned int" into "unsigned long" for the Intercept* fns.
//
// Revision 1.2  2007/10/30 17:14:11  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.2  2006/10/06 14:04:15  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.1.1.1  2006/05/30 11:39:36  cocoa
// Imported files
//
// Revision 1.2  2006/03/27 12:21:26  cocoa
// Minor silly changes to reduce number of complaints from some compiler or other.
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.1.1.1  2005/05/03 15:47:30  cocoa
// Imported files
//
// Revision 1.3  2005/04/19 14:06:04  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.2  2005/04/01 16:18:20  cocoa
// Friday check-in.  Fixed a bug in the ctor for GeneralQuotientRingImpl.
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.3  2004/11/11 13:02:16  cocoa
// -- added \include *.txt for doxygen
//
// Revision 1.2  2004/11/03 17:58:56  cocoa
// -- improved log messages
// -- minor tidying
//
// Revision 1.1.1.1  2003/09/24 12:55:43  cocoa
// Imported files
//
// Revision 1.2  2002/03/08 18:29:51  abbott
// Added new functions InterceptNew and InterceptDelete.
//
// Revision 1.1  2002/01/31 11:32:18  abbott
// Initial revision
//

#endif
