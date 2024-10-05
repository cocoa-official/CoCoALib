#ifndef CoCoA_PPMonoidSparse_H
#define CoCoA_PPMonoidSparse_H

//   Copyright (c)  2005,2007,2015  John Abbott,  Anna M. Bigatti

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


#include <vector>
// using vector;

namespace CoCoA
{

  class symbol;         // fwd decl -- defined in symbol.H
  class PPMonoid;       // fwd decl -- defined in PPMonoid.H
  class PPOrdering;     // fwd decl -- defined in PPOrdering.H
  class PPOrderingCtor; // fwd decl -- defined in PPOrdering.H

  // Pseudo-ctor
  PPMonoid NewPPMonoidSparse(const std::vector<symbol>& IndetNames, const PPOrdering& ord);
  PPMonoid NewPPMonoidSparse(const std::vector<symbol>& IndetNames, const PPOrderingCtor& ord);

} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/PPMonoidSparse.H,v 1.6 2022/02/18 14:11:47 abbott Exp $
// $Log: PPMonoidSparse.H,v $
// Revision 1.6  2022/02/18 14:11:47  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.5  2021/01/07 14:21:05  abbott
// Summary: Corrected copyright
//
// Revision 1.4  2015/11/30 21:53:55  abbott
// Summary: Major update to matrices for orderings (not yet complete, some tests fail)
//
// Revision 1.3  2011/06/23 16:07:13  abbott
// Added incomplete but compiling PPMonoidSparse: first prototype,
// simple rather than fast!
//
// Revision 1.2  2007/10/30 17:14:12  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.5  2007/03/08 17:43:11  cocoa
// Swapped order of args to the NewPPMonoid pseudo ctors.
//
// Revision 1.4  2006/11/27 13:41:32  cocoa
// -- added #include <cstddef> //using std::size_t
// -- removed #include "config.H"
//
// Revision 1.3  2006/11/24 17:32:49  cocoa
// -- reorganized includes of header files
//
// Revision 1.2  2006/10/06 14:04:15  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.6  2006/04/28 16:33:51  cocoa
// Used SmartPtrIRC for PPOrderings.
//
// Revision 1.5  2006/04/27 13:45:30  cocoa
// Changed name of NewIdentityRingHom to NewIdentityHom.
// Changed name of member functions which print out their own object
// into myOutputSelf (to distinguish from "transitive" myOutput fns).
//
// Revision 1.4  2006/03/12 21:28:34  cocoa
// Major check in after many changes
//
// Revision 1.3  2006/02/20 22:41:20  cocoa
// All forms of the log function for power products now return SmallExponent_t
// (instead of int).  exponents now resizes the vector rather than requiring
// the user to pass in the correct size.
//
// Revision 1.2  2006/01/17 10:23:08  cocoa
// Updated DivMask; many consequential changes.
// A few other minor fixes.
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.1  2005/09/22 18:04:17  cocoa
// It compiles; the tests run OK.  The examples compile.
// No documentation -- the mindless eurocrats have rendered
// me mindless too.
//

#endif
