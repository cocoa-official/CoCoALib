#ifndef CoCoA_CanonicalHom_H
#define CoCoA_CanonicalHom_H

//   Copyright (c)  2007,2009  John Abbott,  Anna M. Bigatti

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


namespace CoCoA
{
  class ring;    // fwd decl -- defined in ring.H
  class RingHom; // fwd decl -- defined in RingHom.H

  RingHom CanonicalHom(const ring& domain, const ring& codomain);
  RingHom ChainCanonicalHom(const ring& domain, const ring& codomain);

} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/CanonicalHom.H,v 1.7 2022/02/18 14:11:45 abbott Exp $
// $Log: CanonicalHom.H,v $
// Revision 1.7  2022/02/18 14:11:45  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.6  2021/01/07 14:10:46  abbott
// Summary: Corrected copyright
//
// Revision 1.5  2018/05/17 16:05:43  bigatti
// -- renamed TmpChainCanonicalHom --> ChainCanonicalHom
//
// Revision 1.4  2011/02/18 12:56:08  bigatti
// -- added TmpChainCanonicalHom
//
// Revision 1.3  2009/07/24 14:21:04  abbott
// Cleaned up include directives, and added some fwd decls.
//
// Revision 1.2  2007/10/30 17:14:12  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.1  2007/03/05 21:25:02  cocoa
// New CanonicalHom pseudo-ctor.
//
//

#endif
