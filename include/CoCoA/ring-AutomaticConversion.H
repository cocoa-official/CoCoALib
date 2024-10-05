#ifndef CoCoA_ring_AutomaticConversion_H
#define CoCoA_ring_AutomaticConversion_H

//   Copyright (c)  2020  John Abbott,  Anna M. Bigatti

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

  class ring;
  class RingHom;
  class ErrorInfo;  // fwd decl

  RingHom AutomaticConversionHom(const ring& R1, const ring& R2, const ErrorInfo& err);

} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/ring-AutomaticConversion.H,v 1.4 2022/02/18 14:11:52 abbott Exp $
// $Log: ring-AutomaticConversion.H,v $
// Revision 1.4  2022/02/18 14:11:52  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.3  2020/06/22 15:42:05  abbott
// Summary: AutomaticConversionHom now expects a ErrorInfo
//
// Revision 1.2  2020/06/20 19:12:16  abbott
// Summary: AutomaticConversionHom now requires 3rd arg (FnName of caller)
//
// Revision 1.1  2020/06/19 14:59:55  abbott
// Summary: Added new fn AutomaticConversionHom
//
//

#endif
