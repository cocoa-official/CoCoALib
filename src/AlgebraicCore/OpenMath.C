//   Copyright (c)  2005-2008  John Abbott and Anna M. Bigatti

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


#include "CoCoA/OpenMath.H"
#include "CoCoA/MachineInt.H"
#include "CoCoA/error.H"

#include <iostream>
using std::ostream;

namespace CoCoA
{

  OpenMathSymbol::OpenMathSymbol():
      myCD("--unset--"),
      myName("--unset--")
  {}


  OpenMathSymbol::OpenMathSymbol(const char* const cd, const char* const name):
      myCD(cd),
      myName(name)
  {}

  OpenMathSymbol::OpenMathSymbol(const std::string& cd, const std::string& name):
      myCD(cd),
      myName(name)
  {}


  void OpenMathSymbol::myOutputSelf(std::ostream& out) const
  {
    out << "OpenMathSymbol(" << myCD << ", " << myName << ")";
  }


  const std::string& CD(const OpenMathSymbol& s)
  {
    return s.myCD;
  }

  const std::string& name(const OpenMathSymbol& s)
  {
    return s.myName;
  }


  std::ostream& operator<<(ostream& out, const OpenMathSymbol& oms)
  {
    oms.myOutputSelf(out);
    return out;
  }

  //-------------------------------------------------------

  OpenMathOutputBase::OpenMathOutputBase():
    IntrusiveReferenceCount()
  {}


  OpenMathOutputBase::~OpenMathOutputBase()
  {}


  OpenMathInputBase::OpenMathInputBase():
    IntrusiveReferenceCount()
  {}


  OpenMathInputBase::~OpenMathInputBase()
  {}



// Removed to permit auto conversion from bool to bool3
//   OpenMathOutput& operator<<(OpenMathOutput& OMOut, const MachineInt& n)
//   {
//     OMOut->mySend(n);
//     return OMOut;
//   }
  OpenMathOutput& operator<<(OpenMathOutput& OMOut, int n)
  {
    OMOut->mySend(n);
    return OMOut;
  }
  OpenMathOutput& operator<<(OpenMathOutput& OMOut, unsigned int n)
  {
    OMOut->mySend(n);
    return OMOut;
  }
  OpenMathOutput& operator<<(OpenMathOutput& OMOut, long n)
  {
    OMOut->mySend(n);
    return OMOut;
  }
  OpenMathOutput& operator<<(OpenMathOutput& OMOut, unsigned long n)
  {
    OMOut->mySend(n);
    return OMOut;
  }


  OpenMathOutput& operator<<(OpenMathOutput& OMOut, const OpenMathSymbol& s)
  {
    OMOut->mySend(s);
    return OMOut;
  }

  //---------------------------------------------------------------------------

  OpenMathInput& operator>>(OpenMathInput& OMIn, long n)
  {
    if (!OMIn->myRecv(n))
      CoCoA_THROW_ERROR(ERR::BadOpenMath, "reading machine integer");
    return OMIn;
  }

  OpenMathInput& operator>>(OpenMathInput& OMIn, OpenMathSymbol& s)
  {
    if (!OMIn->myRecv(s))
      CoCoA_THROW_ERROR(ERR::BadOpenMath, "reading OM symbol");
    return OMIn;
  }


} // end of namespace CoCoA

// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/OpenMath.C,v 1.10 2022/02/18 14:11:56 abbott Exp $
// $Log: OpenMath.C,v $
// Revision 1.10  2022/02/18 14:11:56  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.9  2021/01/07 15:16:51  abbott
// Summary: Corrected copyright
//
// Revision 1.8  2020/06/17 15:49:24  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.7  2012/05/04 16:59:43  abbott
// For future compatibility with auto conversion from bool to 3-way bool,
// removed sending of MachineInt; instead added operators for sending
// (unsigned) int and (unsigned) long [but no other integral types for now].
//
// Revision 1.6  2011/11/09 14:09:53  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.5  2011/03/11 14:49:08  abbott
// Changed size_t into long.
//
// Revision 1.4  2010/03/18 13:53:48  abbott
// Minor rationalization to OpenMath implementation: moved op<<(OMOut,ZZ) to ZZ files.
//
// Revision 1.3  2008/12/16 21:10:32  abbott
// Replaced the various output fns for different sort of machine integers by a
// single one for MachineInt.
//
// Revision 1.2  2007/10/30 17:14:08  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.4  2006/11/27 13:04:35  cocoa
// Added explicit OpenMath output operators for all integer types.
//
// Revision 1.3  2006/11/23 17:10:52  cocoa
// -- changed: OpenMathOutput and OpenMathInput are now a class (instead of typedef)
//
// Revision 1.2  2006/11/02 13:25:44  cocoa
// Simplification of header files: the OpenMath classes have been renamed.
// Many minor consequential changes.
//
