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
