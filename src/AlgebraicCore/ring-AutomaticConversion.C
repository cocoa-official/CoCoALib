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


#include "CoCoA/ring-AutomaticConversion.H"
#include "CoCoA/error.H"
#include "CoCoA/CanonicalHom.H"
#include "CoCoA/RingHom.H"


namespace CoCoA
{

  // /// WHAT TO DO HERE???  How much of errto preserve, how much to change???
  RingHom AutomaticConversionHom(const ring& R1, const ring& R2, const ErrorInfo& err)
  {
    if (R1 == R2) CoCoA_THROW_ERROR(ERR::IncompatArgs, "AutomaticConversionHom: rings must be different"); // ??? get FILE/LINE from err???
    try
    {
      if (RingID(R1) < RingID(R2))
        return CanonicalHom(R1,R2);
      return CanonicalHom(R2,R1);
    }
    catch (const ErrorInfo&) // e.g. we do not catch Timeouts & interrupts
    {
      ThrowException(err);
    }
    return IdentityHom(R1); // NEVER EXECUTED, just to keep compiler quiet
  }


  // RingHom AutomaticConversionHom(const ring& R1, const ring& R2, const char* const FnName)
  // {
  //   if (R1 == R2) CoCoA_THROW_ERROR(ERR::IncompatArgs, "AutomaticConversionHom: rings must be different");
  //   try
  //   {
  //     if (RingID(R1) < RingID(R2))
  //       return CanonicalHom(R1,R2);
  //     return CanonicalHom(R2,R1);
  //   }
  //   catch (const ErrorInfo&)
  //   {
  //     CoCoA_THROW_ERROR("Automatic ring conversion not possible", FnName);
  //   }
  //   return IdentityHom(R1); // NEVER EXECUTED, just to keep compiler quiet
  // }


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/ring-AutomaticConversion.C,v 1.6 2023/07/04 09:20:13 abbott Exp $
// $Log: ring-AutomaticConversion.C,v $
// Revision 1.6  2023/07/04 09:20:13  abbott
// Summary: Changed TimeOut to Timeout
//
// Revision 1.5  2022/02/18 14:12:02  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.4  2020/06/22 15:43:11  abbott
// Summary: AutomaticConversionHom now expects a ErrorInfo obj
//
// Revision 1.3  2020/06/20 19:12:17  abbott
// Summary: AutomaticConversionHom now requires 3rd arg (FnName of caller)
//
// Revision 1.2  2020/06/19 19:39:42  abbott
// Summary: Added useles line to keep compiler quiet
//
// Revision 1.1  2020/06/19 14:59:55  abbott
// Summary: Added new fn AutomaticConversionHom
//
//
