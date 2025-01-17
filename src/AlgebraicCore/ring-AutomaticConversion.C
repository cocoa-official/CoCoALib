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
