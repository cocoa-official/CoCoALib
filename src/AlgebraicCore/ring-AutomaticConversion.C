//   Copyright (c)  2020,2025  John Abbott,  Anna M. Bigatti

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

  // Throws ERR::MixedRings with given context if there is no CanonicalHom
  RingHom AutomaticConversionHom(const ring& R1, const ring& R2, const ErrorContext& ErrCtx)
  {
    if (R1 == R2)
      CoCoA_THROW_ERROR_WITH_CONTEXT3(ERR::IncompatArgs, "rings must be different", ErrCtx); // ??? get FILE/LINE from err???
    try
    {
      if (RingID(R1) < RingID(R2))
        return CanonicalHom(R1,R2);
      return CanonicalHom(R2,R1);
    }
    catch (const ErrorInfo&) // e.g. we do not catch Timeouts & interrupts
    {
      CoCoA_THROW_ERROR_WITH_CONTEXT2(ERR::MixedRings, ErrCtx);
    }
    return IdentityHom(R1); // NEVER EXECUTED, just to keep compiler quiet
  }


} // end of namespace CoCoA
