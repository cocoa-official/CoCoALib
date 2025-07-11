#ifndef CoCoA_ring_AutomaticConversion_H
#define CoCoA_ring_AutomaticConversion_H

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


#include "CoCoA/error.H"

namespace CoCoA
{

  class ring;
  class RingHom;
  class ErrorInfo;  // fwd decl

  RingHom AutomaticConversionHom(const ring& R1, const ring& R2, const ErrorContext& ErrCtx);

} // end of namespace CoCoA

#endif
