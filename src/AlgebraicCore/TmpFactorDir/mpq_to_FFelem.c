//   Copyright (c)  1997-2006  John Abbott,  Anna M Bigatti

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

#include "mpq_to_FFelem.h"

int mpq_to_FFelem(FFelem *r, mpq_t q)
{
  FFelem p = CurrentFF.prime;
  FFelem num, den;

  den = mpz_fdiv_ui(mpq_denref(q), p);
  if (den == 0) return 0;
  num = mpz_fdiv_ui(mpq_numref(q), p);
  *r = FFdiv(num, den);
  return 1;
}
