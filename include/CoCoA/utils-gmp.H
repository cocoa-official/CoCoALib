#ifndef CoCoA_utils_gmp_H
#define CoCoA_utils_gmp_H

//   Copyright (c)  2014,2016,2021  John Abbott, and Anna M. Bigatti

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

#include "gmp.h"

namespace CoCoA
{

  int mpq_cmpabs(const mpq_t q1, const mpq_t q2);
  bool mpz_rounddiv(mpz_t Q, const mpz_t N, const mpz_t D); // result == true iff division was exact
  bool mpq_round(mpz_t N, const mpq_t Q); // result == true iff Q is an integer
  double mpq_get_d_2exp(signed long int* exp, const mpq_t Q);

} // end of namespace CoCoA

#endif
