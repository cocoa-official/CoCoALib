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



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/utils-gmp.H,v 1.4 2022/02/18 14:11:52 abbott Exp $
// $Log: utils-gmp.H,v $
// Revision 1.4  2022/02/18 14:11:52  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.3  2021/01/31 09:59:37  abbott
// Summary: Added mpq_get_d_2exp
//
// Revision 1.2  2021/01/07 14:48:34  abbott
// Summary: Corrected copyright
//
// Revision 1.1  2016/03/25 20:39:44  abbott
// Summary: Renamed from utils_gmp to utils-gmp (as otherwise LaTeX crashed when generating doc)
//
// Revision 1.3  2016/03/25 20:01:07  abbott
// Summary: Added new fns mpz_rounddiv & mpq_round
//
// Revision 1.2  2014/06/17 10:04:21  abbott
// Summary: Corrected include directive (gmp.H -> gmp.h)
// Author: JAA
//
// Revision 1.1  2014/06/13 12:08:39  abbott
// Summary: Header for new GMP fn (CmpAbs for rationals)
// Author: JAA
//
//

#endif
