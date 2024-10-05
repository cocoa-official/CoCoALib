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

#include "DMPZ_to_DUPFF.h"
#include "jaaerror.h"

static FFelem FFpow(FFelem b, int e)
{
  if (e == 0) return 1;
  if (e < 0) return FFpow(FFdiv(1, b), -e);
  if (e == 1) return b;
  if (e&1) return FFmul(b, FFpow(FFmul(b, b), (e-1)/2));
  return FFpow(FFmul(b, b), e/2);
}


DUPFF DMPZ_to_DUPFF(DMPZ f, int var, int *substitution)
{
  FFelem p = CurrentFF.prime;
  int df, d, c, i;
  DUPFF ans;
  DMPZ iter;
  
  df = DMPZdeg(f, var);
  ans = DUPFFnew(df);
  for(i=0; i<=df;i++) ans->coeffs[i] = 0;
  for(iter=f;iter;iter = iter->next)
  {
    d = iter->exps[var];
    c = mpz_fdiv_ui(iter->coeff, p);
    for(i=0; c!=0 && i<NVARS;i++)
      if (i != var)
        c = FFmul(c, FFpow(substitution[i], iter->exps[i]));
    ans->coeffs[d] += c;
    if (ans->coeffs[d] >= p) ans->coeffs[d] -= p;
  }
  for(i=df;i>=0;i--) if (ans->coeffs[i] != 0) break;
  ans->deg = i;
  return ans;
}
