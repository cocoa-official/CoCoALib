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

#include "jaaerror.h"
#include "DUPFF.h"
#include "DUPFFmod.h"

void DUPFFmul3mod(DUPFF ans, const DUPFF x, const DUPFF y, const DUPFF m)
{
  if (ans == m) { JERROR(JERROR_ALIASING); return; }
  if (ans->maxdeg < DUPFFdeg(x) + DUPFFdeg(y))
  {
    JERROR(JERROR_DEG_TOO_LOW);
    ans->deg = -1;
    return;
  }
  DUPFFmul3(ans, x, y);
  DUPFFrem2(ans, m);
}


void DUPFFexpt3mod(DUPFF ans, const DUPFF base, const long power, const DUPFF m)
{
  if (power < 1) { JERROR(JERROR_EXPT_MOD); ans->deg = -1; return; }
  if (ans->maxdeg < 2*DUPFFdeg(base)-2)
  {
    JERROR(JERROR_DEG_TOO_LOW);
    ans->deg = -1;
    return;
  }
  if (power == 1)
  {
    DUPFFcopy2(ans, base);
    return;
  }
  DUPFFexpt3mod(ans, base, power/2, m);
  DUPFFmul3mod(ans, ans, ans, m);
  if ((power&1) == 0) return;
  DUPFFmul3mod(ans, ans, base, m);
}
