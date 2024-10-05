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

#include "DUPIsubtract_product.h"
#include "jaaerror.h"

/***************************************************************************/
/* this computes f -= g*h for DUPI f and DUPFFs g,h; for Hensel lifting    */

void DUPIsubtract_product(DUPI f, const DUPFF g, const DUPFF h)
{
  int i, dg, dh, dgh;
  int32 *fc, *fci;
  uint32 *gc, *hc, *gci, *hci;

  dg = DUPFFdeg(g);
  dh = DUPFFdeg(h);
  if (dg < 0 || dh < 0) return;

  dgh = dg+dh;
  if (f->maxdeg < dgh)  /* error case */
  {
    JERROR(JERROR_DEG_TOO_LOW);
    return;
  }

  for (i=DUPIdeg(f) + 1; i <= dgh; i++) f->coeffs[i] = 0;

  fc = f->coeffs + dg;
  gc = g->coeffs;
  hc = h->coeffs;

  for (gci=gc+dg; gci>=gc; gci--, fc--)
  {
    if (!*gci) continue;
    for (hci=hc+dh, fci=fc+dh; hci>=hc; hci--, fci--)
    {
      *fci -= *gci * *hci;
    }
  }

  /* fix degree of f */
  if (DUPIdeg(f) > dgh) dgh = DUPIdeg(f);
  while (dgh >= 0 && f->coeffs[dgh] == 0) dgh--;
  f->deg = dgh;
  return;
}
