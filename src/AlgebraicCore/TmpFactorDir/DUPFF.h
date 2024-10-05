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

#ifndef DUPFF_h
#define DUPFF_h

#include "FF.h"

struct DUPFFstruct
{
  int maxdeg;
  int deg;
  FFelem *coeffs;
};

typedef struct DUPFFstruct *DUPFF;

DUPFF DUPFFnew(const int maxdeg);
void DUPFFfree(DUPFF x);
void DUPFFswap(DUPFF x, DUPFF y);
DUPFF DUPFFcopy(const DUPFF f);
void DUPFFcopy2(DUPFF dest, const DUPFF src);
int DUPFFequal(const DUPFF x, const DUPFF y);

int DUPFFdeg(const DUPFF f);
FFelem DUPFFlc(const DUPFF f);
void DUPFFmake_monic(DUPFF f);
void DUPFFdiv2ff(DUPFF f, const FFelem c);
void DUPFFmul2ff(DUPFF f, const FFelem c);


DUPFF DUPFFadd(const DUPFF x, const DUPFF y);
void DUPFFadd3(DUPFF sum, const DUPFF x, const DUPFF y);
DUPFF DUPFFsub(const DUPFF x, const DUPFF y);
void DUPFFsub3(DUPFF sum, const DUPFF x, const DUPFF y);
DUPFF DUPFFmul(const DUPFF x, const DUPFF y);
void DUPFFmul3(DUPFF ans, const DUPFF x, const DUPFF y);
void DUPFFsquare(DUPFF f);
DUPFF DUPFFexpt(const DUPFF base, const int power);
void DUPFFexpt3(DUPFF ans, const DUPFF base, const int power);


void DUPFFshift_add(DUPFF f, const DUPFF g, int deg, const FFelem coeff);
void DUPFFshift_add_raw(FFelem *dest, const FFelem *src, const FFelem *srclast, const FFelem scalar);

DUPFF DUPFFdiv(const DUPFF num, const DUPFF den);
DUPFF DUPFFrem(const DUPFF num, const DUPFF den);
void DUPFFrem2(DUPFF x, const DUPFF m);
void DUPFFdiv4(DUPFF quot, DUPFF rem, const DUPFF num, const DUPFF den);

DUPFF DUPFFgcd(const DUPFF fin, const DUPFF gin);
void DUPFFgcd2(DUPFF f, DUPFF g); /* DESTRUCTIVE */
DUPFF DUPFFexgcd(DUPFF *fcofac, DUPFF *gcofac, const DUPFF f, const DUPFF g);

FFelem DUPFFeval(const DUPFF f, FFelem x);

/* To avoid problems on Macintoshes etc which don't have printf... */
#ifdef FACTOR_DEBUG
#include "DUPFFprint.h"
#endif

#endif

