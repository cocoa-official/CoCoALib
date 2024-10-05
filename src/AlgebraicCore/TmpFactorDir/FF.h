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

#ifndef FF_h
#define FF_h

#include "jaa.h"

struct FFstruct
{
  FFelem prime;
  unsigned int *LogTable; /* this may be NULL */
  FFelem       *ExpTable; /* this may be NULL */
  FFelem k, shift;
};

typedef struct FFstruct *FF;

extern struct FFstruct CurrentFF;

FF FFctor(FFelem p);
void FFdtor(FF Fp);
void FFselect(const FF Fp);

FFelem FFchar(const FF Fp);

/*inline*/ FFelem FFadd(const FFelem x, const FFelem y);
/*inline*/ FFelem FFsub(const FFelem x, const FFelem y);
/*inline*/ FFelem FFmul(const FFelem x, const FFelem y);
/*inline*/ FFelem FFdiv(const FFelem x, const FFelem y);

/*inline*/ unsigned int FFlog(FFelem x);
/*inline*/ FFelem FFexp(unsigned int x);

FFelem FFpower(const FFelem x, int e); /* computes x to the power e */

#endif
