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

#ifndef jaa_h
#define jaa_h

#include <limits.h>

/* These definitions are PLATFORM DEPENDENT.                        */
/* If your machine does not have 8-bit bytes then you must alter    */
/* the 8 in WGD.c (there does not seem to be a portable way of      */
/* knowing how many bits are in a long, for instance).              */

/* Try to decide automatically whether we are on 32-bit or 64-bit platform. */
/* See also further explanations below.                                     */


/* The obvious test is the line below, but constants bigger than 32767 are   */
/* not portable, so instead we divide by 1000 twice then check the quotient. */
/*#if ULONG_MAX == 4294967295UL*/
#if (ULONG_MAX/1000)/1000 == 4294

/* This is for a 32-bit machine */
typedef unsigned long FFelem;
#define MAX_PRIME 46341UL
#define MAX_FFelem 4294967295UL

#else

/* This is for a 64-bit machine */

typedef unsigned long FFelem;
#define MAX_PRIME 3037000499UL
#define MAX_FFelem 18446744073709551615UL
#endif

#if 0
/* This is for a 32-bit machine with fast double-precision floating point */
THIS DOES NOT CURRENTLY WORK!!  PLENTY OF PROBLEMS LATER ON!!
typedef double FFelem;
#define MAX_PRIME 67108863.0
#define MAX_FFelem 4503599627370495.0
#endif

/* COMMENTS                                                             */
/* MAX_PRIME should always satisfy 2*(MAX_PRIME-1)^2 < MAX_FFelem       */
/* otherwise changes will be needed to several routines (e.g FFdet).    */
/* So that two FFelems can be summed easily we impose the condition that*/
/* twice the modulus can be represented as an FFelem.  So that two such */
/* elements can be multiplied easily we impose the condition that the   */
/* square of the modulus fits inside an FFelem.                         */
/* The GMP package imposes some restrictions on what we can choose as   */
/* the type of FFelem because we want to use functions like mpz_fdiv_ui.*/
/* So the type chosen for an FFelem must be interconvertible with an    */
/* "unsigned long" without any loss of precision for integer values in  */
/* the range zero up to MAX_PRIME (search for functions with names like */
/* mpz_fdiv_ui, mpz_add_ui, etc.).                                      */
/* If MAX_PRIME is set to a value larger than 2^32 then some changes    */
/* will be needed in the file FindPrimRoot.c (not easy if efficiency is */
/* important).                                                          */


#endif
