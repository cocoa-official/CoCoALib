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

#ifndef DUPZFACTOR_COMBINE_H
#define DUPZFACTOR_COMBINE_H

#include "DUPZ.h"
#include "DUPFF.h"
#include "DUPZfactors.h"
#include "DUPZfactor_bound.h"
#include "DUPZfactor_info.h"

struct DUPZfactor_combine_struct
{
  DUPZfactor_info info; /* information about poly being factored       */
  int nfactors_orig;    /* original number of factors modulo p^k.      */
  int nfactors;         /* number of unused factors modulo p^k.        */
  DUPZ **factors;       /* array of pointers to the factors modulo p^k.*/
  int *used;            /* flags saying which factors are "used"       */
  mpz_t modulus;        /* the value of p^k                            */
  int early;            /* boolean, true iff doing early searching     */
  int single_factor;    /* boolean, true iff single factor bound       */
                        /* < modulus, & all factor bound > modulus     */

  int n1coeffd_flag;
  double n1coeffd_bound;
  double *n1coeffd;
  double *n1coeffd_sum;
  
  int evaluation_point; /* randomly chosen evaluation point            */
  double log_ev_pt;     /* log of eval point plus root bound           */
  mpz_t f_value;        /* value of f at evaluation point              */
  mpz_t *factor_value;  /* value of the factors at evaluation point    */

  mpz_t P, Q;           /* num/den bounds for rational reconstruction  */
  DUPZ factor;          /* a putative true factor                      */
  int tuple_deg;        /* degree of the putative factor               */
  double log_tc;        /* log of trailing coeff of factor             */
  DUPZ quot, rem;       /* used if we actually divide factor in f      */
  mpz_t value_at_1;     /* used in a quick test before dividing        */
  mpz_t value_at_m1;    /*                                             */

  int tuple_size;
  int *combination;
  int *tc_combination;  /* to avoid repeated multiplications in tctest */
  mpz_t *tc_product;    /*                                             */
};

typedef struct DUPZfactor_combine_struct *DUPZfactor_combiner;


DUPZfactor_combiner DUPZfactor_combiner_ctor(DUPZfactor_info info);
void DUPZfactor_combiner_dtor(DUPZfactor_combiner THIS);

void DUPZfactor_combine_init();
void DUPZfactor_combine_done();
void DUPZfactor_combine_early(DUPZfactor_combiner THIS);
void DUPZfactor_combine_final(DUPZfactor_combiner THIS);



#endif
