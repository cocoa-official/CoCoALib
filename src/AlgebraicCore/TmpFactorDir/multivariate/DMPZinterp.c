//   Copyright (c)  1997-2006  John Abbott

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

#include <stdio.h>

#include <stdlib.h>
#include "jalloc.h"
#include "DMPZinterp.h"

#include "DMPZeval.h"
#include "DUPZinterp.h"
#include "Qsolve.h"

DMPZ DMPZinterp(void (*evaluate)(mpz_t, const int*), int *degs)
{
  DMPZ ans, newans;
  int range;
  int i, j, k, d, var, *p, nterms, *exps, rank, *X, row, nrows;
  DMPZ iter;
  DUPZ coeff;
  mpz_t **Y, one;
  Qmat M, rhs, soln;

  exps = (int*)MALLOC(NVARS*sizeof(int));
  for (i=0; i < NVARS; i++) exps[i] = 0;
  mpz_init_set_ui(one, 1);
  ans = DMPZctor(one, exps);
  mpz_clear(one);
//  mpz_init(denom);
  
  for (var=0; var < NVARS; var++)
  {
    nterms = DMPZ_nterms(ans);
#ifdef FACTOR_DEBUG
    printf("Number of terms so far: %d\n", nterms);
    printf("Interp variable x(%d) with degree %d\n", var, degs[var]);
#endif
    M=Qmat_ctor(nterms,nterms);
    rhs=Qmat_ctor(nterms,1);
    soln=Qmat_ctor(nterms,1);
    Y = (mpz_t **)MALLOC(nterms*sizeof(mpz_t*));
    for (i=0; i < nterms; i++)
    {
      Y[i] = (mpz_t*)MALLOC((1+degs[var])*sizeof(mpz_t));
      for (j=0; j <= degs[var]; j++) mpz_init(Y[i][j]);
    }
    range = 2*nterms*nterms;    /* gives a good prob of non-singularity */
    if (range < 10) range = 10; /* avoid very small ranges */
    p = (int*)MALLOC(NVARS*sizeof(int));
    for (i=0; i < NVARS; i++) p[i] = 1+(rand()%range);
#ifdef FACTOR_DEBUG
    for (d=0; d <= degs[var]; d++) printf("."); printf("!\n");
#endif
    for (d=0; d <= degs[var]; d++)
    {
#ifdef FACTOR_DEBUG
      printf("%d",d%10);fflush(stdout);
#endif
      nrows = 0;
      row = 0;
      while (1)
      {
        for (i=0; i < var; i++) p[i] = 1+(rand()%range);
        p[var] = d;
#if 0
printf("eval point: ");for(i=0; i < NVARS; i++)printf("%d ",p[i]);printf("\n");
#endif
        evaluate(mpq_numref(rhs->entry[row][0]), p);
#if 0
printf("evaluate gave ");mpz_out_str(stdout,10,rhs[row][0]);printf("\n");
#endif
        for (j=0,iter=ans; iter; j++,iter=iter->next)
          DMPZeval_power_product(mpq_numref(M->entry[row][j]), iter->exps, p);
        row++;
        nrows++;
        if (nrows < nterms) continue;
//#ifdef FACTOR_DEBUG
//        printf("Matrix is (%dx%d):\n",nterms,nterms);
//        for (i=0;i<nterms;i++){for(j=0;j<nterms;j++){mpz_out_str(stdout,10,mpq_numref(M[i][j]));printf("  ");}printf("\n");}
//        printf("rhs is"); for (i=0;i<nterms;i++){mpz_out_str(stdout,10,mpq_numref(rhs[i][0]));printf("  ");}printf("\n");
//#endif

rank = Qsolve(soln, M, rhs);
/*        OK = Zsolve(soln, &denom, nterms, M, 1, rhs);*/
        if (rank == nterms) break;
#ifdef FACTOR_DEBUG
printf("singular\n");
#endif
        row = rand()%nrows;
        nrows--;
      }
#ifdef FACTOR_DEBUG
      for (i=0; i < nterms; i++)
        if (mpz_cmp_ui(mpq_denref(soln->entry[i][0]), 1) != 0) printf("Ooops, not integer\n");
#endif
      for (i=0; i < nterms; i++)
        mpz_set(Y[i][d], mpq_numref(soln->entry[i][0]));
    }
    FREE(p);
Qmat_dtor(M);
Qmat_dtor(rhs);
Qmat_dtor(soln);
//#ifdef FACTOR_DEBUG
//    printf("\n");
//#endif
    newans = NULL;
    X = (int*)MALLOC((1+degs[var])*sizeof(int));
    for (i=0; i <= degs[var]; i++) X[i] = i;
    coeff = DUPZnew(degs[var]);
    for (i=0, iter=ans; iter; i++, iter=iter->next)
    {
      DUPZinterp(coeff, degs[var], X, Y[i]);
#ifdef FACTOR_DEBUG
printf("DUPZinterp gave: ");DUPZprint(coeff);
#endif
      for (j=DUPZdeg(coeff); j >= 0; j--)
      {
        if (mpz_sgn(coeff->coeffs[j]) == 0) continue;
        exps = (int*)MALLOC(NVARS*sizeof(int));
        for(k=0; k < NVARS; k++) exps[k] = iter->exps[k];
        exps[var] = j;
        newans = DMPZprepend(coeff->coeffs[j], exps, newans);
      }
    }
    DUPZfree(coeff);
    for (i=0; i < nterms; i++)
    {
      for (j=0; j <= degs[var]; j++)
        mpz_clear(Y[i][j]);
      FREE(Y[i]);
    }
    FREE(Y);
    FREE(X);
    DMPZdtor(ans);
    ans = DMPZreverse(newans);
  }
  return ans;
}
