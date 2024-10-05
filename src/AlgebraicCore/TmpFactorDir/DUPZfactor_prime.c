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

#include "jalloc.h"
#include "DUPZfactor_prime.h"
#include "DUPZ_DUPFF.h"
#include "DUPFFsqfrd.h"
#include "DUPFFfactor.h"
#include "primes.h"
#include "DUPZfactor_refine_fds.h"


/***************************************************************************/
/* This function is called to create an initial collection of modular      */
/* factorizations, and with it an initial factor degree set.               */
/* Although long, this code is quite simple.                               */

int DUPZfactor_add_prime(DUPZfactor_info THIS)
{
  DUPZ f = THIS->f; /* just an alias */
  FF FFq;
  DUPFFlist qfactors, iter;
  int q, nprimes;
  DUPFF fq;

  nprimes = THIS->nprimes;
  if (nprimes == THIS->max_nprimes) return 0;
  if (nprimes > 0) q = FFchar(THIS->FFq[nprimes-1]);
  else q = 7;

  /* skip primes dividing lc(f), tc(f), or discr(f) */
  while (1)
  {
    q = nextprime(q);
    if (mpz_fdiv_ui(DUPZlc(f), q) == 0) continue;
    if (mpz_fdiv_ui(f->coeffs[0], q) == 0) continue;

    FFq = FFctor(q);
    FFselect(FFq);
    fq = DUPZ_to_DUPFF(f);
    if (DUPFFsqfr(fq)) break;  /* OK, this prime is good */
    DUPFFfree(fq);
    FFdtor(FFq);
  }

  /* make f monic and factorize mod p */
  DUPFFmake_monic(fq);
  qfactors = DUPFFlist_sort(DUPFFfactor(fq));
  DUPFFfree(fq);
  
  /* make all the factors monic */
  for (iter = qfactors; iter; iter = iter->next)
    DUPFFmake_monic(iter->poly);

#ifdef FACTOR_DEBUG
  printf("Degree pattern for prime %d: ", q);
  for (iter=qfactors; iter; iter = iter->next)
    printf("%d ",DUPFFdeg(iter->poly));
  printf("\n");
#endif

  /* store factorization mod p for possible later use */
  THIS->FFq[nprimes]      = FFq;
  THIS->qfactors[nprimes] = qfactors;
  THIS->nprimes++;

  /* Update factor degree set using factorization for this q */
  return DUPZfactor_refine_fds(THIS->fds, qfactors);
}


/***************************************************************************/

int DUPZfactor_pick_prime(DUPZfactor_info THIS)
{
  DUPFF fi; /* only used as an alias */
  DUPFFlist qfactors, iter;
  int i, degfi, nprimes, chosen;
  int *n1count, *delta, *nfactors;
  int df;
  int init_nprimes;

  df = DUPZdeg(THIS->f);
  /* Decide on initial number of primes to try, roughly log2(df)-1 */
  /* This is just a heuristic value with some empirical support.   */
  for (init_nprimes=3; df >> (init_nprimes+1); init_nprimes++) {}
  --init_nprimes;

  for (i=0; i < init_nprimes; i++)
    if (DUPZfactor_add_prime(THIS)) break; /* break if proved irreducible */
  if (i < init_nprimes) return 1;


  for (THIS->dmax = df/2; !THIS->fds[THIS->dmax]; THIS->dmax--) {}

  nprimes = THIS->nprimes;
  n1count = (int*)MALLOC(nprimes*sizeof(int));
  delta = (int*)MALLOC(nprimes*sizeof(int));
  nfactors = (int*)MALLOC(nprimes*sizeof(int));
  for (i=0; i < nprimes; i++)
  {
    qfactors = THIS->qfactors[i];
    nfactors[i] = DUPFFlist_length(qfactors);
    n1count[i] = 0;
    for (delta[i] = 1;; delta[i]++)
    {
      for (iter = qfactors; iter; iter = iter->next)
      {
        fi = iter->poly;
        degfi = DUPFFdeg(fi);
        n1count[i] += (fi->coeffs[degfi - delta[i]] != 0);
      }
      if (n1count[i] > 0) break;
    }
  }

  
  chosen = 0;
  for (i=1; i < nprimes; i++)
  {
    if (nfactors[i] > nfactors[chosen]) continue;
    if (nfactors[i] == nfactors[chosen] && n1count[i] < n1count[chosen]) continue;
    if (nfactors[i] == nfactors[chosen] &&
        n1count[i] == n1count[chosen] &&
        delta[i] < delta[chosen]) continue;
    chosen = i;
  }
  
  FFselect(THIS->FFq[chosen]);
  THIS->p = FFchar(THIS->FFq[chosen]);
  THIS->p_index = chosen;
  THIS->pfactors = THIS->qfactors[chosen];
//  THIS->delta = delta[chosen];
#ifdef FACTOR_DEBUG
  printf("The chosen prime is %d\n", THIS->p);
#endif
  FREE(n1count); FREE(delta); FREE(nfactors);

  /* These three lines should be part of Hensel initialization? */
  THIS->current_height = 1;
  mpz_set_ui(THIS->Q, THIS->p);
  mpz_set_ui(THIS->recip_lcf, FFdiv(1, mpz_fdiv_ui(DUPZlc(THIS->f), THIS->p)));
  return 0;
}

