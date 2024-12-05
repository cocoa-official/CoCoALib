//   Copyright (c)  2012-2013-2019  John Abbott, Anna M. Bigatti

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


#include "CoCoA/SparsePolyOps-ideal-points.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/BigRatOps.H"
#include "CoCoA/CanonicalHom.H"
#include "CoCoA/DenseMatrix.H"
#include "CoCoA/interrupt.H"
#include "CoCoA/MatrixOps.H"
#include "CoCoA/MatrixView.H"
#include "CoCoA/NumTheory-gcd.H"
#include "CoCoA/PPMonoid.H"
#include "CoCoA/QBGenerator.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/SparsePolyOps-ideal.H"
#include "CoCoA/SparsePolyOps-RingElem.H"
//#include "CoCoA/SparsePolyRing.H" // from SparsePolyOps-RingElem.H
#include "CoCoA/ideal.H"
#include "CoCoA/matrix.H"
#include "CoCoA/ring.H"
#include "CoCoA/utils.H"
#include "TmpFactorDir/BM.h"
#include "TmpFactorDir/FF.h"

// The old C4 code still uses malloc/free...
#include <stdlib.h>

//#include <vector>
using std::vector;

#include <iostream>


namespace CoCoA
{

  namespace // anonymous namespace for file local fns
  {
    bool IsValidSolution(const matrix& M)
    {
      return !(NumRows(M) == 0 && NumCols(M) == 0);
    }


  vector<RingElem> eval(const PPMonoidElem& t, const ConstMatrixView& pts)
  {
    const ring k = RingOf(pts);
    const long NumPts = NumRows(pts);
    const long NumVars = NumIndets(owner(t));
    const vector<long> expv = exponents(t);
    vector<RingElem> ans; ans.reserve(NumPts);
    for (long i=0; i < NumPts; ++i)
    {
      RingElem val = one(k);
      for (long j=0; j < NumVars; ++j)
        if (expv[j] > 0) // to avoid accessing non-existent columns in matrix pts!!!
          val *= power(pts(i,j), expv[j]);
      ans.push_back(val);
    }
    return ans;
  }

  }


  // Simple rather than efficient (esp. the call to eval)
  std::vector<RingElem> BM_generic(const SparsePolyRing& P, const ConstMatrixView& pts)
  {
    if (CoeffRing(P) != RingOf(pts))  CoCoA_THROW_ERROR1(ERR::MixedRings);
    if (NumIndets(P) < NumCols(pts))  CoCoA_THROW_ERROR1(ERR::IncompatDims);

    const long NumPts = NumRows(pts);
    const long dim = NumCols(pts);
    const ring k = CoeffRing(P);

    vector<RingElem> GB;
    const PPMonoid TT = PPM(P);
    QBGenerator QBG(TT);
    QBG.myCornerPPIntoQB(one(TT));
    matrix M = NewDenseMat(k, 1, NumPts);
    // Fill first row with 1:
    for (int i=0; i<NumPts; ++i)  SetEntry(M,0,i, 1);

    // The next loop removes the last indets from consideration.
    for (int i=dim; i < NumIndets(TT); ++i)
      QBG.myCornerPPIntoAvoidSet(indet(TT,i));

    while (!QBG.myCorners().empty())
    {
      CheckForInterrupt("BM_generic");
      const PPMonoidElem t = QBG.myCorners().front();
      const vector<RingElem> v = eval(t, pts);
      ConstMatrixView NewRow = RowMat(v);
      const matrix a = LinSolve(transpose(M), transpose(NewRow));
      if (IsValidSolution(a))
      {
        QBG.myCornerPPIntoAvoidSet(t);
        RingElem NewGBElem = monomial(P, t);
        const vector<PPMonoidElem>& QB =  QBG.myQB();
        const long nrows = NumRows(M);
        for (int i=0; i < nrows; ++i)
          NewGBElem -= monomial(P, a(i,0), QB[i]);
        GB.push_back(NewGBElem);
      }
      else
      {
        QBG.myCornerPPIntoQB(t);
        M = NewDenseMat(ConcatVer(M, NewRow));
      }
    }
    return GB;
  }


// TRULY HORRIBLE ANCIENT CODE WRITTEN BY YOURS TRULY... SIGH!!!
static const PPMonoid* pp_cmp_PPM;
/* ASSUMES the two args are never equal */
static int pp_cmp(const void *arg1, const void *arg2)
{
//  int var, ans;
  const int **pp1 = (const int**)arg1;
  const int **pp2 = (const int**)arg2;
//  ring R = pp_cmp_ring;
  const int nvars = NumIndets(*pp_cmp_PPM);
//  term PP1, PP2;

//  PP1 = term_init(R->indetsNo);
  vector<long> expv1(nvars);
  vector<long> expv2(nvars);
  for (int var = 0; var < nvars; ++var)
  {
    expv1[var] = pp1[0][var];
    expv2[var] = pp2[0][var];
  }
  return cmp(PPMonoidElem(*pp_cmp_PPM, expv1), PPMonoidElem(*pp_cmp_PPM, expv2));
//     term_put_nth(PP1, (var+1), pp1[0][var] * iv_get_nth(R->weights, var+1));
//   PP2 = term_init(R->indetsNo);
//   for (var = 0; var < nvars; var++)
//     term_put_nth(PP2, (var+1), pp2[0][var] * iv_get_nth(R->weights, var+1));
//   ans = term_is_gt(PP1, PP2, R);
//   term_free(PP1); term_free(PP2);
//   if (ans) return 1;
//   return -1;
}

  std::vector<RingElem> BM_modp(const SparsePolyRing& P, const ConstMatrixView& pts)
  {
    ring Fp = CoeffRing(P);

    const int NumPts = NumRows(pts);
    const int NumVars = NumCols(pts);
    const long p = ConvertTo<long>(characteristic(Fp));
    FF FFp = FFctor(p);
    FFselect(FFp);
    FFelem** points_p = (FFelem**)malloc(NumPts*sizeof(FFelem*));
    for (int i=0; i < NumPts; ++i)
    {
      points_p[i] = (FFelem*)malloc(NumVars*sizeof(FFelem));
      for (int j=0; j < NumVars; ++j)
      {
        points_p[i][j] = LeastNNegRemainder(ConvertTo<BigInt>(pts(i,j)), p);
      }
    }
    pp_cmp_PPM = &PPM(P);
    const BM modp = BM_affine_mod_p(NumVars, NumPts, points_p, pp_cmp);
    for (int i=0; i < NumPts; ++i)  { free(points_p[i]); }    free(points_p);
    if (modp == nullptr)  { FFdtor(FFp); return std::vector<RingElem>(); } // empty list means error

    const int GBsize = modp->GBsize;
    std::vector<RingElem> GB(GBsize);
    vector<long> expv(NumVars);
    for (int i=0; i < GBsize; ++i)
    {
      CheckForInterrupt("BM_modp");
      for (int var = 0; var < NumVars; ++var)
        expv[var] = modp->pp[modp->GB[i]][var];
      RingElem GBelem = monomial(P, expv);
      for (int j=0; j < NumPts; ++j)
      {
        const int c = modp->M[modp->GB[i]][j+NumPts];
        if (c == 0) continue;
        for (int var = 0; var < NumVars; ++var)
          expv[var] = modp->pp[modp->sep[j]][var];
        GBelem += monomial(P, c, expv);
      }
      GB[i] = GBelem;
    }
    BM_dtor(modp);
    FFdtor(FFp);
    return GB;
  }

  namespace
  {
    int deg(int nvars, int *pp)
    {
      int ans=0; 
      for (int i=0; i < nvars; i++) ans += pp[i];
      return ans;
    }
  }
  

  std::vector<RingElem> BM_projective_modp(const SparsePolyRing& P, const ConstMatrixView& pts)
  {
    ring Fp = CoeffRing(P);

    const int NumPts = NumRows(pts);
    const int NumVars = NumCols(pts);
    const long p = ConvertTo<long>(characteristic(Fp));
    FF FFp = FFctor(p);
    FFselect(FFp);
    FFelem** points_p = (FFelem**)malloc(NumPts*sizeof(FFelem*));
    for (long i=0; i < NumPts; ++i)
    {
      points_p[i] = (FFelem*)malloc(NumVars*sizeof(FFelem));
      for (long j=0; j < NumVars; ++j)
      {
        points_p[i][j] = LeastNNegRemainder(ConvertTo<BigInt>(pts(i,j)), p);
      }
    }
    pp_cmp_PPM = &PPM(P);
    const BM modp = BM_projective_mod_p(NumVars, NumPts, points_p, pp_cmp);
    for (int i=0; i < NumPts; ++i)  { free(points_p[i]); }    free(points_p);
    if (modp == nullptr)  { FFdtor(FFp); return std::vector<RingElem>(); } // empty list means error

    // Now extract the answer...
    const int GBsize = modp->GBsize;
    std::vector<RingElem> GB(GBsize);
    vector<long> expv(NumVars);
    for (int i=0; i < GBsize; ++i)
    {
      CheckForInterrupt("BM_projective_mod_p");
      for (int var = 0; var < NumVars; ++var)
        expv[var] = modp->pp[modp->GB[i]][var];
      RingElem GBelem = monomial(P, expv);
      //      std::cout << expv << std::endl;
      long d = deg(NumVars, modp->pp[modp->GB[i]]);
      long top;
      for (top=modp->npp-1;top>0;top--)
        if (deg(NumVars,modp->pp[top]) == d)  break;
      for (int j=0; j < NumPts; ++j)
      {
        const int c = modp->M[modp->GB[i]][j+NumPts];
        if (c == 0) continue;
        long n;
        for (n=top; n >= 0; n--)  if (modp->issep[n] == j)  break;
        for (int var = 0; var < NumVars; ++var)
          expv[var] = modp->pp[n][var];
        GBelem += monomial(P, c, expv);
        //        std::cout << expv << std::endl;
        }
      GB[i] = GBelem;
      //      std::cout << "GB[i] = " << GB[i] << std::endl;
    }
    BM_dtor(modp);
    FFdtor(FFp);
    return GB;
  }
  
  std::vector<RingElem> BM_QQ(const SparsePolyRing& P, const ConstMatrixView& pts_in)
  {
    const long NumPts = NumRows(pts_in);
    const long dim = NumCols(pts_in);
    matrix pts = NewDenseMat(RingQQ(), NumPts, dim);
    for (long i=0; i < NumPts; ++i)
      for (long j=0; j < dim; ++j)
      {
        BigRat q;
        if (!IsRational(q, pts_in(i,j)))  throw 999;
        SetEntry(pts,i,j, q);
      }

    // Ensure input pts have integer coords by using
    // scale factors for each indet.
    vector<BigInt> ScaleFactor(dim, BigInt(1));
    for (long j=0; j < dim; ++j)
      for (long i=0; i < NumPts; ++i)
        ScaleFactor[j] = lcm(ScaleFactor[j], ConvertTo<BigInt>(den(pts(i,j))));

    mpz_t **points = (mpz_t**)malloc(NumPts*sizeof(mpz_t*));
    for (long i=0; i < NumPts; ++i)
    {
      points[i] = (mpz_t*)malloc(dim*sizeof(mpz_t));
      for (long j=0; j < dim; ++j) mpz_init(points[i][j]);
      for (long j=0; j < dim; ++j)
      {
        mpz_set(points[i][j], mpzref(ConvertTo<BigInt>(ScaleFactor[j]*pts(i,j))));
      }
    }


    BMGB char0; // these will be "filled in" by BM_affine below
    BM modp;    //
            
    pp_cmp_PPM = &PPM(P); // not threadsafe!
    BM_affine(&char0, &modp, dim, NumPts, points, pp_cmp); // THIS CALL DOES THE REAL WORK!!!
    pp_cmp_PPM = nullptr;
    for (long i=NumPts-1; i >=0 ; --i)
    {
      for (long j=0; j < dim; ++j)  mpz_clear(points[i][j]);
      free(points[i]);
    }
    free(points);

    if (modp == nullptr)
    {
      if (char0 != nullptr)  BMGB_dtor(char0);
      CoCoA_THROW_ERROR1(ERR::ShouldNeverGetHere);
    }

    // Now extract the answer...
    const int GBsize = char0->GBsize;
    std::vector<RingElem> GB(GBsize);
    const long NumVars = dim;
    vector<long> expv(NumVars); // buffer for creating monomials
    for (int i=0; i < GBsize; ++i)
    {
      CheckForInterrupt("BM_QQ");
      BigInt denom(1); // scale factor needed to make GB elem monic.
      for (int var = 0; var < NumVars; ++var)
      {
        expv[var] = modp->pp[modp->GB[i]][var];
        denom *= power(ScaleFactor[var], expv[var]);
      }
      RingElem GBelem = monomial(P, expv);

      for (int j=0; j < NumPts; ++j)
      {
        if (mpq_sgn(char0->GB[i][j])==0)  continue;
        BigRat c = BigRatFromMPQ(char0->GB[i][j]);
        for (int var = 0; var < NumVars; ++var)
        {
          expv[var] = modp->pp[modp->sep[j]][var];
          c *= power(ScaleFactor[var], expv[var]);
        }
        GBelem += monomial(P, c/denom, expv);
      }
      GB[i] = GBelem;
    }
    BMGB_dtor(char0);
    BM_dtor(modp);
    return GB;
    // ignoring separators for the moment
  }


  std::vector<RingElem> BM_projective_QQ(const SparsePolyRing& P, const ConstMatrixView& pts_in)
  {
    const long NumPts = NumRows(pts_in);
    const long NumVars = NumCols(pts_in);
    matrix pts = NewDenseMat(RingQQ(), NumPts, NumVars);
    for (long i=0; i < NumPts; ++i)
      for (long j=0; j < NumVars; ++j)
      {
        BigRat q;
        if (!IsRational(q, pts_in(i,j)))  throw 999;
        SetEntry(pts,i,j, q);
      }

    // Ensure input pts have integer coords by using
    // scale factors for each indet.
    vector<BigInt> ScaleFactor(NumVars, BigInt(1));
    for (long j=0; j < NumVars; ++j)
      for (long i=0; i < NumPts; ++i)
        ScaleFactor[j] = lcm(ScaleFactor[j], ConvertTo<BigInt>(den(pts(i,j))));

    mpz_t **points = (mpz_t**)malloc(NumPts*sizeof(mpz_t*));
    for (long i=0; i < NumPts; ++i)
    {
      points[i] = (mpz_t*)malloc(NumVars*sizeof(mpz_t));
      for (long j=0; j < NumVars; ++j)  mpz_init(points[i][j]);
      for (long j=0; j < NumVars; ++j)
      {
        mpz_set(points[i][j], mpzref(ConvertTo<BigInt>(ScaleFactor[j]*pts(i,j))));
      }
    }


    BMGB char0; // these will be "filled in" by BM_projective below
    BM modp;    //
            
    pp_cmp_PPM = &PPM(P); // not threadsafe!
    BM_projective(&char0, &modp, NumVars, NumPts, points, pp_cmp); // THIS CALL DOES THE REAL WORK!!!
    pp_cmp_PPM = nullptr;
    for (long i=NumPts-1; i >=0 ; --i)
    {
      for (long j=0; j < NumVars; ++j)  mpz_clear(points[i][j]);
      free(points[i]);
    }
    free(points);

    if (modp == nullptr)
    {
      if (char0 != nullptr)  BMGB_dtor(char0);
      CoCoA_THROW_ERROR1(ERR::ShouldNeverGetHere);
    }

    // Now extract the answer...
    const int GBsize = char0->GBsize;
    std::vector<RingElem> GB(GBsize);
    vector<long> expv(NumVars); // buffer for creating monomials
    for (int i=0; i < GBsize; ++i)
    {
      BigInt denom(1); // scale factor needed to make GB elem monic.
      for (int var = 0; var < NumVars; ++var)
      {
        expv[var] = modp->pp[modp->GB[i]][var];
        denom *= power(ScaleFactor[var], expv[var]);
      }
      RingElem GBelem = monomial(P, expv);
      //std::cout << "* GBelem is " << GBelem << std::endl;

      long d = deg(NumVars, modp->pp[modp->GB[i]]);
      long top;
      for (top=modp->npp-1; top>0; --top)
        if (deg(NumVars, modp->pp[top])==d)  break;
      for (int j=0; j < NumPts; ++j)
      {
        if (mpq_sgn(char0->GB[i][j])==0)  continue;
        BigRat c = BigRatFromMPQ(char0->GB[i][j]);
        long n;
        for (n=top; n>=0; --n)  if (modp->issep[n]==j)  break;
        for (int var = 0; var < NumVars; ++var)
        {
          expv[var] = modp->pp[n][var];
          c *= power(ScaleFactor[var], expv[var]);
        }
        GBelem += monomial(P, c/denom, expv);
        //std::cout << j << "-GBelem is " << GBelem << std::endl;
      }
      GB[i] = GBelem;
    }
    BMGB_dtor(char0);
    BM_dtor(modp);
    return GB;
    // ignoring separators for the moment
  }


  bool DuplicateRows(const ConstMatrixView& M)
  {
    const long nrows = NumRows(M);
    const long ncols = NumCols(M);
    for (long i1=0; i1 < nrows; ++i1)
      for (long i2=i1+1; i2 < nrows; ++i2)
      {
        bool NoDifference = true;
        for (long j=0; NoDifference && j < ncols; ++j)
          NoDifference = (M(i1,j) == M(i2,j));
        if (NoDifference) return true;
      }
    return false;
  }

  bool ProjDuplicateRows(const ConstMatrixView& M)
  {
    const long nrows = NumRows(M);
    const long ncols = NumCols(M);
    for (long i1=0; i1 < nrows; ++i1)
    {
      long FirstNonZeroIdx = -1;
      for (long j=0; (FirstNonZeroIdx==-1) && j < ncols; ++j)
        if (M(i1,j) != 0)  FirstNonZeroIdx = j;
      if (FirstNonZeroIdx == -1)  return true;
      for (long i2=i1+1; i2 < nrows; ++i2)
      {
        bool NoDifference = true;
        for (long j=0; NoDifference && j < FirstNonZeroIdx; ++j)
          NoDifference = IsZero(M(i2,j));
        if (!NoDifference)  break;
        if (IsZero(M(i2,FirstNonZeroIdx)))  break;
        RingElem scale =  M(i2,FirstNonZeroIdx)/M(i1,FirstNonZeroIdx);
        for (long j=FirstNonZeroIdx; NoDifference && j < ncols; ++j)
          NoDifference = (M(i1,j)*scale == M(i2,j));
        if (NoDifference)  return true;
      }
    }
    return false;
  }

  ideal IdealOfPoints(const SparsePolyRing& P, const ConstMatrixView& M)
  {
    if (NumIndets(P) != NumCols(M)) CoCoA_THROW_ERROR1(ERR::BadMatrixSize);
    if (DuplicateRows(M)) CoCoA_THROW_ERROR2(ERR::BadArg, "points must be distinct");
    if (CoeffRing(P) != RingOf(M))
    {
      //      CoCoA_THROW_ERROR(ERR::MixedRings, "IdealOfPoints");
      return IdealOfPoints(P, CanonicalHom(RingOf(M),CoeffRing(P))( M ));
    }
    
    ideal I(P, std::vector<RingElem>(0));
    if (IsFiniteField(CoeffRing(P)))  I = ideal(P, BM_modp(P, M));
    else if (IsQQ(CoeffRing(P)))  I = ideal(P, BM_QQ(P, M));
    // generic case
    else I = ideal(P, BM_generic(P, M));
    SetGBasisAsGens(I);
    const SparsePolyRingBase::IdealImpl* const ptrI = 
      SparsePolyRingBase::IdealImpl::ourGetPtr(I);
    ptrI->myAssignMaximalFlag(NumRows(M)==1);
    ptrI->myAssignPrimaryFlag(NumRows(M)==1);
    ptrI->myAssignPrimeFlag(NumRows(M)==1);
    ptrI->myAssignRadicalFlag(true);
    return I;
  }


  ideal IdealOfProjectivePoints(const SparsePolyRing& P, const ConstMatrixView& M)
  {
    if (CoeffRing(P) != RingOf(M))  CoCoA_THROW_ERROR1(ERR::MixedRings);
    if (NumIndets(P) != NumCols(M))  CoCoA_THROW_ERROR1(ERR::BadMatrixSize);
    if (ProjDuplicateRows(M))
      CoCoA_THROW_ERROR2(ERR::BadArg, "points must be distinct");
    ideal I(P, std::vector<RingElem>(0));
    if (IsFiniteField(CoeffRing(P)))  I = ideal(P, BM_projective_modp(P, M));
    else if (IsQQ(CoeffRing(P)))  I = ideal(P, BM_projective_QQ(P, M));
    // generic case
    else 
      CoCoA_THROW_ERROR2(ERR::NYI, "general case");
      //I = ideal(P, BM_generic(P, M));
    SetGBasisAsGens(I);
    return I;
  }


} // end of namespace CoCoA
