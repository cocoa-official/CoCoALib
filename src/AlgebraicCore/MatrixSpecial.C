//   Copyright (c)  2010,2016  Anna Bigatti

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


#include "CoCoA/DenseMatrix.H"
#include "CoCoA/MachineInt.H"
#include "CoCoA/BigRat.H"
#include "CoCoA/MatrixSpecial.H"
#include "CoCoA/MatrixOps.H" // for matrix product in RandomUnimodularMat
#include "CoCoA/MatrixView.H" // for BlockMat2x2
#include "CoCoA/PolyRing.H" // for JacobianMat
#include "CoCoA/SparsePolyOps-RingElem.H" // for CoeffVecWRT (in SylvesterMat)
#include "CoCoA/error.H"
#include "CoCoA/interrupt.H"
#include "CoCoA/matrix.H"
#include "CoCoA/random.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/utils.H" // for len 
#include "CoCoA/VectorOps.H" // for HasUniqueOwner

#include <iostream>
using std::ostream;
#include <limits>
using std::numeric_limits;
//#include <vector>
using std::vector;

namespace CoCoA
{


  //---------  jacobian matrix

  matrix JacobianMat_aux(const PolyRing& P, const std::vector<RingElem>& polys, const std::vector<RingElem>& inds)
  {
    const long LenPolys = len(polys);
    const long LenInds = len(inds);
    matrix ans = NewDenseMat(P, LenPolys, LenInds);
    for (long i=0; i < LenPolys; ++i)
      for (long j=0; j < LenInds; ++j)
        SetEntry(ans, i,j, deriv(polys[i], inds[j]));
    return ans;
  }


  matrix JacobianMat(const std::vector<RingElem>& polys)
  {
    if (polys.empty())  CoCoA_THROW_ERROR1(ERR::ReqNonEmpty);
    if (!HasUniqueOwner(polys))  CoCoA_THROW_ERROR1(ERR::MixedRings);
    const PolyRing P = owner(polys[0]);
    return JacobianMat_aux(P, polys, indets(P));
  }
  

  matrix JacobianMat(const std::vector<RingElem>& polys, const std::vector<RingElem>& inds)
  {
    if (len(inds)==0 && len(polys)==0) CoCoA_THROW_ERROR1(ERR::ReqNonEmpty);
    if (!HasUniqueOwner(inds))  CoCoA_THROW_ERROR1(ERR::MixedRings);
    if (!HasUniqueOwner(polys))  CoCoA_THROW_ERROR1(ERR::MixedRings);
    if (len(inds)==0)
      return NewDenseMat(owner(polys[0]), len(polys), 0);
    const PolyRing P(owner(inds[0]));
    if (polys.empty())  return NewDenseMat(P, 0, len(inds));
    if (owner(polys[0]) != P)  CoCoA_THROW_ERROR1(ERR::MixedRings);
    return JacobianMat_aux(P, polys, inds);
  }
  


  matrix LawrenceMat(ConstMatrixView M)
  {
    const long nrows = NumRows(M);
    const long ncols = NumCols(M);
    const ring& R = RingOf(M);
    ConstMatrix Z = ZeroMat(R, nrows, ncols);
    ConstMatrix I = IdentityMat(R, ncols);
    return NewDenseMat(BlockMat2x2(M, Z,
                                   I, I));
  }
  


  //---------  Sylvester matrix
  // Creates a new DenseMat; could also make a ConstMatrixView version...(KISS?)

  matrix SylvesterMat(ConstRefRingElem f, ConstRefRingElem g, ConstRefRingElem x)
  {
    const PolyRing P = owner(x);
    if (owner(f) != P || owner(g) != P)
      CoCoA_THROW_ERROR1(ERR::MixedRings);
    if (IsZero(f) || IsZero(g))
      CoCoA_THROW_ERROR1(ERR::ReqNonZeroRingElem);
    long index;
    if (!IsIndet(index, x))
      CoCoA_THROW_ERROR1(ERR::ReqIndet);

    const long degf = deg(f, index);
    const long degg = deg(g, index);
    const vector<RingElem> cf = CoeffVecWRT(f,x);
    const vector<RingElem> cg = CoeffVecWRT(g,x);

    matrix ans = NewDenseMat(P, degf+degg, degf+degg);
    for (int i=0; i < degg; ++i)
      for (int j=0; j <= degf; ++j)
        SetEntry(ans,i,i+j, cf[degf-j]);

    for (int i=0; i < degf; ++i)
      for (int j=0; j <= degg; ++j)
        SetEntry(ans,i+degg,i+j, cg[degg-j]);

    return ans;
  }
  

// VandermondeMatrix
// HessianMatrix
// HilbertMatrix
// HilbertInverseMatrix
// ToeplitzMatrix
// WronskianMatrix


  matrix RandomUnimodularMat(const ring& R, const MachineInt& N, const MachineInt& Niters /*=0*/)
  {
    if (IsNegative(N) || !IsSignedLong(N))
      CoCoA_THROW_ERROR2(ERR::ReqNonNegative, "dimension");
    if (IsNegative(Niters) || !IsSignedLong(Niters))
      CoCoA_THROW_ERROR2(ERR::ReqNonNegative, "num iters");
    if (IsZero(N)) return NewDenseMat(R,0,0);
    const long n = AsSignedLong(N);
    if (n == 1)  return NewDenseMat(IdentityMat(R,1)); 
    const int niters = IntegerCast<int>(IsZero(Niters)?25*n:AsSignedLong(Niters));
    if (niters > 250*n)
      return RandomUnimodularMat(R, N, niters/2)*RandomUnimodularMat(R, N, niters-niters/2); // Sometimes a SLUG!!!
    vector< vector<BigInt> > VV(n, vector<BigInt>(n));
    for (int i=0; i < n; ++i)
      VV[i][i] = (RandomBool())?1:-1;

    for (int iter=0; iter < niters; ++iter)
    {
      CheckForInterrupt("RandomUnimodularMat: main loop");
      const int i = RandomLong(0, n-1);
      int j = RandomLong(0, n-2);
      if (j == i)  j = n-1;
      // Decide whether to add or subtract row j to row i...
      if (RandomBool())
      {
        for (int k=0; k < n; ++k)
          VV[i][k] += VV[j][k];
      }
      else
      {
        for (int k=0; k < n; ++k)
          VV[i][k] -= VV[j][k];
      }
    }
    matrix ans = NewDenseMat(R, n,n);
    for (int i=0; i < n; ++i)
      for (int j=0; j < n; ++j)
        SetEntry(ans,i,j, VV[i][j]);
    return ans;
  }


  // >>>INEFFICIENT<<<
  // Simple impl, but would be better to use a ConstMatrix (like IdentityMat)
  matrix HilbertMat(const MachineInt& N)
  {
    if (IsNegative(N) || !IsSignedLong(N))
      CoCoA_THROW_ERROR1(ERR::ReqNonNegative);
    const long n = AsSignedLong(N);

    matrix ans = NewDenseMat(RingQQ(), n,n);
    for (int i=0; i < n; ++i)
      for (int j=0; j < n; ++j)
        SetEntry(ans,i,j, BigRat(1, i+j+1));
    return ans;
  }




} // end of namespace CoCoA
