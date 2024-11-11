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
    if (polys.empty()) CoCoA_THROW_ERROR1(ERR::ReqNonEmpty);
    if (!HasUniqueOwner(polys)) CoCoA_THROW_ERROR(ERR::MixedRings, "JacobianMat(polys)");
    const PolyRing P = owner(polys[0]);
    return JacobianMat_aux(P, polys, indets(P));
  }
  

  matrix JacobianMat(const std::vector<RingElem>& polys, const std::vector<RingElem>& inds)
  {
    if (len(inds)==0 && len(polys)==0) CoCoA_THROW_ERROR1(ERR::ReqNonEmpty);
    if (!HasUniqueOwner(inds)) CoCoA_THROW_ERROR(ERR::MixedRings, "JacobianMat");
    if (!HasUniqueOwner(polys)) CoCoA_THROW_ERROR(ERR::MixedRings, "JacobianMat");
    if (len(inds)==0)
      return NewDenseMat(owner(polys[0]), len(polys), 0);
    const PolyRing P(owner(inds[0]));
    if (polys.empty()) return NewDenseMat(P, 0, len(inds));
    if (owner(polys[0]) != P) CoCoA_THROW_ERROR(ERR::MixedRings, "JacobianMat");
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
      CoCoA_THROW_ERROR(ERR::MixedRings, "SylvesterMat(f, g, x)");
    if (IsZero(f) || IsZero(g))
      CoCoA_THROW_ERROR(ERR::ReqNonZeroRingElem, "SylvesterMat(f, g, x)");
    long index;
    if (!IsIndet(index, x))
      CoCoA_THROW_ERROR(ERR::ReqIndet, "SylvesterMat(f, g, x)");

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
      CoCoA_THROW_ERROR(ERR::ReqNonNegative, "RandomUnimodularMat, 2nd arg");
    if (IsNegative(Niters) || !IsSignedLong(Niters))
      CoCoA_THROW_ERROR(ERR::ReqNonNegative, "RandomUnimodularMat, 3rd arg");
    if (IsZero(N)) return NewDenseMat(R,0,0);
    const long n = AsSignedLong(N);
    if (n == 1) return NewDenseMat(IdentityMat(R,1)); 
    const int niters = IntegerCast<int>(IsZero(Niters)?25*n:AsSignedLong(Niters));
    if (niters > 250*n) return RandomUnimodularMat(R, N, niters/2)*RandomUnimodularMat(R, N, niters-niters/2); // Sometimes a SLUG!!!
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
      CoCoA_THROW_ERROR(ERR::ReqNonNegative, "HilbertMat");
    const long n = AsSignedLong(N);

    matrix ans = NewDenseMat(RingQQ(), n,n);
    for (int i=0; i < n; ++i)
      for (int j=0; j < n; ++j)
        SetEntry(ans,i,j, BigRat(1, i+j+1));
    return ans;
  }




} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/MatrixSpecial.C,v 1.26 2024/07/15 17:03:57 bigatti Exp $
// $Log: MatrixSpecial.C,v $
// Revision 1.26  2024/07/15 17:03:57  bigatti
// Summary: Error codes:  added CoCoA_THROW_ERROR1(ErrID);
// Removed BadRowIndex, BadColIndex;  Renamed  ERR::Empty --> ERR::ReqNonEmpty
//
// Revision 1.25  2024/07/02 20:42:07  abbott
// Summary: Renaming of errors (redmine 308)
//
// Revision 1.24  2024/07/02 15:31:12  bigatti
// Summary: Changed error codes: LogZero into ReqNonZero
//
// and Not... into ReqNonNegative, ReqNonNegativeGrading, ReqPositive, ReqPositiveGrading
//
// Revision 1.23  2024/07/02 09:57:57  bigatti
// Summary: error codes, first batch:
// ReqUnivariate, ReqNonZero, ReqNonZeroGradingDim, ReqNonZeroModulus,  ReqNonZeroRingElem
//
// Revision 1.22  2022/02/18 14:11:55  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.21  2020/06/17 15:49:24  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.20  2020/05/26 12:06:18  abbott
// Summary: Renamed TensorMat to KroneckerProd; doc & tests updated
//
// Revision 1.19  2020/03/06 18:39:41  abbott
// Summary: Improved readability of LawrenceMat
//
// Revision 1.18  2020/02/27 17:39:45  bigatti
// -- added LawrenceMat
//
// Revision 1.17  2020/01/18 21:40:11  abbott
// Summary: Cleaned SylvesterMat; added recursive call to RandomUnimodularMat
//
// Revision 1.16  2019/10/11 19:54:28  abbott
// Summary: Renamed jacobian to JacobianMat
//
// Revision 1.15  2019/10/11 12:55:05  abbott
// Summary: Added SylvesterMat
//
// Revision 1.14  2018/05/17 15:37:13  bigatti
// -- renamed VectorOperations --> VectorOps
//
// Revision 1.13  2018/01/25 14:04:11  abbott
// Summary: Changed RandomUnimodularMat to produce det +1 or -1
//
// Revision 1.12  2017/11/08 14:04:39  abbott
// Summary: Added new fn HilbertMat
//
// Revision 1.11  2016/10/27 14:03:04  abbott
// Summary: Added RandomUnimodularMat
//
// Revision 1.10  2014/07/31 14:45:17  abbott
// Summary: Merged io.H and UtilsTemplate.H into new header VectorOperations.H
// Author: JAA
//
// Revision 1.9  2014/07/30 14:06:39  abbott
// Summary: Changed BaseRing into RingOf
// Author: JAA
//
// Revision 1.8  2014/07/14 15:06:52  abbott
// Summary: Added include of UtilsTemplate.H
// Author: JAA
//
// Revision 1.7  2014/07/07 12:28:26  abbott
// Summary: Removed AsPolyRing
// Author: JAA
//
// Revision 1.6  2014/04/30 16:08:58  abbott
// Summary: Removed pointless include
// Author: JAA
//
// Revision 1.5  2011/03/23 17:30:12  bigatti
// -- started Sylvester matrix
//
// Revision 1.4  2011/03/21 07:58:29  bigatti
// -- added TensorMat
// -- changed size into len
//
// Revision 1.3  2011/02/22 17:04:41  bigatti
// -- fixed jacobian
//
// Revision 1.2  2011/02/10 15:46:21  bigatti
// -- removed debug printing
//
// Revision 1.1  2011/02/10 15:30:14  bigatti
// -- first import: jacobian
//
