//   Copyright (c)  2017  John Abbott,  Anna M. Bigatti

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

#include "CoCoA/RandomSparseNonSing01Mat.H"

#include "CoCoA/MachineInt.H"
#include "CoCoA/error.H"
#include "CoCoA/ring.H"
#include "CoCoA/DenseMatrix.H"
#include "CoCoA/interrupt.H"
#include "CoCoA/matrix.H"
#include "CoCoA/MatrixOps.H"
#include "CoCoA/random.H"
#include "CoCoA/verbose.H"

#include <cmath>
using std::log;
#include <iostream>
using std::endl;
#include <vector>
using std::vector;

namespace CoCoA
{

  // Create NxN matrix with 0-1 entries.  Each entry is IID with
  // probability "prob" of being 1.  Value of "prob" is chosen so
  // that we expect to generate only "a few" matrices before getting
  // a non-singular one.
  // MagicFactor was determined empirically: lower values give sparser
  // matrices, but we need more iterations before getting a non-sing one.

  // HINT: create candidate in vector<vector<unsigned char>> then compute
  // its det modulo some "largish" prime; if non-zero return, o/w try again
  // slight risk of discarding a mat whose det is non-zero but div by
  // test prime...  Maybe change prime each time det is computed.

  matrix RandomSparseNonSing01Mat(const ring& R, const MachineInt& N)
  {
    if (IsNegative(N) || !IsSignedLong(N))
      CoCoA_THROW_ERROR(ERR::ReqPositive, "RandomSparseNonSing01Mat");
    const long n = AsSignedLong(N);
    VerboseLog VERBOSE("RandomSparseNonSing01Mat");
    
    constexpr double MagicFactor = 0.80; // empirically determined (see note above)
    double prob = MagicFactor*(std::log(n)/n);
    matrix M = NewDenseMat(R, n,n);
    long NumIters = 0;
    while (true)
    {
      ++NumIters;
      if ((NumIters&7) == 0) prob *= 1.0625; // taking too long, so increase prob
      int NumNZCols = 0;
      vector<bool> ColHasNZEntry(n);
      for (int i=0; i < n; ++i)
      {
        bool RowHasNZEntry = false;
        do
        {
          CheckForInterrupt("RandomSparseNonSing01Mat");
          for (int j=0; j < n; ++j)
          {
            if (RandomBiasedBool(prob))
            {
              if (!ColHasNZEntry[j]) { ColHasNZEntry[j] = true; ++NumNZCols; }
              RowHasNZEntry = true;
              SetEntry(M, i,j, one(R));
            }
          }
        } while (!RowHasNZEntry);
      }
      VERBOSE(80) << "NumNZCols = " << NumNZCols << " out of " << n << endl;
      if (NumNZCols == n && !IsZeroDet(M))
      {
        VERBOSE(70) << "Exit after " << NumIters << " iters  (with prob=" << prob << ")" << std::endl;
        return M;
      }
      // Set all entries of M to zero, and try again.
      for (int i=0; i < n; ++i)
        for (int j=0; j < n; ++j)
          if (!IsZero(M(i,j))) SetEntry(M, i,j, zero(R));
    }
  }


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/RandomSparseNonSing01Mat.C,v 1.9 2024/07/02 15:36:14 bigatti Exp $
// $Log: RandomSparseNonSing01Mat.C,v $
// Revision 1.9  2024/07/02 15:36:14  bigatti
// Summary:Changed error codes: LogZero into ReqNonZero
// and Not... into ReqNonNegative, ReqNonNegativeGrading, ReqPositive, ReqPositiveGrading
//
// Revision 1.8  2024/03/08 19:59:59  abbott
// Summary: Used constexpr (redmine 1511)
//
// Revision 1.7  2022/12/05 21:45:41  abbott
// Summary: Now using IsZeroDet (instead of IsZero(det(...)))
//
// Revision 1.6  2022/02/18 14:11:57  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.5  2020/06/17 15:49:25  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.4  2019/09/23 08:11:26  abbott
// Summary: Now increases prob slowly (o/w takes too long for large matrices)
//
// Revision 1.3  2019/03/04 16:17:17  abbott
// Summary: Made faster; added verbose output (level 70)
//
// Revision 1.2  2018/05/17 15:40:20  bigatti
// -- renamed MatrixOperations --> MatrixOps
//
// Revision 1.1  2017/11/14 14:59:35  abbott
// Summary: New fns: RandomSmallPrime, RandomSparseNonSing01Matrix
//
//
