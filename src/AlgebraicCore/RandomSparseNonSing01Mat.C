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
      CoCoA_THROW_ERROR1(ERR::ReqPositive);
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
