#ifndef CoCoA_SparsePolyOps_eval_H
#define CoCoA_SparsePolyOps_eval_H

//   Copyright (c)  2023  John Abbott,  Anna M. Bigatti

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

#include "CoCoA/BigInt.H"
#include "CoCoA/ring.H"

#include <iosfwd>
//using std::ostream;
#include <vector>

namespace CoCoA
{

  // Functions to evaluate a UNIVARIATE POLYNOMIAL over ZZ at an integer/rational point
  // !! If the evaluation point is rational then only the numerator is returned !!
  
  // Use EvalUPoly if you plan to evaluate f more than once
  class EvalUPoly
  {
  public:
    EvalUPoly(ConstRefRingElem f);
    BigInt operator()(long a) const;             // returns f(a)
    BigInt operator()(const BigInt& a) const;    // returns f(a)
    BigInt operator()(long n, long d) const;                    // returns "NUMERATOR", i.e. d^deg(f)*f(n/d)
    BigInt operator()(const BigInt& n, const BigInt& d) const;  // returns "NUMERATOR", i.e. d^deg(f)*f(n/d)
    std::ostream& myOutput(std::ostream& out) const;
  private: // data members
    std::vector<BigInt> myCoeffs;
  };


  // Evaluate at integer
  BigInt EvalAt(ConstRefRingElem f, long a);
  BigInt EvalAt(ConstRefRingElem f, const BigInt& a);

  // These return just the NUMERATOR:  namely d^deg(f) * f(n/d)
  BigInt EvalAt(ConstRefRingElem f, long n, long d);
  BigInt EvalAt(ConstRefRingElem f, const BigInt& n, const BigInt& d);
    

  std::ostream& operator<<(std::ostream& out, const EvalUPoly& eval);

} // end of namespace CoCoA

#endif
