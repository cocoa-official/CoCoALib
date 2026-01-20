#ifndef CoCoA_obsolescent_H
#define CoCoA_obsolescent_H

//   Copyright (c)  2016-2018,2022  John Abbott, Anna M. Bigatti

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


// ***** !!! IMPORTANT !!! *****
// This file is not included automatically in CoCoA/library.H;
// you must use a separate #include "CoCoA/obsolescent.H" directive.

// ***** OBSOLESCENT FUNCTIONS *****

#include "CoCoA/factorization.H"

namespace CoCoA
{

  class ideal;                // fwd decl -- defined in ideal.H
  class ConstRefPPMonoidElem; // fwd decl -- defined in PPMonoid.H
  class MachineInt;           // fwd decl -- defined in MachineInt.H
  class BigInt;               // fwd decl -- defined in BigInt.H
  class PPOrdering;           // fwd decl -- defined in PPOrdering.H
  class ring;                 // fwd decl -- defined in ring.H
  class SparsePolyRing;       // fwd decl -- defined in SparsePolyRing.H

  // Added 2016-11-03
  bool IsRadical(ConstRefPPMonoidElem pp);  // RENAMED to IsSqFree

  // Added 2016-11-07
  bool AreGensSquareFreeMonomial(const ideal& I);  // renamed to AreGensSqFreeMonomial

  // Added 2017-11-10
  PPOrdering NewLexOrdering(const MachineInt& NumIndets);          // use pseudo-ctor lex
  PPOrdering NewStdDegLexOrdering(const MachineInt& NumIndets);    // use pseudo-ctor StdDegLex
  PPOrdering NewStdDegRevLexOrdering(const MachineInt& NumIndets); // use pseudo-ctor StdDegRevLex

  // Added 2017-11-20
  ideal minimalize(const ideal& I);       // use IdealOfMinGens
  FGModule minimalize(const FGModule& M); // use SubmoduleOfMinGens
  

  // Added 2018-09-29
  SparsePolyRing NewPolyRing(const ring& CoeffRing, long NumIndets); ///< Default ordering is StdDegRevLex
  SparsePolyRing NewPolyRing(const ring& CoeffRing, long NumIndets, const PPOrdering& ord);
  SparsePolyRing NewPolyRing(const ring& CoeffRing, long NumIndets, const PPOrderingCtor& ord);


  BigInt iroot(const MachineInt& n, const MachineInt& r);
  BigInt iroot(const MachineInt& n, const BigInt& R);
  BigInt iroot(const BigInt& N,     const MachineInt& r);
  BigInt iroot(const BigInt& N,     const BigInt& R);

  factorization<long>   SmoothFactor(const MachineInt& N, const MachineInt& TrialLimit);
  factorization<BigInt> SmoothFactor(const BigInt& N,     const MachineInt& TrialLimit);
  factorization<BigInt> SmoothFactor(const BigInt& N,     const BigInt& TrialLimit);

//   bool IsExactIRoot(long& ans,   const MachineInt& n, const MachineInt& r);
//   bool IsExactIRoot(BigInt& ans, const MachineInt& n, const MachineInt& r);
//   bool IsExactIRoot(long& ans,   const MachineInt& n, const BigInt& R);
//   bool IsExactIRoot(BigInt& ans, const MachineInt& n, const BigInt& R);
// //???  bool IsExactIRoot(long& ans, const BigInt& N, const MachineInt& r);
//   bool IsExactIRoot(BigInt& ans, const BigInt& N,     const MachineInt& r);
// //???  bool IsExactIRoot(long& ans, const BigInt& N, const BigInt& R);
//   bool IsExactIRoot(BigInt& ans, const BigInt& N,     const BigInt& R);


  matrix jacobian(const std::vector<RingElem>& polys);
  matrix jacobian(const std::vector<RingElem>& polys, const std::vector<RingElem>& inds);

  matrix TensorMat(ConstMatrixView A, ConstMatrixView B); // or Kronecker product


  // Added 2022-02 (0.99800)
  
  matrix MakeTermOrd(ConstMatrixView M);
  matrix MakeTermOrd(ConstMatrixView M, const MachineInt& GrDim);

  
  // Added 2022-08 (0.99802)
  std::size_t SizeInBase(const BigInt& N, long base);


// Added 2026-01-20 (0.99880)
  RingElem cyclotomic(long n, ConstRefRingElem x);

} // end of namespace CoCoA

#endif
