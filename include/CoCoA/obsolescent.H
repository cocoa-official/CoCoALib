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


} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/obsolescent.H,v 1.15 2022/08/08 17:53:45 abbott Exp $
// $Log: obsolescent.H,v $
// Revision 1.15  2022/08/08 17:53:45  abbott
// Summary: SizeInBase is obsolete (redmine 1693)
//
// Revision 1.14  2022/03/04 11:40:34  abbott
// Summary: Corrected comment about cocoalib version for obsoleting MakeTermOrd
//
// Revision 1.13  2022/03/04 11:31:28  bigatti
// Summary: added MakeTermOrd
//
// Revision 1.12  2022/02/18 14:11:52  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.11  2022/02/02 09:28:55  abbott
// Summary: Renamed SmoothFactor to factor_TrialDiv (redmine 950)
//
// Revision 1.10  2021/08/04 19:06:58  abbott
// Summary: Removed const (redmine 1606)
//
// Revision 1.9  2020/05/26 12:06:18  abbott
// Summary: Renamed TensorMat to KroneckerProd; doc & tests updated
//
// Revision 1.8  2019/10/11 19:54:28  abbott
// Summary: Renamed jacobian to JacobianMat
//
// Revision 1.7  2019/09/16 17:44:00  abbott
// Summary: Added iroot
//
// Revision 1.6  2018/10/02 09:45:30  abbott
// Summary: Moved pseudo-ctors NewPolyRing(CoeffRing, NumIndets, ...) to obsolescent
//
// Revision 1.5  2017/11/20 20:38:27  bigatti
// -- added minimalized
//
// Revision 1.4  2017/11/10 16:02:27  abbott
// Summary: Removed NewLexOrdering, NewStdDegLexOrdering, NewStdDegRevLexOrdering; consequential changes
//
// Revision 1.3  2016/11/07 14:16:51  bigatti
// -- added AreGensSquareFreeMonomial
//
// Revision 1.2  2016/11/04 20:40:47  abbott
// Summary: Cleaning
//
// Revision 1.1  2016/11/03 12:29:58  abbott
// Summary: Added file for obsolescent fns; also there is a global flag saying whether to give error if calling one.
//
//

#endif
