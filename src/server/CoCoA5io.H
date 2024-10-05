#ifndef CoCoA_CoCoA5io_H
#define CoCoA_CoCoA5io_H

//   Copyright (c)  2019  Anna Bigatti,  John Abbott

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

//#include "CoCoA/library.H"

#include "CoCoA/TmpGTypes.H"  // PolyList and VectorList

#include <string>

namespace CoCoA
{
  class BigInt;           // forward declaration -- CoCoA/BigInt.H
  class RingElemAlias; // forward declaration -- CoCoA/ring.H
  typedef const RingElemAlias& ConstRefRingElem; // forward declaration -- CoCoA/ring.H
  class FreeModule;       // forward declaration -- CoCoA/FreeModule.H
  class PPMonoid;         // forward declaration -- CoCoA/PPMonoid.H
  class PPMonoidElem;     // forward declaration -- CoCoA/PPMonoid.H
  class PPOrdering;       // forward declaration -- CoCoA/PPOrdering.H
  class SparsePolyRing;   // forward declaration -- CoCoA/SparsePolyRing.H
  class ideal;            // forward declaration -- CoCoA/ideal.H
  class matrix;           // forward declaration -- CoCoA/matrix.H
  //----------------------------------------------------------------------
  enum SkipTagType {TagWasRead, GetTag};  // used in functions Read...(...)

  //----------------------------------------------------------------------
  std::string ReadOperationString(std::istream& in);
  void ReadPPs(std::istream& in, std::vector<PPMonoidElem>& PPs, const PPMonoid& PPM, SkipTagType ST);
  void ReadPolyList(std::istream& in, PolyList& PL, const SparsePolyRing& P, SkipTagType ST);
  void ReadVectorList(std::istream& in, VectorList& L, const FreeModule& FM, SkipTagType ST);

  // cannot be "void Read..(res, ...)" because there is no empty ctor
  // anyway, they are reference counted objects
  int ReadVerbosityLevel(std::istream& in, SkipTagType ST);
  matrix ReadIntegerMatrix(std::istream& in, SkipTagType ST);
  matrix ReadRationalMatrix(std::istream& in, SkipTagType ST);
  FreeModule ReadFreeModule(std::istream& in, SkipTagType ST);
  SparsePolyRing ReadPolyRing(std::istream& in, SkipTagType ST);
  PPOrdering ReadPolyRingOrdering(std::istream& in, long NumInds, SkipTagType ST);

  //------------------------------
  void PrintTimeToLog(double T);
  void PrintTimeToCoCoA5(double T);
  void PrintVersionToCoCoA5();
  void EndOfTransmissionToCoCoA5();

  void WritePolyRing(std::ostream& out, const SparsePolyRing& P);
  void WriteMatrix(std::ostream& out, const matrix &M);
  void WritePoly(std::ostream& out, ConstRefRingElem f);
  void WritePolyList(std::ostream& out, const PolyList &PL);
  void WriteIdeal(std::ostream& out, const ideal &I);

  void WriteMatrixInVar(std::ostream& out, const std::string& VarName, const matrix &M);
  void WritePolyListInVar(std::ostream& out, const std::string& VarName, const PolyList& PL);
  void WriteVectorListInVar(std::ostream& out, const std::string& VarName, const VectorList& VL);
  //----------------------------------------------------------------------
  void AssertTag(const std::string& input_tag, const std::string& expected_tag);
  void SkipTag(std::istream& in, const std::string& expected_tag);
  void ThrowInputError(const std::string& unknown_tag);

  //----------------------------------------------------------------------
  SparsePolyRing NewPolyRingServer(const BigInt& charact, long FloatPrecision, long NumParams, const PPOrdering& O);
}


// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/server/CoCoA5io.H,v 1.2 2022/02/22 20:40:04 abbott Exp $
// $Log: CoCoA5io.H,v $
// Revision 1.2  2022/02/22 20:40:04  abbott
// Summary: Updated copyright message (redmine 1555)
//
// Revision 1.1  2019/03/20 16:22:15  bigatti
// -- first import
//

#endif
