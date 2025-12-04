#ifndef CoCoA5_CoCoALibSupplement_H
#define CoCoA5_CoCoALibSupplement_H

//   Copyright (c)  2011-2017,2020,2022  John Abbott,  Anna M. Bigatti

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

#include <vector>
#include "CoCoA/factorization.H"
#include "CoCoA/error.H"

namespace CoCoA
{
// fwd decl -- defined in library.H
	class BigInt;
	class BigRat;
  struct CoeffPP;
	class ConstMatrixView;
	class RingElemAlias;
  typedef const RingElemAlias& ConstRefRingElem;
	class ConstMatrix;
	class MatrixView;
	class RingElem;
	class RingHom;
	class degree;
	class ideal;
	class matrix;
	class FinGenModule;
	class module;
	class ModuleElem;
	class ring;
        class SparsePolyRing;

  std::vector<long> VectorLong(const std::vector<BigInt>& BigIntVec, const std::string& FuncName);
  std::vector<long> VectorLongDecr1(const std::vector<BigInt>& BigIntVec, const ERR::ID& ErrID, const std::string& FuncName);

  std::string fold_forC5(const std::string& str, const BigInt& MaxLineLen);

  std::vector<BigInt> DegreeToVec(degree d);  // for Value::from

  RingElem LPP_forC5(ConstRefRingElem f); ///< NB result belongs to owner(f)
  RingElem LPP_forC5(const ModuleElem& f);
  RingElem LM_forC5(ConstRefRingElem f);  ///< NB result belongs to owner(f)
  ModuleElem LM_forC5(const ModuleElem& f);
  bool IsCoprime_forC5(ConstRefRingElem f, ConstRefRingElem g);
  bool IsTerm_forC5(ConstRefRingElem f);
  long UnivariateIndetIndex_forC5(ConstRefRingElem f);
  RingElem LT_forC5(ConstRefRingElem f); ///< NB result belongs to owner(f)
  ModuleElem LT_forC5(const ModuleElem& f);
  long LPosn_forC5(const ModuleElem& f);
  long FirstNonZeroPosn_forC5(const ModuleElem& f);
  RingElem CoeffOfTerm_forC5(ConstRefRingElem f, ConstRefRingElem t); ///< NB result belongs to CoeffRing(owner(f))
  std::vector<RingElem> coefficients_forC5(ConstRefRingElem f);  ///< NB result belongs to CoeffRing(owner(f))
  std::vector<RingElem> monomials_forC5(ConstRefRingElem f);  ///< NB result belongs to owner(f)
  std::vector<RingElem> support_forC5(ConstRefRingElem f);  ///< NB result belongs to owner(f)
  bool IsFactorClosed_forC5(const std::vector<RingElem>& v);

  RingElem DensePoly_forC5(const ring& P, const BigInt& D);

  BigInt NextPrime_forC5(const BigInt& N);
  BigInt PrevPrime_forC5(const BigInt& N);
  factorization<BigInt> FactorINT_forC5(const BigInt& N);
// 2023-12-22  factorization<BigInt> FactorINT_timeout_forC5(const BigInt& N, const BigInt& TimeLimit);
  factorization<BigInt> FactorINT_TrialDiv_forC5(const BigInt& N, const BigInt& TrialLimit);
  factorization<BigInt> FactorINT_PollardRho_forC5(const BigInt& N, const BigInt& Niters);
  std::vector<long> AllFactors_forC5(const BigInt& N);
  BigRat CpuTime_forC5();
  BigRat ElapsedTime_forC5();
  //  std::string date_forC5();
  BigInt date_forC5();
  BigInt TimeOfDay_forC5();
  void reseed_forC5(const BigInt& seed);
  const BigInt random_forC5(const BigInt& lo, const BigInt& hi);
  BigInt CertifyNotPower_forC5(BigInt N, BigInt UPBexp);
  std::vector<long> InvTotient_forC5(const BigInt& N); // N must fit into long!!

  BigInt lcm_forC5(const std::vector<BigInt>& v);
  BigInt gcd_forC5(const std::vector<BigInt>& v);
  std::vector<BigInt> CoprimeFactorBasis_forC5(const std::vector<BigInt>& v);
  std::vector<RingElem> CoprimeFactorBasis_forC5(const std::vector<RingElem>& v);
  RingElem lcm_forC5(const std::vector<RingElem>& v);
  RingElem gcd_forC5(const std::vector<RingElem>& v);
  RingElem ContentWRT_forC5(ConstRefRingElem f, const std::vector<RingElem>& v);
  std::vector<CoeffPP> CoefficientsWRT_forC5(ConstRefRingElem f, const std::vector<RingElem>& v);
  std::vector<RingElem> CoeffListWRT_forC5(ConstRefRingElem f, ConstRefRingElem x);
  std::vector<RingElem> CoeffListWRTSupport_forC5(ConstRefRingElem f, ConstRefRingElem basis);
  RingElem HomogCompt_forC5(ConstRefRingElem f, const BigInt& D);
///???  RingElem resultant_forC5(ConstRefRingElem f, ConstRefRingElem g);
  RingElem resultant_forC5(ConstRefRingElem f, ConstRefRingElem g, ConstRefRingElem x);
///???  RingElem discriminant_forC5(ConstRefRingElem f);
  RingElem discriminant_forC5(ConstRefRingElem f, ConstRefRingElem x);

  std::vector<RingElem> QuotientBasis_forC5(const ideal& I);
  std::vector<RingElem> QuotientBasisSorted_forC5(const ideal& I);
  std::vector<RingElem> JanetBasis_forC5(const ideal& I);
  RingElem HilbertNumQuot_forC5(const ideal& I);
  std::vector<BigInt> ContFrac_forC5(const BigRat& q);
  std::vector<BigRat> CFApproximants_forC5(const BigRat& q);
  std::vector<BigInt> BinomialRepr_forC5(const BigInt& N, const BigInt& r);
  BigInt BinomialReprShift_forC5(const BigInt& N, const BigInt& r, const BigInt& shift1, const BigInt& shift2);
  BigInt NumPartitions_forC5(const BigInt& N);
  BigInt CatalanNumber_forC5(const BigInt& N);

  ring RingQQt_forC5(const BigInt& NumIndets);

  //----- matrix --------------------------------------------------
  matrix NewDenseMat_forC5(const ring& R, const BigInt& NR, const BigInt& NC);
//   matrix HomogElimMat_forC5(ConstMatrixView M, const std::vector<BigInt>& ElimInds);
//   matrix ElimMat_forC5(ConstMatrixView M, const std::vector<BigInt>& ElimInds);
//   matrix ElimMat_forC5(const BigInt& n, const std::vector<BigInt>& ElimInds);
  matrix LexMat_forC5(const BigInt& N);
  matrix RevLexMat_forC5(const BigInt& N);
  matrix StdDegLexMat_forC5(const BigInt& N);
  matrix StdDegRevLexMat_forC5(const BigInt& N);
  matrix XelMat_forC5(const BigInt& N);

  void SetEntry_forC5(MatrixView& M, const BigInt& I, const BigInt& J, ConstRefRingElem x);
  std::vector<RingElem> GetRow_forC5(ConstMatrixView M, const BigInt& I);
  std::vector<RingElem> GetCol_forC5(ConstMatrixView M, const BigInt& I);
  void SetRow_forC5(MatrixView& M, const BigInt& I, const std::vector<RingElem>& v);
  void SetCol_forC5(MatrixView& M, const BigInt& I, const std::vector<RingElem>& v);
  void SwapRows_forC5(matrix& M, const BigInt& row1, const BigInt& row2);
  void SwapCols_forC5(matrix& M, const BigInt& col1, const BigInt& col2);
  ConstMatrix IdentityMat_forC5(const ring& R, const BigInt& N);
  ConstMatrix ZeroMat_forC5(const ring& R, const BigInt& NRows, const BigInt& NCols);
  ConstMatrixView transposed_forC5(ConstMatrixView M);

  //---- points ----
  void ApproxPointsNBM_forC5(std::vector<RingElem>& QB, std::vector<RingElem>& BB, std::vector<RingElem>& AV, const SparsePolyRing& P, ConstMatrixView OrigPts, ConstMatrixView OrigTolerance);
  void ApproxPointsSOI_forC5(std::vector<RingElem>& QB, std::vector<RingElem>& BB, std::vector<RingElem>& AV, const SparsePolyRing& P, ConstMatrixView OrigPts, ConstMatrixView OrigTolerance);
  RingElem ClosePassingPoly_forC5(const ring& P, ConstMatrixView OrigPts, ConstMatrixView OrigTolerance);
  void PreprocessPts_forC5(const std::string& WhichAlgm, std::vector< std::vector<RingElem> >& NewPts, std::vector<long>& weights, ConstMatrixView OrigPts, ConstMatrixView epsilon);

  //---- module ----
  module NewFreeModule_forC5(const ring& R, ConstMatrixView M);
  long NumCompts_forC5(const module& M);
  ModuleElem NewFreeModuleElem(const module& M, const std::vector<RingElem>& v);
  const std::vector<ModuleElem>& CanonicalBasis_forC5(const module& F);
  FinGenModule SubmoduleCols_forC5(const module& F, ConstMatrixView M);
  FinGenModule SubmoduleRows_forC5(const module& F, ConstMatrixView M);

  // Delete me if you want to (by Nico Mexis)
  RingElem GraeffeN_forC5(ConstRefRingElem f, const BigInt& N);
  RingElem cyclotomic_forC5(const BigInt& N, ConstRefRingElem x); // x must be indet
  BigInt CyclotomicIndex_forC5(ConstRefRingElem f);
  BigInt CyclotomicTest_forC5(ConstRefRingElem f);
  BigInt LRSDegeneracyOrder_forC5(ConstRefRingElem f);
  std::vector<BigInt> LRSDegeneracyOrders_forC5(ConstRefRingElem f);
  std::vector<BigInt> LRSDegeneracyOrdersUnchecked_forC5(ConstRefRingElem f, const BigInt& VerifLevel);
  bool IsLRSDegenerate_forC5(ConstRefRingElem f);
  bool IsLRSDegenerateOrder_forC5(ConstRefRingElem f, const BigInt& N);
  BigInt MoebiusFn_forC5(const BigInt &N);

// Moved to Interpreter.C  (used only there)
/////  //---- utils ----
/////  std::string message_forC5(const exception& exc);

}

// Close include guard:
#endif
