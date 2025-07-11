//   Copyright (c)  2016-2018  John Abbott,  Anna M. Bigatti

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


#include "CoCoA/library.H"

#include <string>
#include<iostream>
using std::endl;

namespace CoCoA
{

  namespace // anonymous for file local fn
  {
    // Procedure to give error if obsolescent fns are forbidden, and otherwise print out a warning on CoCoA::LogStream.
    void LogObsolescentFn(const char* const FnName, const char* const UsefulAdvice)
    {
      if (!IsAllowedObsolescentFnCall())
        CoCoA_THROW_ERROR2(ERR::OBSOLESCENT, FnName + std::string(" -- ") + UsefulAdvice);
      LogStream() << "WARNING: called obsolescent fn `" << FnName << "' -- " << UsefulAdvice << endl;
    }

  } // end of namespace anonymous

  
  ///////////////////////////////////////////////////////
  // The obsolescent fns below are ALWAYS defined.
  
  std::size_t SizeInBase(const BigInt& N, long base)  // REMOVED 2022-08-08
  {
    if (base < 2 || base > 36)
      CoCoA_THROW_ERROR1(ERR::BadNumBase);
    return mpz_sizeinbase(mpzref(N), base);
  }




  bool IsRadical(ConstRefPPMonoidElem pp)  // RENAMED to IsSqFree
  {
    LogObsolescentFn("IsRadical(ConstRefPPMonoidElem)", "renamed to IsSqFree");
    return IsSqFree(pp);
  }


  bool AreGensSquareFreeMonomial(const ideal& I)
  {
    LogObsolescentFn("AreGensSquareFreeMonomial(ideal)",
                     "renamed to AreGensSqFreeMonomial");
    return AreGensSqFreeMonomial(I);
  }

  
  PPOrdering NewLexOrdering(const MachineInt& NumIndets)
  {
    LogObsolescentFn("NewLexOrdering", "use pseudo-ctor `lex'");
    return lex(NumIndets);
  }

  PPOrdering NewStdDegLexOrdering(const MachineInt& NumIndets)
  {
    LogObsolescentFn("NewStdDegLexOrdering", "use pseudo-ctor `StdDegLex'");
    return StdDegLex(NumIndets);
  }
  
  PPOrdering NewStdDegRevLexOrdering(const MachineInt& NumIndets)
  {
    LogObsolescentFn("NewStdDegRevLexOrdering", "use pseudo-ctor `StdDegRevLex'");
    return StdDegRevLex(NumIndets);
  }

  ideal minimalize(const ideal& I)
  {
    LogObsolescentFn("minimalize", "use \"IdealOfMinGens\"");
    return IdealOfMinGens(I);
  }

  FGModule minimalize(const FGModule& M)
  {
    LogObsolescentFn("minimalize", "use \"SubmoduleOfMinGens\"");
    return SubmoduleOfMinGens(M);
  }


  SparsePolyRing NewPolyRing(const ring& CoeffRing, long NumIndets, const PPOrdering& ord)
  {
    LogObsolescentFn("NewPolyRing (without symbol names)", "use \"SymbolRange\"");
    return NewPolyRing(CoeffRing, SymbolRange("x", 0, NumIndets-1), ord);
    // if (IsRingFp(CoeffRing))
    //   return NewPolyRing_DMPII(CoeffRing, NumIndets, OrdCtor);
    // return NewPolyRing_DMPI(CoeffRing, NumIndets, OrdCtor);
  }

  SparsePolyRing NewPolyRing(const ring& CoeffRing, long NumIndets, const PPOrderingCtor& OrdCtor)
  {
    LogObsolescentFn("NewPolyRing (without symbol names)", "use \"SymbolRange\"");
    return NewPolyRing(CoeffRing, SymbolRange("x", 0, NumIndets-1), OrdCtor);
    // if (IsRingFp(CoeffRing))
    //   return NewPolyRing_DMPII(CoeffRing, NumIndets, OrdCtor);
    // return NewPolyRing_DMPI(CoeffRing, NumIndets, OrdCtor);
  }

  SparsePolyRing NewPolyRing(const ring& CoeffRing, long NumIndets)
  {
    LogObsolescentFn("NewPolyRing (without symbol names)", "use \"SymbolRange\"");
    return NewPolyRing(CoeffRing, SymbolRange("x", 0, NumIndets-1), StdDegRevLex);
  }


  BigInt iroot(const MachineInt& n, const MachineInt& r)
  {
    LogObsolescentFn("iroot", "use \"FloorRoot\"");
    if (IsNegative(n))  CoCoA_THROW_ERROR2(ERR::ReqNonNegative, "1st arg");
    return FloorRoot(n,r);
  }
  
  BigInt iroot(const MachineInt& n, const BigInt& R)
  {
    LogObsolescentFn("iroot", "use \"FloorRoot\"");
    if (IsNegative(n))  CoCoA_THROW_ERROR2(ERR::ReqNonNegative, "1st arg");
    return FloorRoot(n,R);
  }
    
  BigInt iroot(const BigInt& N,     const MachineInt& r)
  {
    LogObsolescentFn("iroot", "use \"FloorRoot\"");
    if (N < 0)  CoCoA_THROW_ERROR1(ERR::ReqNonNegative);
    return FloorRoot(N,r);
  }

  BigInt iroot(const BigInt& N,     const BigInt& R)
  {
    LogObsolescentFn("iroot", "use \"FloorRoot\"");
    if (N < 0)  CoCoA_THROW_ERROR1(ERR::ReqNonNegative);
    return FloorRoot(N,R);
  }


  factorization<long>   SmoothFactor(const MachineInt& N, const MachineInt& TrialLimit)
  {
    LogObsolescentFn("SmoothFactor", "use \"factor_TrialDiv\"");
    return factor_TrialDiv(N,TrialLimit);
  }

  factorization<BigInt> SmoothFactor(const BigInt& N,     const MachineInt& TrialLimit)
  {
    LogObsolescentFn("SmoothFactor", "use \"factor_TrialDiv\"");
    return factor_TrialDiv(N,TrialLimit);
  }

  factorization<BigInt> SmoothFactor(const BigInt& N,     const BigInt& TrialLimit)
  {
    LogObsolescentFn("SmoothFactor", "use \"factor_TrialDiv\"");
    return factor_TrialDiv(N,TrialLimit);
  }




  matrix jacobian(const std::vector<RingElem>& polys)
  {
    LogObsolescentFn("jacobian", "use \"JacobianMat\"");
    return JacobianMat(polys);
  }
  
  matrix jacobian(const std::vector<RingElem>& polys, const std::vector<RingElem>& inds)
  {
    LogObsolescentFn("jacobian", "use \"JacobianMat\"");
    return JacobianMat(polys,inds);
  }

  matrix TensorMat(ConstMatrixView A, ConstMatrixView B)
  {
    LogObsolescentFn("TensorMat", "use \"KroneckerProd\"");
    return KroneckerProd(A,B);
  }


  // 2022
  
  matrix MakeTermOrd(ConstMatrixView M)
  {
    LogObsolescentFn("MakeTermOrd", "use \"MakeTermOrdMat\"");
    return MakeTermOrdMat(M);
  }

  
  matrix MakeTermOrd(ConstMatrixView M, const MachineInt& GrDim)
  {
    LogObsolescentFn("MakeTermOrd", "use \"MakeTermOrdMat\"");
    return MakeTermOrdMat(M, GrDim);
  }


} // end of namespace CoCoA
