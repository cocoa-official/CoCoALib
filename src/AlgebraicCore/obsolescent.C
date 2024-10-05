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
        CoCoA_THROW_ERROR(ERR::OBSOLESCENT, FnName + std::string(" -- ") + UsefulAdvice);
      LogStream() << "WARNING: called obsolescent fn `" << FnName << "' -- " << UsefulAdvice << endl;
    }

  } // end of namespace anonymous

  
  ///////////////////////////////////////////////////////
  // The obsolescent fns below are ALWAYS defined.
  
  std::size_t SizeInBase(const BigInt& N, long base)  // REMOVED 2022-08-08
  {
    if (base < 2 || base > 36)
      CoCoA_THROW_ERROR(ERR::BadNumBase, "SizeInBase(BigInt,base)");
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
    if (IsNegative(n)) CoCoA_THROW_ERROR("iroot","1st arg must be non-negative");
    return FloorRoot(n,r);
  }
  
  BigInt iroot(const MachineInt& n, const BigInt& R)
  {
    LogObsolescentFn("iroot", "use \"FloorRoot\"");
    if (IsNegative(n)) CoCoA_THROW_ERROR("iroot","1st arg must be non-negative");
    return FloorRoot(n,R);
  }
    
  BigInt iroot(const BigInt& N,     const MachineInt& r)
  {
    LogObsolescentFn("iroot", "use \"FloorRoot\"");
    if (N < 0) CoCoA_THROW_ERROR("iroot","1st arg must be non-negative");
    return FloorRoot(N,r);
  }

  BigInt iroot(const BigInt& N,     const BigInt& R)
  {
    LogObsolescentFn("iroot", "use \"FloorRoot\"");
    if (N < 0) CoCoA_THROW_ERROR("iroot","1st arg must be non-negative");
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


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/obsolescent.C,v 1.18 2022/08/08 17:53:45 abbott Exp $
// $Log: obsolescent.C,v $
// Revision 1.18  2022/08/08 17:53:45  abbott
// Summary: SizeInBase is obsolete (redmine 1693)
//
// Revision 1.17  2022/02/18 14:12:02  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.16  2022/02/14 14:52:09  bigatti
// Summary: added MakeTermOrd
//
// Revision 1.15  2022/02/02 09:28:55  abbott
// Summary: Renamed SmoothFactor to factor_TrialDiv (redmine 950)
//
// Revision 1.14  2021/08/04 19:09:13  abbott
// Summary: Removed const (redmine 1606)
//
// Revision 1.13  2020/06/17 15:49:30  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.12  2020/05/26 12:06:19  abbott
// Summary: Renamed TensorMat to KroneckerProd; doc & tests updated
//
// Revision 1.11  2019/10/11 19:54:29  abbott
// Summary: Renamed jacobian to JacobianMat
//
// Revision 1.10  2019/09/16 14:36:30  abbott
// Summary: NewPolyRing no calls LogObsolescentFn; added iroot & IsExactIroot
//
// Revision 1.9  2018/10/02 09:45:30  abbott
// Summary: Moved pseudo-ctors NewPolyRing(CoeffRing, NumIndets, ...) to obsolescent
//
// Revision 1.8  2017/11/20 20:38:27  bigatti
// -- added minimalized
//
// Revision 1.7  2017/11/10 16:02:27  abbott
// Summary: Removed NewLexOrdering, NewStdDegLexOrdering, NewStdDegRevLexOrdering; consequential changes
//
// Revision 1.6  2017/03/29 15:40:40  abbott
// Summary: Now prints UsefulAdvice in log message
//
// Revision 1.5  2017/01/25 13:01:44  abbott
// Summary: Warning message now output to CoCoA::LogStream (instead of clog)
//
// Revision 1.4  2016/11/07 14:16:50  bigatti
// -- added AreGensSquareFreeMonomial
//
// Revision 1.3  2016/11/05 16:34:17  abbott
// Summary: Put LogObsolescentFn into anon namespace
//
// Revision 1.2  2016/11/04 20:44:07  abbott
// Summary: Cleaned and simplified
//
// Revision 1.1  2016/11/03 12:29:58  abbott
// Summary: Added file for obsolescent fns; also there is a global flag saying whether to give error if calling one.
//
//
