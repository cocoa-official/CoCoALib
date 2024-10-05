//   Copyright (c)  2011  Anna Bigatti

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


// Source code for toric ideals

#include "CoCoA/MatrixView.H"
#include "CoCoA/PPMonoid.H"
#include "CoCoA/SparsePolyOps-RingElem.H"
//#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/SparsePolyOps-hilbert.H"
#include "CoCoA/TmpToric.H"
#include "CoCoA/BigInt.H"
#include "CoCoA/ideal.H"
#include "CoCoA/matrix.H"
#include "CoCoA/ring.H"
#include "CoCoA/utils.H"
#include "CoCoA/verbose.H"  // for VerboseLog
#include "TmpHilbertDir/AnnaUtils.h"
#include "TmpHilbertDir/IVectors.h"
#include "TmpHilbertDir/eterms.h"
#include "TmpHilbertDir/unipoly.h" // for PoincareMaxPower
#include "TmpHilbertDir/poincare.h"
#include "TmpHilbertDir/toric.h"

#include <iostream> // for debugging only
#include <vector>
using std::vector;

namespace CoCoA
{

  //----------------------------------------------------------------------
  namespace // anonymous
  {
    ints IndicesForToric(const std::vector<long> v);
    ideal NewIdeal(const SparsePolyRing& P, const BitermList& BL);
    ideal NewIdeal(const SparsePolyRing& P, const BitermList& BL, const ints Weights);
    biterm PPs2Binom(ConstRefPPMonoidElem t1, ConstRefPPMonoidElem t2);
    BitermList NewBitermList(const ideal& I);
    void MatKerToBListAndIndices(ConstMatrixView TrM, BitermList *BL, ints *Indices, ints *Weights);
  }
  //----------------------------------------------------------------------


  ideal SequentialToric_C(const ideal& I, const std::vector<long>& indices)
  {
    // check input
    const SparsePolyRing P = RingOf(I);
    //    ::StartPoincare(0);
    ::StartToric(NumIndets(P));
    VerboseLog VERBOSE("SequentialToric_C(I)");
    VERBOSE(99) << "I := " << I << std::endl;
    BitermList BLIn = NewBitermList(I);
    VERBOSE(99) << "conv " << NewIdeal(P, BLIn) << std::endl;
    ints IndicesIn = IndicesForToric(indices);
    VERBOSE(90) << IntsGetLen(IndicesIn) << " indices:" << std::endl;
    for (long i=1; i<=IntsGetLen(IndicesIn); ++i)
      VERBOSE(90) << IndicesIn[i] << ", " << std::endl;
    VERBOSE(90) << "------------" << std::endl;
    BitermList BLOut = SequentialToric(BLIn, IndicesIn);
    // BLIn and IndicesIn are freed by SequentialToric
    ideal T = NewIdeal(P, BLOut);
    EraseAndFreeBList(BLOut);
    return T;
  }


  ideal SequentialToric_C(const SparsePolyRing& P, ConstMatrixView M)
  {
    // check input
    //    ::StartPoincare(0);
    ::StartToric(NumIndets(P));
    VerboseLog VERBOSE("SequentialToric_C(M)");
    BitermList BLIn;
    ints IndicesIn;
    ints WeightsIn;
    MatKerToBListAndIndices(transpose(M), &BLIn, &IndicesIn, &WeightsIn);    
    VERBOSE(99) << "I := " << NewIdeal(P, BLIn, WeightsIn) << std::endl;
    VERBOSE(90) << IntsGetLen(IndicesIn) << " indices:" << std::endl;
    for (long i=1; i<=IntsGetLen(IndicesIn); ++i)
      VERBOSE(90) << IndicesIn[i] << ", " << std::endl;
    VERBOSE(90) << "------------" << std::endl;
    BitermList BLOut = SequentialToric(BLIn, IndicesIn);
    // BLIn and IndicesIn are freed by SequentialToric
    ideal T = NewIdeal(P, BLOut, WeightsIn);
    EraseAndFreeBList(BLOut);
    ints_free(WeightsIn);
    return T;
  }


  void EndToric_C()
  {
    EndPoincare_C(); // clear *C* global variable
  }
  


  namespace // anonymous
  {

    ints IndicesForToric(const std::vector<long> v)
    {
      long l = len(v);
      ints res=ints_malloc(l);
      IntsSetSize(res,l);
      IntsSetLen(res,0);
      for (long i=0; i<l; ++i)
        IntsPutLast(res, v[i]+1);
      return res;
    }


    RingElem BitermToRingElem(const SparsePolyRing& P, biterm b, ints weights)
    {
      ivec v = b->Vect;
      long n = ivec_len(v);
      if (n>NumIndets(P))
        CoCoA_THROW_ERROR("wrong number of indeterminates", "BitermToRingElem");
      PPMonoid M = PPM(P);
//   if ( ToricAlg == EATI )

      PPMonoidElem pp1(M);
      PPMonoidElem pp2(M);
      for ( long i=1 ; i<=n ; ++i )
        if (ivec_nth(v,i) > 0)
          M->myMulIndetPower(raw(pp1), i-1, ivec_nth(v,i)/weights[i]);
        else if (ivec_nth(v,i) < 0)
          M->myMulIndetPower(raw(pp2), i-1, -ivec_nth(v,i)/weights[i]);
      return monomial(P, pp1) - monomial(P, pp2);
    }


    ideal NewIdeal(const SparsePolyRing& P, const BitermList& BL)
    {
      biterms Bs = Biterms(BL);
      const long n = BListLen(BL);
      std::vector<RingElem> res;

      ints weights1;
      weights1 = ints_malloc(NumIndets(P));
      IntsSetSize(weights1, NumIndets(P));
      IntsSetLen(weights1, NumIndets(P));
      for (long i=1 ; i<=NumIndets(P) ; i++ )  weights1[i] = 1;
      for (long i=1 ; i<=n ; i++ )
        //        if ( ToricAlg != EATI || ivec_nth(Bs[i]->Vect, IndNo)==0 )
        res.push_back(BitermToRingElem(P, Bs[i], weights1));
      ints_free(weights1);
      return ideal(P, res);
    }
    
    
    ideal NewIdeal(const SparsePolyRing& P, const BitermList& BL, const ints weights)
    {
      biterms Bs = Biterms(BL);
      const long n = BListLen(BL);
      std::vector<RingElem> res;
      for (long i=1 ; i<=n ; i++ )
        //        if ( ToricAlg != EATI || ivec_nth(Bs[i]->Vect, IndNo)==0 )
        res.push_back(BitermToRingElem(P, Bs[i], weights));
      return ideal(P, res);
    }
    
    
    biterm PPs2Binom(ConstRefPPMonoidElem t1, ConstRefPPMonoidElem t2)
    {
      ivec Vect;
      ivec_elem n=NumIndets(owner(t1));
      int ord;
 
//   if ( ToricAlg == EATI )
//   {  
//     Vect = ivec_init(i+1);
//     ivec_set_nth(Vect, i+1, 0);
//     ord = i+1;
//   }
//   else
      {
        Vect = ivec_init(n);
        ord = 1;
      }
      const vector<long> exps1 = exponents(t1);
      const vector<long> exps2 = exponents(t2);
      const long MaxExp = std::numeric_limits<int>::max()/n;
      for ( long i=0 ; i<n ; ++i )
      {
        if (exps1[i] >= MaxExp || exps2[i] >= MaxExp)
          CoCoA_THROW_ERROR("exponents too big in binomial", "PPs2Binom");
        ivec_set_nth(Vect, (unsigned long)i+1, exps1[i]-exps2[i]);
//       if ( term_degree(t1)!=term_degree(t2) )
//       {
//         ivec_free(Vect);
//         return nullptr;
//       }
      }
      return BitermMake(Vect, ord);
    }


    BitermList NewBitermList(const ideal& I)
    {
      const SparsePolyRing P = RingOf(I);
      const std::vector<RingElem>& GensI = gens(I);
      BitermList BL = BLNew(len(GensI), NumIndets(P));;
      biterms Bs = Biterms(BL);
      long MaxDeg = 0, BLLen = 0;

      RingElem g(P);
      for ( long i=0 ; i<len(GensI) ; ++i )
      {
        g = GensI[i];
        Bs[++BLLen] = PPs2Binom(LPP(g), LPP(g-monomial(P, LC(g), LPP(g))));
//       if ( Bs[BLLen]==nullptr )
//       {
// 	EraseAndFreeBList(BL);
// 	math_error = cocoa_error("BAD_PARAMS");
// 	return nullptr;
//       }
        BListSetLen(BL, BLLen);
        if ( BitermLtDeg(Bs[BLLen]) > MaxDeg)  MaxDeg = BitermLtDeg(Bs[BLLen]);
      }
      BListSetMaxDeg (BL, MaxDeg);
      
      return BL;
    }
    
    
    void MatKerToBListAndIndices(ConstMatrixView M, BitermList *BL, ints *Indices, ints *Weights)
    {
      VerboseLog VERBOSE("MatKerToBListAndIndices");
      int dim;
      int **small_basis, matrices =8, CurrLen =0, IndNo;
      ints Indices1, weights;  
      biterm B /*, MaxDegB*/; // MaxDegB set but never used!!!
      biterms Bs;
      ivec aux_V; 
      int Index, MaxDeg=0;

      const long nrows = NumRows(M);
      const long ncols = NumCols(M);
      
//   if ( ToricAlg == EATI )
//     IndNo = nrows+1;
//   else
    IndNo = nrows;

  /*  TODO CONTROLLA!!!!!!!!!!!!!!!!!!!!!!!!!!!
  rum_init(sizeof(expType)*(IndNo+1), RUM_STD_SIZE);
  rum_init(sizeof(expType)*(IndNo+2), RUM_STD_SIZE);
  */

  /* INDICES */
  Indices1 = ints_malloc(ncols);
  IntsSetSize(Indices1, ncols);
  IntsSetLen(Indices1, 0);
  *Indices = Indices1;
  INDICES = Indices1;
  
  /* M --> small_matrix */
  int **Msmall;
  Msmall = (int**)malloc(nrows*sizeof(int*));
  int tmp_int;
  BigInt tmp_BigInt;
  for (long i=0; i<nrows; ++i)
  {
    Msmall[i] = (int*)malloc(ncols*sizeof(int));
    for (long j=0; j<ncols; ++j)
    {
      if (!IsInteger(tmp_BigInt, M(i,j)))
        CoCoA_THROW_ERROR(ERR::BadArg, "Toric/MatKerToBListAndIndices");
      if (!IsConvertible(tmp_int, tmp_BigInt))
        CoCoA_THROW_ERROR(ERR::BadArg, "Toric/MatKerToBListAndIndices");
      Msmall[i][j] = tmp_int;
    }
  }
  /* WEIGHTS */
  weights = ints_init(nrows);
  IntsSetLen(weights, nrows);
  int i, j;
  for (i=0; i<nrows; ++i)
  {
    weights[i+1]=0;
    for (j=0; j<ncols; ++j)  weights[i+1] += Msmall[i][j];
  }
  i = 1;
  while ( weights[i+1]==weights[i] )  if ((++i)==nrows) break;
  if ((i)==nrows)
    for (j=1; j <= nrows; j++)  weights[j] =1;
  else
  {
    i = 0;
    while ( Msmall[i][0]!=0 )  if ((++i)==nrows) break;    
    if ((i)==nrows)
      for (j=0; j<nrows; j++)  weights[j+1] = Msmall[j][0];
  }
  *Weights = weights;  
  /* NULL_SPACE TO BLIST */
  dim = null_space( &small_basis, Msmall, nrows, ncols);  
//   for (i=0 ; i<dim ; i++) // print null_space
//   {
//     for (Index=0 ; Index<IndNo ; Index++ )
//       std::cout << small_basis[i][Index] << "\t ";
//     std::cout << std::endl;    
//   }
  
  *BL = BLNew (matrices*dim, IndNo);
  BListSetLen (*BL, matrices*dim);
  Bs = Biterms(*BL);
  int MaxExp = std::numeric_limits<int>::max()/(IndNo); // guarantee conv binom
  MaxExp /= IndNo; // *extra margin* for GBasis/toric computation
  for ( i=0 ; i<dim ; i++)
  {
    aux_V =ivec_init(IndNo);
    for ( Index=0 ; Index<IndNo ; Index++ )
      if ( Index>=nrows )
        ivec_set_nth(aux_V, Index+1, 0);
      else
      {
        if (std::abs(small_basis[i][Index]) >= MaxExp/weights[Index+1])
          CoCoA_THROW_ERROR("exponents too big in binomial", "MatKerToBListAndIndices");
        ivec_set_nth(aux_V, Index+1, small_basis[i][Index]*weights[Index+1]);
      }
    B = BitermMake(aux_V,1); // frees aux_V
    //VERBOSE(99) << "+++++  BitermLtDeg(B) = " << BitermLtDeg(B) << std::endl;
    Bs[++CurrLen] = B;
    if ( BitermLtDeg(B) > MaxDeg) { MaxDeg = BitermLtDeg(B); /*MaxDegB = B*/; }
  }
  for (i=0; i<dim; i++) free(small_basis[i]);
  free(small_basis);
  // "matrices"-1 more new null_spaces of M
  for ( --matrices ; matrices>0 ; --matrices )
  {
    dim = null_space( &small_basis, Msmall, nrows, ncols);  
//     for (i=0 ; i<dim ; i++)
//     {
//       for (Index=0 ; Index<IndNo ; Index++ )
//         std::cout << small_basis[i][Index] << "\t ";
//       std::cout << std::endl;    
//     }
    for ( i=0 ; i<dim ; i++ )
    {
      aux_V =ivec_init(IndNo);
      for ( Index=0 ; Index<IndNo ; Index++ )
        if ( Index>=nrows )
          ivec_set_nth(aux_V, Index+1, 0);
        else
        {
          if (std::abs(small_basis[i][Index]) >= MaxExp/weights[Index+1])
            CoCoA_THROW_ERROR("exponents too big in binomial", "MatKerToBListAndIndices");
          ivec_set_nth(aux_V, Index+1, small_basis[i][Index]*weights[Index+1]);
        }
      B = BitermMake(aux_V,1); // frees aux_V
      Bs[++CurrLen] = B; 
      if ( BitermLtDeg(B) > MaxDeg) { MaxDeg = BitermLtDeg(B); /*MaxDegB = B;*/ }
    }
    for (i=0; i<dim; i++) free(small_basis[i]);
    free(small_basis);
  }
  for (i=0; i<nrows ; i++) free(Msmall[i]);
  free(Msmall);
  BListSetMaxDeg (*BL, MaxDeg);

  /*
  if (trouble) return -dim;
  return dim;
  */
}

  } // anonymous namespace



} // end of namespace CoCoA

// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/TmpToric.C,v 1.27 2023/03/09 22:41:07 abbott Exp $
// $Log: TmpToric.C,v $
// Revision 1.27  2023/03/09 22:41:07  abbott
// Summary: Tidied use of exponents (redmine 1694)
//
// Revision 1.26  2022/08/05 13:14:22  abbott
// Summary: Replaced NumRows/NumCols in for-lop exit cond by local vars
//
// Revision 1.25  2022/03/07 10:58:54  bigatti
// Summary: commented out MaxDegB (unused)
//
// Revision 1.24  2022/02/18 14:12:01  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.23  2020/06/17 15:49:29  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.22  2019/03/19 11:07:08  abbott
// Summary: Replaced 0 by nullptr where appropriate
//
// Revision 1.21  2019/03/04 11:13:20  abbott
// Summary: Removed redundant include directive
//
// Revision 1.20  2018/08/06 12:50:43  bigatti
// -- stricter limit for exponents in MatKerToBListAndIndices
//
// Revision 1.19  2018/06/25 09:28:08  bigatti
// -- set limit to exponents in input: MaxDeg/numindets
//
// Revision 1.18  2018/05/18 16:38:51  bigatti
// -- added include SparsePolyOps-RingElem.H
//
// Revision 1.17  2018/04/10 14:30:49  bigatti
// -- fixed includes
//
// Revision 1.16  2017/03/29 14:43:56  bigatti
// -- modified verbosity
//
// Revision 1.15  2017/02/14 17:34:34  bigatti
// -- added verbosity
//
// Revision 1.14  2017/02/08 17:03:09  abbott
// Summary: Changed indentation to avoid compiler warnings
//
// Revision 1.13  2015/06/11 16:57:24  bigatti
// -- using new functions monomial(ring, pp) and monomial(ring, expv)
//
// Revision 1.12  2014/07/30 14:11:54  abbott
// Summary: Changed name AmbientRing --> RingOf
// Author: JAA
//
// Revision 1.11  2014/07/07 13:16:30  abbott
// Summary: Removed AsSparsePolyRing
// Author: JAA
//
// Revision 1.10  2014/04/17 13:39:39  bigatti
// -- MatrixViews --> MatrixView
//
// Revision 1.9  2013/02/04 14:25:35  bigatti
// -- free weights1 (memory leak found with valgrind)
//
// Revision 1.8  2013/02/01 17:57:29  bigatti
// -- added EraseAndFreeBList in SequentialToric_C
//
// Revision 1.7  2013/01/31 17:11:37  bigatti
// -- renamed a variable
//
// Revision 1.6  2013/01/31 13:14:02  bigatti
// -- removed comments
//
// Revision 1.5  2013/01/31 13:13:32  bigatti
// -- added free
//
// Revision 1.4  2013/01/31 12:57:01  bigatti
// -- investigating memory leak - clarifying code
//
// Revision 1.3  2013/01/25 13:49:06  bigatti
// -- nrows->ncols in  Msmall[i] = (int*)malloc(ncols*sizeof(int));  (segv)
//
// Revision 1.2  2011/08/14 15:52:16  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.1  2011/05/09 14:48:11  bigatti
// -- first import
//
