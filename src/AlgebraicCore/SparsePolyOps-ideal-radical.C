//   Copyright (c)  2023  John Abbott,  Anna M. Bigatti
//   Original authors: Nicolas Jagersma, Alice Moallemy

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


#include "CoCoA/SparsePolyOps-ideal.H"

#include "CoCoA/DenseMatrix.H"
#include "CoCoA/PPMonoid.H" // for PPMonoid
#include "CoCoA/RingZZ.H" // for ideal
#include "CoCoA/SparsePolyIter.H" // for SparsePolyIter
#include "CoCoA/SparsePolyOps-RingElem.H"
#include "CoCoA/SparsePolyOps-hilbert.H"
#include "CoCoA/TmpGReductor.H" // for IsHomog(vector)
#include "CoCoA/TmpPPVector.H" // for ideal
#include "CoCoA/VectorOps.H"  // for len
#include "CoCoA/factor.H"
#include "CoCoA/ideal.H" // for ideal
#include "CoCoA/matrix.H" // for ideal
#include "CoCoA/time.H"
#include "CoCoA/verbose.H" // for VerboseLog
// #include "CoCoA/interrupt.H"

using std::cout;
using std::endl;
using std::min;
using std::vector;


namespace CoCoA
{

//   const string SYNTAX =
//   "SYNTAX \n"
//   "    radical(I: IDEAL): IDEAL \n"
//   "    EquiIsoDec(I: IDEAL): List of IDEALs \n"
//   "    RadicalOfUnmixed(I: IDEAL): IDEAL \n";

// const string DESCRIPTION =
//   "DESCRIPTION \n"
//   "Let P be a polynomial ring and I an ideal. This package computes the \n"
//   "radical of I using the algorithm described in the paper \n"
//   " \n"
//   "M. Caboara, P.Conti and C. Traverso:  Yet Another Ideal Decomposition \n"
//   "Algorithm. Proc. AAECC-12, pp 39-54, 1997, Lecture Notes in Computer \n"
//   "Science, n.1255 Springer-Verlag. \n"
//   "This package has been used as prototype in the development of the radical \n"
//   "implementation for the FRISCO project. \n"
//   " \n"
//   "Coefficient fields : at the moment this implementation works only for the \n"
//   "rationals or if the characteristic is big enough. \n"
//   " \n"
//   "    There are three main functions : \n"
//   " \n"
//   "            EquiIsoDec(I); \n"
//   "    computes an equimensional isoradical decomposition of I, i.e. \n"
//   "    a list of unmixed ideals I_1,..,I_k such that the radical of I is the \n"
//   "    intersection of the radicals of I_1,..,I_k. Redundancies are possible. \n"
//   " \n"
//   "            RadicalOfUnmixed(I); \n"
//   "    computes the radical of an unmixed ideal. \n"
//   " \n"
//   "            radical(I); \n"
//   "    computes the radical of an ideal I. \n"
//   " \n"
//   "In later releases of the package some of the functions will be ported \n"
//   "to C and the algorithms for intersection, transporter and saturation \n"
//   "will be optimized. \n"
//   " \n"
//   ">EXAMPLES< \n"
//   " \n"
//   "    Use R ::= QQ[x,y,z]; \n"
//   "    I:=intersection(ideal(x-1,y-1,z-1),ideal(x-2,y-2)^2,ideal(x)^3); \n"
//   " \n"
//   "    H:=EquiIsoDec(I); \n"
//   "    H; \n"
//   " \n"
//   "    [ideal(x), \n"
//   "     ideal(z - 1, y - 1, x - 1), \n"
//   "     ideal(xy - y^2 - 2x + 2y, x^2 - y^2 - 4x + 4y, \n"
//   "           y^2z - y^2 - 4yz + 4y + 4z - 4, y^3 - 5y^2 + 8y - 4, \n"
//   "	   x - 2)] \n"
//   " \n"
//   "     [RadicalOfUnmixed(I)|I In H]; \n"
//   " \n"
//   " \n"
//   "     [ideal(x), \n"
//   "      ideal(z - 1, y - 1, x - 1), \n"
//   "      ideal(xy - y^2 - 2x + 2y, x^2 - y^2 - 4x + 4y, \n"
//   "            y^2z - y^2 - 4yz + 4y + 4z - 4, y^3 - 5y^2 + 8y - 4, \n"
//   "	    x - 2, y - 2)] \n"
//   " \n"
//   "     Note : the ideals are not necessarily presented in the simplest form. \n"
//   "            The third ideal is in fact ideal(y - 2, x - 2); \n"
//   " \n"
//   " \n"
//   "     radical(I); \n"
//   " \n"
//   "     ideal(xyz - xy - 2xz + 2x, xy^2 - 3xy + 2x, x^2 - xy) \n";

  namespace // anonymous
  {
    // from myInterreduce.C

  /*   bool IsSubset(const vector<RingElem>& L1, const vector<RingElem>& L2){ */
  /*   long i = 0; */
  /*   long j = 0; */
  /*   //int c = 0; */
  /*   for (i = 0; i < len(L1); ++i){ */
  /*     for (j = 0; j < len(L2); ++j){ */
  /*       if (L1[i] == L2[j]) break; */
  /*     } */
  /*     if (j == len(L2)) return false; */
  /*   } */
  /*   return true; */
  /* } // IsSubset(L1, L2) */


  /* // important for diff */
  /* bool IsSubset(const vector<vector<PPMonoidElem> >& LP1, const vector<vector<PPMonoidElem> >& LP2){ */
  /*   long i = 0; */
  /*   long j = 0; */
  /*   for (i = 0; i < len(LP1); ++i){ */
  /*     for (j = 0; j < len(LP2); ++j){ */
  /*       if (LP1[i] == LP2[j]) break; */
  /*     } */
  /*     if (j == len(LP2)) return false; */
  /*   } */
  /*   return true; */
  /* } // IsSubset(LP1, LP2) */


  /* bool IsSubset(const vector<PPMonoidElem>& P1, const vector<PPMonoidElem>& P2){ */
  /*   long i = 0; */
  /*   long j = 0; */
  /*   for (i = 0; i < len(P1); ++i){ */
  /*     for (j = 0; j < len(P2); ++j){ */
  /*       if (P1[i] == P2[j]) break; */
  /*     } */
  /*     if (j == len(P2)) return false; */
  /*   } */
  /*   return true; */
  /* } // IsSubset(P1, P2) */


  // prova queste e cancella IsSubset
    bool IsIn(const RingElem& f, const vector<RingElem>& L)
    {
      for (const auto& g:L)  if (f == g) return true;
      return false;
    }


    bool IsIn(const PPMonoidElem& f, const vector<PPMonoidElem>& L)
    {
      for (const auto& g:L)  if (f == g) return true;
      return false;
    }
    

    // important for MonomialCleaning
    bool EqSet(const vector<RingElem>& L1, const vector<RingElem>& L2)
    {
      for (const auto& f:L1)  if (!IsIn(f,L2)) return false;
      for (const auto& f:L2)  if (!IsIn(f,L1)) return false;
      return true;
    } // EqSet(L1, L2)
    
    
    // important for MakeVarsList 
    bool EqSet(const vector<PPMonoidElem>& L1, const vector<PPMonoidElem>& L2)
    {
      for (const auto& f:L1)  if (!IsIn(f,L2)) return false;
      for (const auto& f:L2)  if (!IsIn(f,L1)) return false;
      return true;
    } // EqSet(P1, P2)
    

    bool IsIn(const vector<PPMonoidElem>& v, const vector<vector<PPMonoidElem>>& L)
    {
      for (const auto& V:L)  if (EqSet(v, V)) return true;
      return false;
    }



  // diff(R1, R2) -> R1\R2
  vector<RingElem>  DIFF(const vector<RingElem>& R1, const vector<RingElem>& R2)
  {
    vector<RingElem>  K;
    for (const auto& f:R1)  if (!IsIn(f, R2)) K.push_back(f);
    return K;
  } // diff(R1, R2)


  // diff(T1, T2) -> T1\T2
  vector<PPMonoidElem>  DIFF(const vector<PPMonoidElem>& T1, const vector<PPMonoidElem>& T2){
    vector<PPMonoidElem>  K;
    for (const auto& t:T1)  if (!IsIn(t, T2)) K.push_back(t);
    return K;
  } // diff(T1, T2)


  vector<long> exponents(const RingElem& f){
    PPMonoidElem Lppf = LPP(f);
    vector<long> expoVec;
    exponents(expoVec, Lppf);
    return expoVec;
  } // exponents(f)


  // important for T2TT and T2FT
  RingElem MakeTerm(const ring& P, const vector<long>& expo)
  {
    if (len(expo) != NumIndets(P))
    {
      cout << "MakeTerm: exponent list has the wrong length" << endl;
      return zero(P); //just to keep the compiler quiet
    }
    RingElem result = one(P);
    for (long i = 0; i < len(expo); ++i)
      result *= power(indet(P, i), expo[i]);
    return result;
  } // MakeTerm(P, expo)


  RingElem LT(const RingElem& f)
  {
    return MakeTerm(owner(f), exponents(f));
  } // LT(f)


  RingElem First(const vector<RingElem>& L)
  {
    if (L.empty())
    {
      RingElem emptyL;
      return emptyL;
    }
    return L.front();
  } // First(L)


  PPMonoidElem First(const vector<PPMonoidElem>& P)
  {
    return P.front();
  } // First(P)


  vector<RingElem> Tail(const vector<RingElem>& L)
  {
    vector<RingElem> result;
    for (long i=1; i < len(L); ++i)  result.push_back(L[i]);
    return result;
  } // Tail(L)


  vector<PPMonoidElem> Tail(const vector<PPMonoidElem>& P)
  {
    vector<PPMonoidElem> result;
    for (long i=1; i < len(P); ++i)  result.push_back(P[i]);
    return result;
  } // Tail(P)



  } // end of namespace anonymous
  


    //  namespace { // anonymous
    ideal radical_0dimDRL(const ideal& I)
    {
//      const SparsePolyRingBase::IdealImpl* const ptrI = 
      const auto ptrI = SparsePolyRingBase::IdealImpl::ourGetPtr(I);
      return ptrI->myRadical_0dimDRL(); // behaves differently from other memb fns
    }
  
  
  ideal radical_0dim(const ideal& I)
  {
    const PolyRing& P = RingOf(I);
    if (!IsField(CoeffRing(P)))  CoCoA_THROW_ERROR1(ERR::ReqField);
    //    if (HasStdDegRevLex(P))  return radical_0dimDRL(I);
    if (HasGBasis(I)) return radical_0dimDRL(I);
    PolyRing P_drl = NewPolyRing(CoeffRing(P), NewSymbols(NumIndets(P)));
    const RingHom phi = PolyAlgebraHom(P, P_drl, indets(P_drl));
    const RingHom psi = PolyAlgebraHom(P_drl, P, indets(P));
    const ideal RadI = radical_0dimDRL(ideal(phi(gens(I))));
    return ideal(psi(gens(RadI)));
  } // radical_0dim(I)


  ideal radical_MonId(const ideal& I)
  {
    CoCoA_ASSERT(IsField(CoeffRing(RingOf(I))));
    VerboseLog VERBOSE("radical_MonId");
    VERBOSE(1000) << " starting " << std::endl;
//    const SparsePolyRingBase::IdealImpl* const ptrI =
    const auto ptrI = SparsePolyRingBase::IdealImpl::ourGetPtr(I);
    return ptrI->myRadical_MonId(); // behaves differently from other memb fns
  }


  //  } // anonymous namespace


////////////  radical general ---------------------------------------
  // global variables:
  static long RADICAL_NEWRINGNUM = 0;
  static bool RADICAL_FULLSqFreeSPLIT = false;
  static bool RADICAL_FULLGISPLIT = false;
  static bool RADICAL_FULLGCDSPLIT = false;
  static bool RADICAL_SPLITTINGGCD = false;
  static bool RADICAL_BISATSPLITTING = false;
  static long RADICAL_SAT = 1;
  static long RADICAL_BRANCHING = 0;

  void PrintSettings(){
    cout << "RADICAL_NEWRINGNUM " << RADICAL_NEWRINGNUM << endl;
    cout << "RADICAL_FULLSqFreeSPLIT " << RADICAL_FULLSqFreeSPLIT << endl;
    cout << "RADICAL_FULLGISPLIT " << RADICAL_FULLGISPLIT << endl;
    cout << "RADICAL_FULLGCDSPLIT " << RADICAL_FULLGCDSPLIT << endl;
    cout << "RADICAL_SPLITTINGGCD " << RADICAL_SPLITTINGGCD << endl;
    cout << "RADICAL_BISATSPLITTING " << RADICAL_BISATSPLITTING << endl;
    cout << "RADICAL_SAT " << RADICAL_SAT << endl;
    cout << "RADICAL_BRANCHING " << RADICAL_BRANCHING << endl;
  }
  
  // important for EqSet
  /* bool IsSubset(const vector<RingElem>& L1, const vector<RingElem>& L2){
     long i = 0;
     long j = 0;
     for (i = 0; i < len(L1); ++i){ 
     for (j = 0; j < len(L2); ++j){
     if (L1[i] == L2[j]) break;
     }
     if (j == len(L2)) return false;
     }
     return true;
     } // IsSubset(L1, L2)
  */

  // important for EqSet
  /*bool IsSubset(const vector<PPMonoidElem>& P1, const vector<PPMonoidElem>& P2){
    long i = 0;
    long j = 0;
    for (i = 0; i < len(P1); ++i){
    for (j = 0; j < len(P2); ++j){
    if (P1[i] == P2[j]) break;
    }
    if (j == len(P2)) return false;
    }
    return true;
    } // IsSubset(P1, P2)

  */
  // important for diff
  /* bool IsSubset(const vector<vector<PPMonoidElem> >& LP1, const vector<vector<PPMonoidElem> >& LP2){
     long i = 0;
     long j = 0;
     for (i = 0; i < len(LP1); ++i){
     for (j = 0; j < len(LP2); ++j){
     if (LP1[i] == LP2[j]) break;
     }
     if (j == len(LP2)) return false;
     }
     return true;
     } // IsSubset(LP1, LP2)
  */



  // important for T2TT and T2FT
  /*vector<long> exponents(const RingElem& f){
    PPMonoidElem Lppf = LPP(f);
    vector<long> expoVec;
    exponents(expoVec, Lppf);
    return expoVec;
    } // exponents(f)


    // important for T2TT and T2FT
    RingElem MakeTerm(const ring& P, vector<long> expo){
    if (len(expo) != NumIndets(P)){
    // CoCoA_THROW_ERROR2(ERR::BadArg, "MakeTerm: exponent list has the wrong length");
    cout << "MakeTerm: exponent list has the wrong length" << endl;
    return zero(P); //just to keep the compiler quiet
    }
    RingElem result = one(P);
    for (long i = 0; i < len(expo); ++i){
    result *= power(indet(P, i), expo[i]);
    }
    return result;
    } // MakeTerm(P, expo)


    RingElem LT(RingElem f){
    return MakeTerm(owner(f), exponents(f));
    } // LT(f)
  


    // Functions for help
    RingElem First(vector<RingElem> L){
    if (L.empty()){
    RingElem emptyL;
    return emptyL;
    }
    return L.front();
    } // First(L)
  */


  /* PPMonoidElem First(vector<PPMonoidElem> P){
     return P.front();
     } // First(P)
  */


  ideal First(const vector<ideal>& I){
    return I.front();
  } // First(I)


  vector<PPMonoidElem> First(const vector<vector<PPMonoidElem> >& PP){
    if (PP.empty()){
      vector<PPMonoidElem> emptyPP;
      return emptyPP;
    }
    return PP.front();
  } // First(PP)


  /*vector<RingElem> Tail(vector<RingElem> L){
    vector<RingElem> result;
    for (long i=1; i < len(L); ++i){
    result.push_back(L[i]);
    }
    return result;
    } // Tail(L)
  */


  /*  vector<PPMonoidElem> Tail(vector<PPMonoidElem> P){
      vector<PPMonoidElem> result;
      for (long i=1; i < len(P); ++i){
      result.push_back(P[i]);
      }
      return result;
      } // Tail(P)
  */


  vector<ideal> Tail(const vector<ideal>& I){
    vector<ideal> result;
    for (long i=1; i < len(I); ++i){
      result.push_back(I[i]);
    }
    return result;
  } // Tail(L)


  vector<vector<PPMonoidElem> > Tail(const vector<vector<PPMonoidElem> >& PP){
    vector<vector<PPMonoidElem> > result;
    for (long i=1; i < len(PP); ++i){
      result.push_back(PP[i]);
    }
    return result;
  } // Tail(PP)


  //    questa DIFF serve? -----------------------
  // diff(L1, L2) -> L1\L2
  /* vector<vector<PPMonoidElem> > DIFF(vector<vector<PPMonoidElem> > L1, vector<vector<PPMonoidElem> > L2){ */
  /*   vector<vector<PPMonoidElem> > K; */
  /*   for (long i=0; i < len(L1); ++i){ */
  /*     vector<vector<PPMonoidElem> > tmpElem; */
  /*     tmpElem.push_back(L1[i]); */
  /*     if (!IsSubset(tmpElem, L2)) K.push_back(L1[i]); */
  /*     tmpElem.clear(); */
  /*   } */
  /*   return K; */
  /* } // diff(L1, L2) */


  // diff(T1, T2) -> T1\T2
  /* vector<PPMonoidElem>  diff(vector<PPMonoidElem> T1, vector<PPMonoidElem> T2){
     vector<PPMonoidElem>  K;
     for (long i=0; i < len(T1); ++i){
     vector<PPMonoidElem> tmpElem;
     tmpElem.push_back(T1[i]);
     if (!IsSubset(tmpElem, T2)) K.push_back(T1[i]);
     tmpElem.clear();
     }
     return K;
     } // diff(T1, T2)
  */


  // diff(R1, R2) -> R1\R2
  /*vector<RingElem>  diff(vector<RingElem> R1, vector<RingElem> R2){
    vector<RingElem>  K;
    for (long i=0; i < len(R1); ++i){
    vector<RingElem> tmpElem;
    tmpElem.push_back(R1[i]);
    if (!IsSubset(tmpElem, R2)) K.push_back(R1[i]);
    tmpElem.clear();
    }
    return K;
    } // diff(R1, R2)
  */


  ideal MonomialCleaning(const ideal& I){
    VerboseLog VERBOSE("MonomialCleaning: ");
    bool IsFirstLoop = true;
    ideal result = I;
    vector<RingElem> G = gens(result);
    do{
      G = gens(result);
      bool MonFound = false;
      for (long j = 0; j < len(G); ++j){
        //if (IsMonomial(G[j]) && !IsSqFree(LPP(G[j]))){
        if (IsMonomial(G[j])){
          MonFound = true;
          G[j] = radical(LT(G[j]));
        }
      }
      if (MonFound || IsFirstLoop){
        vector<RingElem> W = interreduced(G);
        //	cout << "In MonomialCleaning G = " << G << endl;
        result = ideal(gens(ideal(W)));
      }
      IsFirstLoop = false;
      
    }while(!EqSet(G, gens(result)));
    return result;
  } // MonomialCleaning(I)


  // important for WhichIndets
  PPMonoidElem ProductOfSupport(const RingElem& f){
    const PPMonoid PMonoid = PPM(RingOf(ideal(f)));
    PPMonoidElem pElem = one(PMonoid);
    for (SparsePolyIter iter=BeginIter(f); !IsEnded(iter); ++iter){
      pElem *= PP(iter);
    }
    return pElem;
  } // ProductOfSupport(f)


  // in package called WhichVars
  vector<PPMonoidElem> WhichIndets(const PPMonoidElem& p){
    vector<PPMonoidElem> result;
    vector<PPMonoidElem> x = indets(owner(p));
    for (long i = 0; i < len(x); ++i){
      if (IsDivisible(p, x[i])) result.push_back(x[i]);
    }
    return result;
  } // WhichIndets(p)


  // in package called WhichVars
  vector<PPMonoidElem> WhichIndets(const RingElem& f){
    return WhichIndets(ProductOfSupport(f)); //WhichVars(product(support(f)));
  } // WhichVars(f)


  // in package called WhichVars
  vector<PPMonoidElem> WhichIndets(const ideal& I){
    RingElem f = one(RingOf(I));
    vector<RingElem> gensI = gens(I);
    for (long i = 0; i < len(gensI); ++i){
      gensI[i] = monomial(RingOf(I), ProductOfSupport(gensI[i]));
      f *= gensI[i];
    }
    return WhichIndets(f);
  } // WhichVars(L)


  ideal SqFreeMonId(const ideal& I){
    vector<RingElem> T = gens(I);
    for (int i = 0; i < len(T); ++i){
      T[i] = radical(T[i]);
    }
    return IdealOfGBasis(ideal(T));
  } // SqFreeMonId(I)


  // important for SPLITTING
  // ideal SaturationIJ(const ideal& I, const ideal& J){
  //   return saturate(I,J);
  // } // SaturationIJ(I, J)


  vector<ideal> NormalSplitting(const ideal& I, const RingElem& f){
    vector<ideal> result;
    result.push_back(saturate(I, ideal(f)));
    result.push_back(I+ideal(f));
    //     cout << "result of Normal splitting: " << result << endl;
    return result;
  } // NormalSplitting(I, f)

  
  vector<ideal> Splitting(const ideal& I, const RingElem& f)
  {
    if (RingOf(I) == RingOf(ideal(f))) return NormalSplitting(I,f);
    CoCoA_THROW_ERROR1(ERR::MixedRings);
    vector<ideal> result; //just to keep the compiler quiet
    result.push_back(I);
    result.push_back(I);
    return result;
  } // Splitting(I, f)

  
  vector<PPMonoidElem> MakeSet(vector<PPMonoidElem> P)
  {
    vector<PPMonoidElem> result;
//    while (!P.empty()){
      /* PPMonoidElem firstElem = First(P); */
      /* vector<PPMonoidElem> tmp; */
      /* tmp.push_back(firstElem); */
      /* P = Tail(P); */
      /* if (!IsSubset(tmp, P)) result.push_back(firstElem); */
      /* //firstElem.clear(); */
      /* tmp.clear(); */
      for (const auto& f:P)  if (!IsIn(f, result))  result.push_back(f);
//    }
    //std::sort(result.begin(), result.end());
    return result;
  }


  matrix MakeBlockOrder(const vector<PPMonoidElem>& BiggerVars)
  {
    VerboseLog VERBOSE("MakeBlockOrder: ");
    VERBOSE(40) << "BiggerVars: " << BiggerVars << endl;
    vector<PPMonoidElem> BiggerVars2;
    for (long i=0; i<len(BiggerVars); ++i)  BiggerVars2.push_back(BiggerVars[i]);
    BiggerVars2 = MakeSet(BiggerVars2);  // to eliminate duplicates
    VERBOSE(50) << "BiggerVars2: " << BiggerVars2 << endl;
	   
    const long N = NumIndets(owner(BiggerVars[0]));
    matrix OrdM = NewDenseMat(RingZZ(), 2, N);
    // This loop does not check whether there are repeated indets
    for (long i = 0; i < len(BiggerVars2); ++i)
    {
      long j;
      if (!IsIndet(j, BiggerVars2[i]))  CoCoA_THROW_ERROR1(ERR::ReqIndet);
      SetEntry(OrdM, 0, j, 1); // L[IndetIndex(BiggerVars[i])] = 1;
    }
    for (long i=0; i < N; ++i)  SetEntry(OrdM, 1, i, 1-OrdM(0,i));
    return MakeTermOrdMat(OrdM);
  } // MakeBlockOrder(BiggerVars)


  // important for PutPolyLInSubstring / instead of diff(..); 
  bool IsIntersectionEmpty(const vector<PPMonoidElem>& P1, const vector<PPMonoidElem>& P2)
  {
    if (P1.empty() || P2.empty()) return true; // ?
    if (owner(P1[0]) != owner(P2[0]))
    {
      CoCoA_THROW_ERROR1(ERR::MixedRings);
      return false; //just to keep the compiler quite
    }
    for (long i = 0; i < len(P1); ++i)
      for (long j = 0; j < len(P2); ++j)
        if (P1[i] == P2[j]) return false; 
    return true;
  }

  
  vector<RingElem> PutPolyLInSubring(const vector<RingElem>& L, const vector<PPMonoidElem> Inds)
  {
    vector<RingElem> result;
    for (long i = 0; i < len(L); ++i)
      if (IsIntersectionEmpty(Inds, WhichIndets(L[i]))) result.push_back(L[i]);
    return result;
  }


  vector<ideal> PutPolyInSubring(const vector<ideal>& K, const vector<PPMonoidElem> Inds)
  {
    vector<ideal> result;
    for (long i = 0; i < len(K); i++)
      if (IsIntersectionEmpty(Inds, WhichIndets(K[i])))  result.push_back(K[i]);
    return result;
  } // PutPolyInSubring(L, Inds)


  /*
    vector<PPMonoidElem> AlgIndVarsRec(const vector<PPMonoidElem>& L1, long codim){
    vector<PPMonoidElem> L;
    if (L1.empty() == true) return L;
    L = L1;
    if (codim == 0) return L;
    vector<PPMonoidElem> K = WhichIndets(First(L));
    vector<PPMonoidElem> T = Tail(L);
    vector<PPMonoidElem> Vars;
    PPMonoidElem X = L1[0];
    do{
    X = First(K);
    vector<PPMonoidElem> XList;
    XList.push_back(X);
    K = Tail(K);
    vector<PPMonoidElem> L2;
    for (long i = 0; i < len(T); i++){
    if(!IsDivisible(T[i], X)){
	  L2.push_back(T[i]);
    }
    }
    Vars.push_back(X);
    vector<PPMonoidElem> tmp = AlgIndVarsRec(L2, codim-1);
    for (long j = 0; j < len(tmp); j++){
    Vars.push_back(tmp[j]);
    }


    if(len(Vars) == 0) break;
    }while(K.empty() == false);
   
    return Vars;
   
    } //AlgIndVarsRec(PPVector L, long codim)
  */


  vector<PPMonoidElem> AlgIndVarsRec(const vector<PPMonoidElem>& L, long codim){
    VerboseLog VERBOSE("AlgIndVarsRec: ");
    VERBOSE(40) << "L = " << L << endl;
    VERBOSE(40) << "codim = " << codim << endl;
    if (L.empty()) return L;

    PPMonoidElem uno = one(owner(L[0]));
    vector<PPMonoidElem> deadEnd;
    deadEnd.push_back(uno);
    
    if (codim == 0){
      VERBOSE(40) << "dead end" << endl;
      return deadEnd;
    }
    
    PPMonoidElem FirstOfL = First(L);
    vector<PPMonoidElem> X = WhichIndets(FirstOfL);
    for (long i = 0; i < len(X); i++){
      vector<PPMonoidElem> T;
      for (long j = 1; j < len(L); j++)
        if (!IsDivisible(L[j], X[i]))
          T.push_back(L[j]);
      //std::sort (T.begin(), T.end());
      
      if (len(T) == 0 && codim == 1)   // do we have to check codim=1?
      {
        vector<PPMonoidElem> result;
        result.push_back(X[i]);
        VERBOSE(40) << "len(T)=0: " << result << endl;
        return result;                               
      }                                              
      
      vector<PPMonoidElem> RecResult = AlgIndVarsRec(T, codim-1);
      if (RecResult[0] != uno){
        RecResult.push_back(X[i]);
        VERBOSE(40) << "len(T)>0: " << RecResult << endl;
        return RecResult;
      }
      /*
        if (RecResult[0] == impossible){
        VERBOSE(40) << "RecResult[0] = impossible" << endl;
        return RecResult;
        }
        RecResult.push_back(X[i]);
        VERBOSE(40) << "len(T)>0: " << RecResult << endl;
        return RecResult;
      */
    }
    VERBOSE(40) << "dead end" << endl;
    return  deadEnd;
  } //AlgIndVarsRec(PPVector L, long codim)
  
  
  /*
    vector<PPMonoidElem> AlgIndVarsRec(const PPVector& L1, long codim){
    vector<PPMonoidElem> L;
    if (IsEmpty(L1)) return L;
    for (long i=0; i<len(L1); ++i)  L.push_back(PP(L1[i]));


    PPMonoidElem impossible = one(owner(L[0]));
    vector<PPMonoidElem> impossibleList;
    impossibleList.push_back(impossible);


    PPMonoidElem FirstOfL = First(L);
    vector<PPMonoidElem> TailOfL = Tail(L);
    
    if (codim <= 0){
    cout << "AlgIndVarsRec -- codim <= 0" << endl;
    //CoCoA_THROW_ERROR1(ERR::BadArg);
    return impossibleList; // just to keep the compiler quiet
    }


    vector<PPMonoidElem> result;
    for (long i = 0; i < len(WhichIndets(FirstOfL)); i++){
    vector<PPMonoidElem> T;
    for (long j = 1; j < len(L); j++){
	
    if (!IsDivisible(L[j], WhichIndets(FirstOfL)[i])){
	  T.push_back(L[j]);
    }
    }
    if (len(T) != 0){
    // convert vec<PPMonoidElem> to  PPVector 
    const PPMonoid PPM1 = owner(T[0]);
    const DivMaskRule DMR1 = NewDivMaskEvenPowers();
    PPVector TPPVec(PPM1, DMR1);
    for (long k = 0; k < len(T); k++){
    PushBack(TPPVec, T[k]);
    }
    vector<PPMonoidElem> Indets;
    vector<PPMonoidElem> tmp = AlgIndVarsRec(TPPVec, codim-1);
    for (long l = 0; l < len(tmp); ++l){
    Indets.push_back(tmp[l]);
    }
    //      if (len(tmp) == 1 && tmp[0] == impossible){
    result.push_back(WhichIndets(FirstOfL)[i]);
    for (long l = 0; l < len(tmp); ++l){
	  result.push_back(tmp[l]);
    }
    }
    else{
    result.push_back(WhichIndets(FirstOfL)[i]);
    }
    //      }
    }
    return  result; 
    } //AlgIndVarsRec(PPVector L, long codim)
  
  */

  /*
    vector<PPMonoidElem> AlgIndVars(vector<PPMonoidElem> L, long codim){
    //cout << "L: " << L << endl;
    //first convert vec<PPMonidElem> to PPVector
    VerboseLog VERBOSE("AlgIndVars: ");
    VERBOSE(40) << "L: " << L << endl;
    const PPMonoid PPM1 = owner(First(L));
    const DivMaskRule DMR1 = NewDivMaskEvenPowers();
    PPVector LAsPPVec(PPM1, DMR1);
    
    for (long i = 0; i < len(L); ++i){
    PushBack(LAsPPVec, radical(L[i]));
    }
   
    //cout << "LAsPPVec: " << LAsPPVec << endl;
    interreduce(LAsPPVec); //interreduce(ref L);


    vector<PPMonoidElem> LInter;
    for (long j = 0; j < len(LAsPPVec); j++){
    LInter.push_back(PP(LAsPPVec[j]));
    }
    //cout << "Linter: " << LInter << endl;


    vector<PPMonoidElem> Inds;
    Inds = AlgIndVarsRec(LInter, codim);

    VERBOSE(40) << "codim: " << codim << endl;
    VERBOSE(40) << "Inds: " << Inds << endl;

    //    cout << "Inds in AlgIndVars = " << Inds << endl;
    if (len(Inds) != codim){
    //cout << "AlgIndVars -- len(Inds) != codim " << endl;
    CoCoA_THROW_ERROR2(ERR::BadArraySize, "AlgIndVars; len(Inds) != codim");
    return L; // just to keep compiler quiet
    }

    
    //if ((len(Inds) != codim) && (len(Inds) != codim - 1)){
    // CoCoA_THROW_ERROR2(ERR::BadArraySize, "AlgIndVars; len(Inds) != codim && len(Inds) != codim-1");
    // return L; // just to keep compiler quiet
    //}
    

    //std::sort(Inds.begin(), Inds.end());
    return Inds;
    } //AlgIndVars(vec<PPMonoidElem> L, long codim)
  */

  vector<PPMonoidElem> AlgIndVars(const vector<PPMonoidElem>& L, long codim){
    //first convert vec<PPMonidElem> to PPVector
    VerboseLog VERBOSE("AlgIndVars: ");
    VERBOSE(40) << "L: " << L << endl;
    const PPMonoid PPM1 = owner(First(L));
    const DivMaskRule DMR1 = NewDivMaskEvenPowers();
    PPVector LAsPPVec(PPM1, DMR1);
    
    for (long i = 0; i < len(L); ++i)  PushBack(LAsPPVec, radical(L[i]));
   
    interreduce(LAsPPVec); //interreduce(ref L);

    vector<PPMonoidElem> LInter;
    vector<PPMonoidElem> Inds;
    for (long j = 0; j < len(LAsPPVec); j++){
      if (IsIndet(PP(LAsPPVec[j]))){
        Inds.push_back(PP(LAsPPVec[j]));
      }
      else{
        LInter.push_back(PP(LAsPPVec[j]));
      }
    }
    
    VERBOSE(40) << "LInter: " << LInter << endl;
    VERBOSE(40) << "Inds: " << Inds << endl;
    
    vector<PPMonoidElem> Rest;
    Rest = AlgIndVarsRec(LInter, codim-len(Inds));
    for (long i=0; i<len(Rest); ++i){
      Inds.push_back(Rest[i]);
    }

    VERBOSE(40) << "codim: " << codim << endl;
    VERBOSE(40) << "Inds: " << Inds << endl;

    if (len(Inds) != codim)  CoCoA_THROW_ERROR2(ERR::BadArg, "len(Inds) != codim");

    //std::sort(Inds.begin(), Inds.end());
    return Inds;
  } //AlgIndVars(vec<PPMonoidElem> L, long codim)
  
  
  RingElem T2TT(const ring& Kx, const RingElem& T, const vector<PPMonoidElem>& TrueInds)
  {
    vector<long> L = exponents(T);
    vector<long> L1(NumIndets(Kx), 0);
    long index;
    for (long i = 0; i < len(TrueInds); ++i)
      if (IsIndet(index, TrueInds[i]))  L1[index] = L[index];
    return MakeTerm(Kx, L1);
  } // T2TT(Kx, T, TrueInds)
   
  RingElem T2FT(const ring& Kx, const RingElem& T, const vector<PPMonoidElem>& TrueInds){
    vector<long> L = exponents(T);
    long index;
    for (long i = 0; i < len(TrueInds); ++i){
      //L[IndetIndex(TrueInds[i])] = 0;
      if(IsIndet(index, TrueInds[i])) L[index] = 0;
    }
    return MakeTerm(Kx, L);
  } // T2FT(Kx, T, TrueInds)
  
 
  vector<RingElem> LeadTerms(const vector<RingElem>& L, const vector<PPMonoidElem>& Inds){
    VerboseLog VERBOSE("LeadTerms: ");
    ring Kx = RingOf(ideal(L[0]));
    ring CR = CoeffRing(Kx);
    ring NewRing = NewPolyRing(CR, symbols(Kx), NewMatrixOrdering(MakeBlockOrder(Inds), 0));
    RingHom Old2New = PolyAlgebraHom(Kx, NewRing, indets(NewRing));
    vector<RingElem> image = Old2New(L);
    vector<RingElem> K;
    for (long i = 0; i < len(L); ++i){
      K.push_back(LT(image[i]));
    }
    RingHom New2Old = PolyAlgebraHom(NewRing, Kx, indets(Kx));
    K = New2Old(K);
    VERBOSE(50) << "K: " << K << endl;
    return K;
  } // LeadTerms(L, Inds)


  vector<RingElem> Extr(vector<RingElem> G, vector<PPMonoidElem> Inds){
    ring Kx = RingOf(ideal(G[0]));
    ring CR = CoeffRing(Kx);
    ring NewRing = NewPolyRing(CR, symbols(Kx), NewMatrixOrdering(MakeBlockOrder(Inds), 0));
    RingHom Old2New = PolyAlgebraHom(Kx, NewRing, indets(NewRing));
    G = Old2New(G);
    //cout << "G in Extr nachher = " << G << endl;
    for (long i = 0; i < len(G); ++i){
      G[i] = monic(G[i]);
    }
    vector<RingElem> G1;


    while (!G.empty()){
      G1.push_back(First(G)); //instead of append(ref G1, first(G));
      RingElem H = T2TT(NewRing, LT(First(G)), Inds);
      vector<RingElem> HInList;
      HInList.push_back(H);
      //G:=[P In tail(G) | NR(Rad.T2TT(NewRing, LT(P),Inds),[H])<>0];
      for (long i = 0; i < len(Tail(G)); ++i){
        if (NR(T2TT(NewRing, LT(Tail(G)[i]), Inds), HInList) != 0){
          G[i] = Tail(G)[i];
        }
      } //here G filled
      G = Tail(G);
    } // end while


    vector<RingElem> K;
    for (long i = 0; i < len(G1); ++i){
      RingElem Q = zero(NewRing); 
      RingElem T = T2TT(NewRing, LT(G1[i]), Inds);
      RingHom phi = CoeffEmbeddingHom(NewRing);


      while ((G1[i] != (G1[i]-G1[i])) && (T2TT(NewRing, LT(G1[i]), Inds) == T)){
        Q = Q + phi(LC(G1[i]))*T2FT(NewRing, LT(G1[i]), Inds);
        G1[i] = G1[i] - phi(LC(G1[i]))*LT(G1[i]);
      }
      K.push_back(Q);
    }
    RingHom New2Old = PolyAlgebraHom(NewRing, Kx, indets(Kx));
    K = New2Old(K);
    return K;
  } //Extr(G, Inds)


  ideal ElimGB(const ideal& I, const vector<PPMonoidElem>& Inds){
    ring Kx = RingOf(I);
    int GD = GradingDim(Kx);
    if (GD > 0 && IsHomog(gens(I))){
      vector<long> index;
      long index_tmp;
      for (long i = 0; i < len(Inds); ++i){
        if(IsIndet(index_tmp, Inds[i])) index.push_back(index_tmp);
      }
      matrix M = ElimHomogMat(index, GradingMat(Kx));
      PolyRing R_tmp = NewPolyRing(CoeffRing(Kx), symbols(Kx), NewMatrixOrdering(M, GD));


      RingHom phi = PolyAlgebraHom(R_tmp, Kx, indets(Kx));
      RingHom psi = PolyAlgebraHom(Kx, R_tmp, indets(R_tmp));
      vector<RingElem> GB = GBasis(ideal(psi(gens(I))));
      return ideal(phi(GB));
    } // end-if
    else{
      vector<long> index;
      long index_tmp;
      for (long i = 0; i < len(Inds); ++i){
        if(IsIndet(index_tmp, Inds[i])) index.push_back(index_tmp);
      }
      matrix M = ElimMat(index, GradingMat(Kx));
      PolyRing R_tmp = NewPolyRing(CoeffRing(Kx), symbols(Kx), NewMatrixOrdering(M, 1));


      RingHom phi = PolyAlgebraHom(R_tmp, Kx, indets(Kx));
      RingHom psi = PolyAlgebraHom(Kx, R_tmp, indets(R_tmp));
      vector<RingElem> GB = GBasis(ideal(psi(gens(I))));
      return ideal(phi(GB));
    } // end-else
    // last 4 rows
    // vector<RingElem> GB = GBasis(ideal(BringIn(R_tmp, gens(I)))); //BRINGIN??
    // return ideal(BringIn(Kx, GB));
  } // ElimGB(I, Inds)


  vector<vector<PPMonoidElem> > AlgIndVarsListRec(const vector<PPMonoidElem>& L, long codim){
    vector<vector<PPMonoidElem> > result;
    if (L.empty()){
      return result;
    }
    if (codim == 0)
    {
      //  CoCoA_THROW_ERROR(ERR::BadArg, "AlgIndVarsListRec; codim = 0");
      result.push_back(L);
      return result;
    }
    vector<PPMonoidElem> K = WhichIndets(First(L));
    vector<vector<PPMonoidElem> > V;
    while (true){
      PPMonoidElem X = First(K);
      vector<PPMonoidElem> XAsList;
      XAsList.push_back(X);
      K = Tail(K);
      
      vector<PPMonoidElem> L1;
      for (long i = 0; i < len(L); ++i){
        // fill L1
        if (!IsDivisible(L[i], X)) L1.push_back(L[i]);
      }
    
      vector<vector<PPMonoidElem> > V1 = AlgIndVarsListRec(L1, codim-1);
      if (V1.empty()){
        V.push_back(XAsList);
      }
      else{
        for (long i=0; i<len(V1); ++i){
          V1[i].push_back(X);
          V.push_back(V1[i]);
        }
      }
      if (K.empty()) break; //until
    } // end-while
    for (long i=0; i < len(V); ++i){
      if (len(V[i]) == codim) result.push_back(V[i]);
    }
    return result;
  } //AlgIndVarsListRec(L, codim)


  vector<vector<PPMonoidElem> > AlgIndVarsList(vector<PPMonoidElem> L1, long codim, const long& N){
    vector<vector<PPMonoidElem> > result;
    if (L1.empty()) return result;
    if (codim == 0)
    {
      CoCoA_THROW_ERROR2(ERR::BadArg, "codim == 0");
      return result; //just to keep the compiler quiet
    }
    
    for (long i=0; i < len(L1); ++i){
      L1[i] = radical(L1[i]);
    }


    const PPMonoid PPM1 = owner(First(L1));
    const DivMaskRule DMR1 = NewDivMaskEvenPowers();
    PPVector LAsPPVec(PPM1, DMR1);
    for (long i = 0; i < len(L1); ++i){
      PushBack(LAsPPVec, L1[i]);
    }
    interreduce(LAsPPVec);
    
    vector<PPMonoidElem> L;
    for (long i = 0; i < len(LAsPPVec); ++i){
      L.push_back(PP(LAsPPVec[i]));
    }
    
    vector<PPMonoidElem> K = WhichIndets(First(L));
    
    vector<vector<PPMonoidElem> > V;
    while (true){
      PPMonoidElem X = First(K);
      vector<PPMonoidElem> XInList;
      XInList.push_back(X);
      K = Tail(K);
      vector<PPMonoidElem> L1;
      for(long i=0; i < len(L); ++i){
        if (!IsDivisible(L[i], X)) L1.push_back(L[i]);
      }
      vector<vector<PPMonoidElem> > V1 = AlgIndVarsListRec(L1, codim-1);
      if (V1.empty()){
        V.push_back(XInList);
      }
      else{
        for (long i=0; i<len(V1); ++i){
          V1[i].push_back(X);
          V.push_back(V1[i]);
        }
      }
      if (K.empty() || len(V) >= N) break;
    }
    for (long i = 0; i < len(V); ++i){
      if (len(V[i]) == codim) result.push_back(V[i]);
    }
    return result;
  } // AlgIndVarsList(L, codim, N)


  vector<vector<PPMonoidElem> > MakeAdjVarLists(vector<PPMonoidElem> Inds,const long& N){
    VerboseLog VERBOSE("MakeAdjVarLists: ");
    vector<vector<PPMonoidElem> > L;
    vector<PPMonoidElem> firstN;
    vector<PPMonoidElem> rest;
    VERBOSE(60) << "Inds: " << Inds << endl;
    VERBOSE(60) << "N: " << N << endl;
    while (len(Inds) >= N){
      for (long i=0; i < N; ++i){
        // first N elements
        firstN.push_back(Inds[i]); //append(ref L, first(Inds, N))
      }
      //sort(firstN.begin(), firstN.end());
      L.push_back(firstN);
      Inds = Tail(Inds);
    } // end-while
    VERBOSE(60) << "after while" << endl;
    //append(L, concat(Inds,[L[1,1]]));
    Inds.push_back(L[0][0]);
    L.push_back(Inds);
    return L;
  } // MakeAdjVarLists(Inds, N)


  // important for MakeIndsList
  vector<PPMonoidElem> First(const vector<PPMonoidElem>& P, const long& J){
    vector<PPMonoidElem> result;
    for (long i = 0; i < J; ++i){
      result.push_back(P[i]);
    }
    return result;
  } // First(P, J)


  // important for MakeIndsList
  vector<PPMonoidElem> Last(const vector<PPMonoidElem>& P, const long& J){
    vector<PPMonoidElem> result;
    for (long i = len(P) - J; i < len(P); ++i){
      result.push_back(P[i]);
    }
    return result;
  } // Last(P, J)


  vector<vector<PPMonoidElem> > MakeSet(vector<vector<PPMonoidElem> > P)
  {
    vector<vector<PPMonoidElem> > result;
    for (const auto& v:P)  if (!IsIn(v, result))  result.push_back(v);
    // while (!P.empty()){
    //   vector<PPMonoidElem> firstElem = First(P);
    //   vector<vector<PPMonoidElem> > tmp;
    //   tmp.push_back(firstElem);
    //   P = Tail(P);
    //   if (!IsSubset(tmp, P)) result.push_back(firstElem);
    //   firstElem.clear();
    //   tmp.clear();
    // }
    /*for (long i = 0; i < len(result); i++){
      std::sort(result[i].begin(), result[i].end());
      }*/
    return result;
  }

  /*
    vector<PPMonoidElem> MakeSet(vector<PPMonoidElem> P){
    vector<PPMonoidElem> result;
    while (!P.empty()){
    PPMonoidElem firstElem = First(P);
    vector<PPMonoidElem> tmp;
    tmp.push_back(firstElem);
    P = Tail(P);
    if (!IsSubset(tmp, P)) result.push_back(firstElem);
    //firstElem.clear();
    tmp.clear();
    }
    //std::sort(result.begin(), result.end());
    return result;
    }
  */  

  
  vector<vector<PPMonoidElem> > MakeIndsList(const vector<vector<PPMonoidElem> >& VList, const vector<PPMonoidElem>& Frozen, long CDim){
    VerboseLog VERBOSE("MakeIndsList: ");
    VERBOSE(50) << "VList: " << VList << endl;
    VERBOSE(50) << "Frozen: " << Frozen << endl;
    VERBOSE(50) << "CDim: " << CDim << endl;
    
    vector<vector<PPMonoidElem> > K;
    for (long i = 0; i < len(VList); i++){
      K.push_back(DIFF(VList[i], Frozen));
    }
    vector<vector<PPMonoidElem> > K2;
    for (long i = 0; i< len(K); ++i){
      if (len(K[i]) == CDim || len(K[i]) == CDim-1){
        //sort(K[i].begin(), K[i].end());
        K2.push_back(K[i]);
      }
    }
    vector<vector<PPMonoidElem> > LL;
    for (long i = 0; i < len(K); ++i){
    	if (len(K[i]) == CDim){
        for (long J = 0; J < len(K[i]); J++){
          vector<PPMonoidElem> tmp;
          tmp.push_back(K[i][J]);
          LL.push_back(DIFF(K[i], tmp));
        }
      }
      if(len(K[i]) == CDim-1){
        LL.push_back(K[i]);
      }
    }
    // mal versuchen zu sorten
    /*for (long m = 0; m < len(LL); ++m){
      sort(LL[m].begin(), LL[m].end());
      }*/
    return MakeSet(LL);
  } // MakeIndsList(VList, Frozen, CDim)


  vector<vector<PPMonoidElem> > MakeVarsList(const vector<PPMonoidElem>& IM, long CDim, const vector<PPMonoidElem>& Inds){
    VerboseLog VERBOSE("MakeVarsList: ");
    vector<vector<PPMonoidElem> > NAIVarLists = AlgIndVarsList(IM, CDim, ConvertTo<long>(binomial(ConvertTo<BigInt>(len(Inds)), ConvertTo<BigInt>(CDim))));
    vector<PPMonoidElem>  tmp;
    vector<vector<PPMonoidElem> > L1 = MakeIndsList(NAIVarLists, tmp, CDim);
    VERBOSE(60) << "L1 = " << L1 << endl;
    //sorted(L1); implizit in MakeIndsList
    for(long i = 0; i < len(L1); i++){
      sort(L1[i].begin(), L1[i].end());
    }
    VERBOSE(60) << "Inds: " << Inds << endl;
    VERBOSE(60) << "CDim-1: " << CDim-1 << endl;
    vector<vector<PPMonoidElem> > L2 = MakeAdjVarLists(Inds, CDim-1);
    //sorted(L2); implizit in MakeAdjVarLists
    VERBOSE(60) << "L2: " << L2 << endl;
    for(long i = 0; i < len(L2); i++){
      sort(L2[i].begin(), L2[i].end());
    }
 
    vector<vector<PPMonoidElem> > L;
    for (long i = 0; i < len(L1); ++i){
      L.push_back(L1[i]); // MakeSet
    }
    for (long i = 0; i < len(L2); ++i){
      L.push_back(L2[i]); // MakeSet
    }
    L = MakeSet(L);
    // L filled
    vector<PPMonoidElem> FrozenVars;
    vector<vector<PPMonoidElem> > result;
    vector<PPMonoidElem> VarList;
    while (!L.empty() && !EqSet(FrozenVars, Inds)){
      VarList = First(L);
      L = Tail(L);
      result.push_back(VarList);
      vector<PPMonoidElem> tmp2 = DIFF(Inds, VarList);
      for (long i = 0; i < len(tmp2); ++i){
        FrozenVars.push_back(tmp2[i]); //MakeSet(concat(FrozenVars,diff(Inds,VarList)));
        FrozenVars = MakeSet(FrozenVars);
      }
    }
    return result;
  } // MakeVarsList(IM, CDim, Inds)

  
  RingElem FindL(const ideal& I){
    VerboseLog VERBOSE("FindL: ");
    RingElem result;
    vector<RingElem> H = GBasis(I);
    long codim = NumIndets(RingOf(I)) - DimQuot(LT(I));
    //if (codim == 1) return gens(I);
    //check this always before calling FindL!!
    vector<vector<PPMonoidElem> > ToDoVarLists;
    // converting gg = gens(LT(I)) in vec<PPMonoidElem>
    vector<RingElem> gg = gens(LT(I));
    const PPMonoid PPM1 = PPM(RingOf(ideal(First(gg))));
    const DivMaskRule DMR1 = NewDivMaskEvenPowers();
    PPVector gg2(PPM1, DMR1);
    convert(gg2, gg); // gg PPVector
    vector<PPMonoidElem> gensOfLTI;
    for (long i = 0; i < len(gg2); ++i)  gensOfLTI.push_back(PP(gg2[i]));
    ToDoVarLists = MakeVarsList(gensOfLTI, codim, WhichIndets(LT(I)));

    while (!ToDoVarLists.empty())
    {
      VERBOSE(30) << "Main Cycle - F " << len(ToDoVarLists) << endl;
      vector<PPMonoidElem>  GBInds;
      GBInds = First(ToDoVarLists);  
      ToDoVarLists = Tail(ToDoVarLists);
      ideal K1 = ElimGB(I, GBInds);
      vector<RingElem> K = gens(K1);
      vector<RingElem> L = PutPolyLInSubring(K, GBInds);
      if (len(L) > 1)
      {
        VERBOSE(30) << "++++++++++++++++++++++++++ len(L1) > 1" << endl;
        result = FindL(ideal(L));
      }
      if (len(L) == 1)  result = First(L);
      return result;
    } // end-while
    VERBOSE(30) << "ERROR" << endl;
    CoCoA_THROW_ERROR1(ERR::ShouldNeverGetHere); // usually impossible
    return result; //just to keep the compiler quiet
  } // FindL(I)


  vector<ideal> ElimRedundant(vector<ideal> E)
  {
    for (long i=0; i < len(E); ++i)
      for (long j=i+1; j < len(E); ++j)
      {
        if (IsContained(E[j], E[i]))  E[i] = ideal(one(RingOf(E[i])));
        else
          if (IsContained(E[i], E[j]))  E[j] = ideal(one(RingOf(E[j]))); 
      }
    vector<ideal> result;
    for (long i=0; i < len(E); ++i)
      if (!IsElem(one(RingOf(E[i])), E[i]))  result.push_back(E[i]);
    return result;
  } // ElimRedundant(E)


  // important for TrueEquiDec
  RingElem gcd(const vector<RingElem>& g)
  {
    RingElem res = First(g);
    for (long i = 0; i < len(g); ++i)  res = gcd(res, g[i]);
    return res;
  } // gcd(g)



  // ----------- Big Function: TrueEquiDec --------------------


  vector<ideal> TrueEquiDec(ideal I, const long& V)
  {
    VerboseLog VERBOSE("TrueEquiDec: ");
    vector<ideal> result;
    
    ring Kx = RingOf(I);
    //    cout << "TrueEquiDec I = " << I << endl;
    I = MonomialCleaning(I);
    VERBOSE(30) << "-- after MonomialCleaning" << endl;

    //    cout << "TrueEquiDec nach MonomialCleaning I = " << I << endl;
    
    if (IsOne(I))
    {
      RADICAL_BRANCHING--;
      VERBOSE(30) << "-- ideal(1) " << RADICAL_BRANCHING << endl;
      result.push_back(I);
      return result;
    }
    //    cout << "hier" << endl;
    if (AreGensMonomial(I))
    {
      RADICAL_BRANCHING--;
      long D;
      //      cout << "hier2" << endl;
      ///      if (IsHomog(I)) {D = DimQuot(I); //instead of dim(Kx/I);
        //	cout << "hier3" << endl;
      ///      }
      ///      else
      D = DimQuot(LT(I));
      VERBOSE(30) << "-- Mon-Dim Pura = " << D << RADICAL_BRANCHING << endl;
      //      cout << "hier4" << endl;
      result.push_back(SqFreeMonId(I));
      //      cout << "hier5" << endl;
      return result;
    }
    // cout << "hier6" << endl;
    if (len(gens(I)) == 1)
    {
      RADICAL_BRANCHING--;
      long D;
      //      if (IsHomog(I)) D = DimQuot(I); //instead of dim(Kx/I);
      //      else
      D = DimQuot(LT(I));
      VERBOSE(30) << "-- Principale-Dim Pura = " << D << RADICAL_BRANCHING << endl;
      result.push_back(I);
      return result;
    }

    vector<RingElem> H = GBasis(I);
    
    if (AreGensSqFreeMonomial(LT(I)))
    {
      RADICAL_BRANCHING--;
      VERBOSE(30) << "-- LT(I) Radicale " << RADICAL_BRANCHING << endl;
      result.push_back(I);
      return result;
    }

    long D;
    //    if (IsHomog(I)) {D = DimQuot(I);} //instead of dim(Kx/I);
    //    else 
    D = DimQuot(LT(I));
    long CDim = NumIndets(Kx) - D;

    if (CDim > min(V, len(gens(I))))
    {
      RADICAL_BRANCHING--;
      VERBOSE(30) << "-- CDim > Ht " << RADICAL_BRANCHING << endl;
      result.push_back(ideal(one(Kx)));
      return result;
    }

    vector<PPMonoidElem> Inds;

    // converting gg = gens(LT(I)) in vec<PPMonoidElem>
    vector<RingElem> gg = gens(LT(I));
    
    vector<PPMonoidElem> gensOfLTI;
    for (long i = 0; i < len(gg); ++i)  gensOfLTI.push_back(LPP(gg[i]));

    Inds = AlgIndVars(gensOfLTI, CDim); // PPMonoidElem needed
    VERBOSE(50) << "Inds after AlgIndVars: " << Inds << endl;
    
    vector<PPMonoidElem>  emptyList;
    vector<vector<PPMonoidElem> > IndsList;
    IndsList.push_back(Inds);

    Inds = First(MakeIndsList(IndsList,emptyList,CDim)); //List of Inds needed
    VERBOSE(50) << "Inds after MakeIndsList: " << Inds << endl;
    
    if (D - NumIndets(Kx) + len(WhichIndets(I)) == 0)
    {
      RADICAL_BRANCHING--;
      VERBOSE(30) << "-- 0-Dim (local) " << RADICAL_BRANCHING << endl;
      result.push_back(I);
      return result;
    }

    vector<RingElem> K; // to be seen
    if (CDim != 1)
    {
      VERBOSE(30) << "-- codim != 1 " << RADICAL_BRANCHING << endl;
      vector<RingElem> tmp = gens(ElimGB(I,Inds));
      for (long i=0; i < len(tmp); i++)  K.push_back(tmp[i]);
      VERBOSE(30) << " after ElimGB " << RADICAL_BRANCHING << endl;
    }
    else{
      RADICAL_BRANCHING--;
      VERBOSE(30) << "codim = 1 Splitting " << RADICAL_BRANCHING << endl;
      ideal I1 = ideal(gcd(gens(I)));
      ideal I2 = saturate(I, I1);

      vector<ideal> I1V = TrueEquiDec(I1, V);
      vector<ideal> I2V = TrueEquiDec(I2, V);
      vector<ideal> res = I1V;
      for (long i = 0; i < len(I2V); ++i)  res.push_back(I2V[i]);
      
      return res; // Return concat(TrueEquiDec(I1,V), TrueEquiDec(I2,V));
    }

    if (AreGensSqFreeMonomial(ideal(LeadTerms(K, Inds))))
    {
      RADICAL_BRANCHING--;
      VERBOSE(30) << "-- Lt(I) Radicale " << RADICAL_BRANCHING << endl;
      result.push_back(I);
      return result;
    }
    vector<RingElem> L = PutPolyLInSubring(K, Inds);
    // vector<RingElem> L = PutPolyInSubring(gens(First(K)),Inds);
    //    cout << "TrueEquiDec L = " << L << " K= " << K << endl;
    
    VERBOSE(30) << "-- after PutPolyLInSubring" << endl;

    RingElem G; //to be seen
    vector<ideal> Kideal;
    if (len(L) > 1)
    {
      VERBOSE(30) << "GCD splitting " << endl;
      G = gcd(L);
      // if(IsOne(ideal(G))){return result;}
      //      cout << "TrueEquiDec ggt = " << G << endl;
      if (RADICAL_FULLGCDSPLIT){
        VERBOSE(39) << "-- Should not happen" << endl;
        // is always false
        //RADICAL_BRANCHING++;
        //VERBOSE(30) << "[ " << len(L) << " ]" << RADICAL_BRANCHING << endl;
        //RingElem K1 = radical(G);
        //RingElem K2 = radical(P/G);   //P never defined????? 
        //K =concat([I+ideal(K1)],[I+ideal(K2)|PP In L]); WHAT IS THIS??
      }
      else{
        RADICAL_BRANCHING++;
        if (RADICAL_SPLITTINGGCD){
          VERBOSE(30) << RADICAL_BRANCHING << endl;
          RingElem K1 = radical(G);
          Kideal.clear();
          Kideal = Splitting(I, K1);
        }
        else{
          VERBOSE(30) << " Factoring " << RADICAL_BRANCHING << endl;
          RingElem K1 = radical(G);
          RingElem K2 = radical(L[0]/G);
          Kideal.clear();
          Kideal.push_back(I+ideal(K1));
          Kideal.push_back(I+ideal(K2));
        }
      }
      //for loop which contains all results of
      //TrueEquiDec in a list and return...
      //KKK = ConcatLists([TrueEquiDec(JJ,V) | JJ In K]);
      //Return flatten(KKK,1);

      vector<ideal> KKK;
      for (long i = 0; i<len(Kideal); ++i)
      {
        vector<ideal> tmpTrueEquiDec = TrueEquiDec(Kideal[i], V);
        for (long j = 0; j < len(tmpTrueEquiDec); ++j)
          KKK.push_back(tmpTrueEquiDec[j]);
      }
      return KKK;
    } // end-if 

    VERBOSE(30) << "len(L) <= 1: " << len(L) << endl;
    G = First(L);

// factorization<RingElem> Q = factor(one(RingOf(I))); // sichtbar machen
    factorization<RingElem> Q = SqFreeFactor(G); // redmine #1779.#3

    ////  Anna: 2024-02  what is this blabla??
    // if (IsZZ(owner(G))) {
    //   /*  cout << "blabla" << endl;
    //       Q.myAppend(G, 1);  // G^1
    //       cout << "blabla2" << endl;*/
    // }
    // else
    {
//      Q = SqFreeFactor(G);  // done above
      VERBOSE(30) << "-- after SqFreeFactor" << endl;}
    
    vector<RingElem> QFs = Q.myFactors();
    vector<ideal> Kid2;
    if (len(QFs) != 1 || Q.myMultiplicities()[0] != 1)
    {
      VERBOSE(30) << "new SqFree" << endl;

      if (RADICAL_FULLSqFreeSPLIT || IsElem(product(QFs), I))
      {
        RADICAL_BRANCHING = RADICAL_BRANCHING + len(QFs) -1;
        VERBOSE(30) << " Splitting [ " << len(QFs) << " ]" << RADICAL_BRANCHING << endl;

        //K.clear();
        for (long i=0; i < len(QFs); ++i)  Kid2.push_back(I+ideal(QFs[i]));

        //same as above
        //KKK = ConcatLists([TrueEquiDec(JJ,V) | JJ In K]); // K nochmal umfüllen
        //return flatten(KKK, 1); // FLATTEN?? nach konstruktion unnötig
        
        vector<ideal> KKK;
        for (long i = 0; i<len(Kid2); ++i)
        {
          vector<ideal> tmpTrueEquiDec = TrueEquiDec(Kid2[i], V);
          for (long j = 0; j < len(tmpTrueEquiDec); ++j)
            KKK.push_back(tmpTrueEquiDec[j]);
        }
        return KKK;
      }
      else
      {
        RingElem P = product(QFs);
        VERBOSE(30) << " in else  " << endl;
        return TrueEquiDec(I+ideal(P), V);
      }
    }

    //    cout << "nach dem" << endl;

    vector<RingElem> GI;
    //cout << "K = " << K << "L = " << L << "Inds = " << Inds << endl;
    GI = Extr(DIFF(K, L), Inds);
    //cout << "GI = " << GI << endl;
    
    VERBOSE(30) << "-- after Extr" << endl;
    //    long J = 2;
    long J = 1;
    RingElem GIGCD;
    GIGCD = gcd(G, GI[0]);
    while (IsOne(GIGCD) && J < len(GI))
    {
      GIGCD = gcd(G, GI[J]);
      J++;
    }
    
    if (!IsOne(GIGCD))
    {
      RADICAL_BRANCHING++;
      VERBOSE(30) << "Gi Splitting" << RADICAL_BRANCHING << endl;
      vector<ideal> Id = Splitting(I, GI[J-1]);


      vector<ideal> ID1V = TrueEquiDec(Id[0], V);
      VERBOSE(60) << "after ID1V" << endl;
      vector<ideal> ID2V = TrueEquiDec(Id[1], V);
      VERBOSE(60) << "after ID2V" << endl;
      vector<ideal> res = ID1V;
      for (long i = 0; i < len(ID2V); ++i){
        res.push_back(ID2V[i]);
      }
      return res; // Return concat(TrueEquiDec(Id[1],V),TrueEquiDec(Id[2],V));
    }


    RingElem G1 = radical(product(GI));
  
    double t0 = CpuTime();
    ideal I2 = saturate(I, ideal(G1));
    VERBOSE(30) << "-- after saturate: " << CpuTime() -t0 << endl;


    vector<ideal> Kid3;
    if (I2 != I){
      VERBOSE(30) << "I2 Splitting" << endl;
      if (RADICAL_FULLGISPLIT){
        //GI:=[PP In GI| PP<>1];


        vector<RingElem> PP;
        for (long i = 0; i < len(GI); ++i){
          if (!IsOne(GI[i])) PP.push_back(GI[i]); 
        }
        GI = PP; //GI from above
	
        RADICAL_BRANCHING = RADICAL_BRANCHING + len(GI) -1;
        VERBOSE(30) << "[ " << len(GI) << " ]" << RADICAL_BRANCHING;
        //K = concat([I2], flatten([[I+ideal(F)]|F In GI], 1));
        //K.clear();
        Kid3.push_back(I2);
        for (long i = 0; i < len(GI); ++i){
          Kid3.push_back(I+ideal(GI[i]));
        }
       
      }
      else{
        RADICAL_BRANCHING++;
        VERBOSE(30) << RADICAL_BRANCHING << endl;
        if (RADICAL_BISATSPLITTING){
          //K.clear();
          //K = [I2,saturate(I,I2)];
          Kid3.push_back(I2);
          Kid3.push_back(saturate(I, I2));
        }
        else{
          //K.clear();
          //K = [I2,I+ideal(G1)];
          Kid3.push_back(I2);
          Kid3.push_back(I+ideal(G1));
        }
      }
      //same as above
      //K = ConcatLists([TrueEquiDec(JJ,V) | JJ In K]);
      //Return flatten(K,1);


      vector<ideal> KKK;
      for (long i = 0; i<len(Kid3); ++i){
        vector<ideal> tmpTrueEquiDec = TrueEquiDec(Kid3[i], V);
        for (long j = 0; j < len(tmpTrueEquiDec); ++j){
          KKK.push_back(tmpTrueEquiDec[j]);
        }
      }
      return KKK;
    }


    RADICAL_BRANCHING--;
    VERBOSE(30) << "-- Finale - Dim Pura " << D << RADICAL_BRANCHING << endl;
    // l'ultimo return [I];
    result.push_back(I);
    return result;
  } // TrueEquiDec(I, V)



  
  // -----------Big Function: TrueEquiRad ---------------------


  vector<ideal> TrueEquiRad(ideal I)
  {
    VerboseLog VERBOSE("TrueEquiRad: ");
    I = MonomialCleaning(I);
    vector<ideal> result;

    if (AreGensMonomial(I))
    {
      VERBOSE(30) << "-- Monomial ideal" << endl;
      result.push_back(SqFreeMonId(I));
      return result;
    } // ideale monomiale

    if (NumGens(I) == 1)
    {
      VERBOSE(30) << "-- Principal ideal" << endl;
      result.push_back(ideal(radical(First(gens(I)))));
      return result;
    }

    // filling set P
    vector<RingElem> P;
    for (long i=0; i < len(gens(I)); ++i){
      if (gens(I)[i] != zero(RingOf(I)) && deg(gens(I)[i]) != 1){
        P.push_back(gens(I)[i]);
      }
    }
    if (P.empty()){
      VERBOSE(30) << "-- Linear ideal" << endl;
      result.push_back(I);
      return result;
    }

    if (AreGensSqFreeMonomial(LT(I)))
    {
      VERBOSE(30) << "-- Monomial ideal" << endl;
      result.push_back(I);
      return result;
    }
    long D;
    //    if (IsHomog(I)) D = DimQuot(I); //instead of dim(Kx/I);
    //    else
    D = DimQuot(LT(I));
    long codim = NumIndets(RingOf(I)) - D; //dim(RingOf(I)/I);

    if (codim == 1)
    {
      if (NumGens(I) == 1)
      {
        VERBOSE(30) << "-- Hypersurface" << endl;
        RingElem K1 = radical(First(gens(I)));
        result.push_back(ideal(K1));
        return result;
      }
      else
      {
        RingElem H = gcd(gens(I));
        RingElem K1 = radical(H);
        ideal tmp = intersect(ideal(K1), radical(saturate(I, ideal(H))));
        result.push_back(tmp);
        return result;
      }
    }

    // converting gg = gens(LT(I)) in vec<PPMonoidElem>
    vector<RingElem> gg = gens(LT(I));
    const PPMonoid PPM1 = PPM(RingOf(ideal(First(gg))));
    const DivMaskRule DMR1 = NewDivMaskEvenPowers();
    PPVector gg2(PPM1, DMR1);
    convert(gg2, gg); // ora gg PPVector
    vector<PPMonoidElem> gensOfLTI;
    for (long i = 0; i < len(gg2); ++i){
      gensOfLTI.push_back(PP(gg2[i]));
    }

    vector<vector<PPMonoidElem> > ToDoVarLists = MakeVarsList(gensOfLTI, codim, WhichIndets(I));
    VERBOSE(60) << "ToDoVarLists: " << ToDoVarLists << endl;
    
    while (!ToDoVarLists.empty()){
      VERBOSE(30) << "Main Cycle " << len(ToDoVarLists) << endl;

      vector<PPMonoidElem> GBInds = First(ToDoVarLists);
      ToDoVarLists = Tail(ToDoVarLists);
      vector<RingElem> K = gens(ElimGB(I, GBInds));
      VERBOSE(30) << "-- after ElimGB" << endl;

      if (AreGensSqFreeMonomial(ideal(LeadTerms(K, GBInds)))){
        VERBOSE(30) << "-- Lt(I) radical" << endl;
        result.push_back(I);
        return result;
      }

      VERBOSE(60) << "GBInds: " << GBInds << endl;
      VERBOSE(60) << "K: " << K << endl;
      vector<RingElem> L = PutPolyLInSubring(K, GBInds);
 
      bool ToProceed = true;
      if (len(L) > 1){
        long codimI = NumIndets(RingOf(I)) - DimQuot(LT(I));
        VERBOSE(60) << "codimI: " << codimI << endl;
        if (codimI == 1){
          ToProceed = false;
          if (NumGens(I) == 1){
            VERBOSE(30) << "-- Hypersurface" << endl;
            RingElem K1 = radical(First(gens(I)));
            vector<ideal> forReturn;
            forReturn.push_back(ideal(K1));
            return forReturn;
          }
          else{
            RingElem H = gcd(gens(I)); //because return I in Findl in this case
            RingElem K1 = radical(H);
            vector<ideal> forReturn;
            forReturn.push_back(I+intersect(ideal(K1), radical(saturate(I, ideal(H)))));
            return forReturn;
          }
        }
	
        VERBOSE(60) << "Right before FindL" << endl;
        VERBOSE(60) << "L: " << L << endl;
        VERBOSE(60) << "codim(ideal(L)): " << NumIndets(RingOf(ideal(L))) - DimQuot(LT(ideal(L))) << endl;

        // Problem: Never ending loop when codim of ideal(L) is 1. The function FindL calls MakeVarsList which calls MakeAdjVarLists.
        // In MakeAdjVarLists there is a while loop that stops when the length of a specific list is smaller than codim-1.
        // But if codim = 1, that means the while loop can never stop.
        // The following is an attempt to fix this (see CoCoA5 code as well, Line 1076).
        if ((NumIndets(RingOf(ideal(L))) - DimQuot(LT(ideal(L)))) == 1){
          ToProceed = false;
          if (NumGens(I) == 1){
            VERBOSE(30) << "-- Hypersurface" << endl;
            RingElem K1 = radical(First(gens(I)));
            vector<ideal> forReturn;
            forReturn.push_back(ideal(K1));
            return forReturn;
          }
          else{
            RingElem H = gcd(gens(ideal(L))); //because return I in Findl in this case
            RingElem K1 = radical(H);
            vector<ideal> forReturn;
            forReturn.push_back(I+intersect(ideal(K1), radical(saturate(I, ideal(H)))));
            return forReturn;
          }
        }
        // end of newly added stuff

        RingElem L1 = FindL(ideal(L));
        vector<RingElem> L;
        L.push_back(L1);
        //...................................
      }
      
      if (ToProceed && len(L) == 1){
        RingElem G = radical(First(L));
        if (NF(G, I) != 0){
          ideal J = I; // to be seen
          VERBOSE(30) << "Found Useful SqFree" << endl;
          if (false){
            vector<ideal> L1 = TrueEquiRad(I+ideal(G));
            J = First(L1);
          }
          else{
            J = radical(I+ideal(G));
          }
          result.push_back(J);
          return result;
        }
	
      }
    } // end-while


    VERBOSE(30) << "-- radical by exhaustion" << endl;
    result.push_back(I);
    return result;
  } // TrueEquiRad(I)
  


  //------------------ Jacobians by dedo -------------------


  // ideal IdealOfMinors(const matrix& M, const int& N){
  //   return ideal(minors(M, N));                                //MINORS not in the scope
  // }


  // ideal LowerJac(const ideal& I, const long& A){
  //   if (NumIndets(RingOf(I)) == A) return ideal(one(RingOf(I)));
  //   return I + IdealOfMinors(jacobian(gens(I)), NumIndets(RingOf(I))-A);
  // }


  // long LJDim(const ideal& I){
  //   if(IsOne(I)) return -1;
  //   long D;
  //   if (IsHomog(I)) D = DimQuot(I); //instead of dim(Kx/I);
  //   else D = DimQuot(LT(I));
  //   return D; //instead of dim(RingOf(I)/I)
  // }


  // ideal RadicalOfUM(const ideal& I){
  //   long D = len(QuotientBasis(I));
  //   ideal LoJ = LowerJac(I, D);
  //   if (LJDim(LoJ) < D) return I;
  //   long A = D;
  //   ideal SLJ = ideal(one(RingOf(I))); //defined here is ok??
  //   while (LJDim(LoJ) == D){
  //     SLJ = LoJ;
  //     A++;
  //     LoJ = LowerJac(I, A);
  //   }
  //   ideal NI = colon(I, SLJ); // I:SLJ ???
  //   return RadicalOfUM(NI);
  // }


  //---------------------- end --------------------------------


  vector<ideal> EquiIsoDec(const ideal& I)
  {
    ring Kx = RingOf(I);
    RADICAL_BRANCHING = 1;
    return TrueEquiDec(I, len(gens(I)));
  } // EquiIsoDec(I)


  ideal RadicalOfUnmixed(const ideal& I)
  {
    ring P = RingOf(I);
    vector<RingElem> K = gens(I);
    for (long i = 0; i < len(K); ++i)  K[i] = radical(K[i]);
//    ideal J = ideal(K);
    return First(TrueEquiRad(ideal(K)));
  } // RadicalOfUnmixed(I)



  namespace // anonymous
  {
    // Avoid computing "speculative radical" if multivariate in positive char
  RingElem radical_quick(const RingElem& f)
  {
    if (!IsZero(characteristic(owner(f))) && UnivariateIndetIndex(f) < 0)
      return f;
    return radical(f);
  }
  } // end of namespace anonymous
  
  ideal radical_general_compute(const ideal& I)
  {
    VerboseLog VERBOSE("radical_general_compute");
    ring Kx = RingOf(I);
    if (GradingDim(Kx)==0)
    {
      SparsePolyRing R = NewPolyRing(CoeffRing(Kx), NewSymbols(NumIndets(Kx)));
      RingHom phi = PolyAlgebraHom(Kx, R, indets(R));
      RingHom psi = PolyAlgebraHom(R, Kx, indets(Kx));
      return ideal(psi(gens(radical_general_compute(ideal(phi(gens(I)))))));
    }
    // vector<RingElem> K1 = gens(I);
    // for (long i = 0; i < len(K1); ++i) K1[i] = radical(K1[i]);
    // vector<RingElem> K2 = GBasis(I);
    // for (long i = 0; i < len(K2); ++i) K2[i] = radical(K2[i]);
    // // hier versuche aus jedem I ein J
    // ideal J = ideal(K1) + ideal(K2);
    vector<RingElem> K1;
    for (const auto& f:gens(I))  K1.push_back(radical_quick(f));
    for (const auto& f:GBasis(I))  K1.push_back(radical_quick(f));
    ideal J = ideal(K1);
    
    VERBOSE(30) << "after rad gens and GB" << endl;
    ring NewRing = NewPolyRing(CoeffRing(Kx), symbols(Kx), NewMatrixOrdering(OrdMat(Kx), GradingDim(Kx)));
    RADICAL_BRANCHING = 1;
   
    vector<ideal> E = TrueEquiDec(J, len(gens(J))); //
    //cout << "got past TrueEquiDec" << endl;
    VERBOSE(30) << "-- after TrueEquiDec" << endl;
    E = ElimRedundant(E);
    VERBOSE(30) << "-- after ElimRedundant" << endl;
    VERBOSE(30) << "++++ " << len(E) << "Components" << endl;
    //REPAIR
    //vector<ideal> H = First(TrueEquiRad(E)); //instead of first(TrueEquiRad(J))  J in E
    vector<ideal> H;
    for (long i = 0; i < len(E); ++i){
      H.push_back(First(TrueEquiRad(E[i])));
    }
    VERBOSE(30) << "++++ Intersecting" << endl;
    // newdefined E as an ideal?
    //ideal IntersecH := IntersectionList(H);
    ideal IntersecH = First(H);
    for (long i = 1; i < len(H); ++i){
      IntersecH = intersect(IntersecH, H[i]);
    }
    //not needed
    //IntersecH := image(E,RMap(indets(Kx)));
    //RingHom homom = PolyAlgebraHom(Kx, NewRing, indets(NewRing));
    //RingHom homom = PolyAlgebraHom(NewRing, Kx, indets(Kx)); //das gehts nicht
    //vector<RingElem> gensIntersecH = apply(homom, gens(IntersecH));
    
    //return ideal(gensIntersecH);
    return IntersecH;  //A// x
  } // radical_general_compute(I)
  


  

  //}  // ^^^^^ anonymous namespace ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  ideal radical(const ideal& I)
  {
    VerboseLog VERBOSE("radical");
    VERBOSE(1000) << " starting " << std::endl;
    const bool OverField = IsField(CoeffRing(RingOf(I)));
    if (IsTrue3(IsRadical3(I))) return I;
    if (IsZero(I)) return I; // ideal(Kx,[]);  // to avoid having a zero gen
    if (OverField && AreGensMonomial(I)) return radical_MonId(I); // monomial ideal
    if (IsOne(I)) return ideal(one(RingOf(I)));
    if (IsZeroDim(I)) return radical_0dim(I); // 0-dim ideal
    ideal radI = radical_general_compute(I);
    const auto ptrI = SparsePolyRingBase::IdealImpl::ourGetPtr(radI);
    ptrI->myAssignRadicalFlag(true);
    return radI;
  } // radical



} // end of namespace CoCoA
