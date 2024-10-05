//   Copyright (c)  2007  John Abbott,  Anna M. Bigatti

//   This file is part of the source of CoCoALib, the CoCoA Library.

//   CoCoALib is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.

//   CoCoALib is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.

//   You should have received a copy of the GNU General Public License
//   along with CoCoALib.  If not, see <http://www.gnu.org/licenses/>.


// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// [JAA 20141029] What the devil is this file???
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#include "CoCoA/library.H"


#include <iostream>
using std::cout;
using std::cout;
using std::endl;
using std::ostream;
#include <vector>
using std::vector;
using std::pair;
#include <string>
using std::string;
#include <functional>
#include <algorithm>
using std::bind2nd;
using std::not_equal_to;

//----------------------------------------------------------------------
// Test for Dynamic Buchberger Algorithm 
//----------------------------------------------------------------------

namespace CoCoA
{

/*
  Things to do: 
  passing the_OldPP as ideal is there a simplification?
  Use the colon of Old:*it, bigger in lieu of Old+*it, lesser
  understand the use of function adapters et al to fully exploit the STL
  function objects with class to have a status
  smarter way to erase elements of the_CurrentPP:
  use lists, not vectors
  every time there is a comparison, 
  if you find lesser, erase all previous
  if you find bigger, erase it
*/

///////////////////////////  UTILITIES

  void AssignMatrix(matrix& a,const matrix& b)
  {
    CoCoA_ASSERT_ALWAYS(RingOf(a) == RingOf(b));
    CoCoA_ASSERT_ALWAYS(NumCols(a) == NumCols(b));
    CoCoA_ASSERT_ALWAYS(NumRows(a) == NumRows(b));

    for (unsigned int i = 0;i != NumRows(a);++i)
      for (unsigned int j = 0;j != NumCols(a);++j)
        SetEntry(a,i,j,b(i,j));
  }


  // returns the support of the poly f
  vector<PPMonoidElem> supp(ConstRefRingElem the_f)
  {
    vector<PPMonoidElem> v;
    for (SparsePolyIter it = BeginIter(the_f);!IsEnded(it);++it)
      v.push_back(PP(it));
    return v;
  }//supp


  // returns the coeffs of the poly f
  vector<RingElem> coefficients(ConstRefRingElem the_f)
  {
    vector<RingElem> v;
    for (SparsePolyIter it = BeginIter(the_f);!IsEnded(it);++it)
      v.push_back(coeff(it));
    return v;
  }//supp


// ??? copied from GReductor.C
  bool BoolCmpLPP(ConstRefRingElem f, ConstRefRingElem g)
  {
    CoCoA_ASSERT_ALWAYS(owner(f) == owner(g));//BUG HUNTING  ???
    CoCoA_ASSERT_ALWAYS(!IsZero(f));//BUG HUNTING  ???
    CoCoA_ASSERT_ALWAYS(!IsZero(g));//BUG HUNTING  ???
  
    const SparsePolyRing& P  =  owner(f);
    return P->myCmpLPP(raw(f), raw(g)) > 0;
  }


// ??? copied from ex-QuotientBasis.C
  bool IsDivisibleL(ConstRefPPMonoidElem pp, const vector<PPMonoidElem>& ByL)
  {
    const long n  =  len(ByL);
    for (long i = 0; i < n; ++i)
      if ( IsDivisible(pp, ByL[i]) ) return true;
    return false;
  }//IsDivisible


  bool IsStrictlyDivisible(ConstRefPPMonoidElem pp, const vector<PPMonoidElem>& ByL)
  {
    const long n  =  len(ByL);
    for (long i = 0; i < n; ++i)
      if ( (pp != ByL[i]) && IsDivisible(pp, ByL[i]) ) return true;
    return false;
  }//IsStrictlyDivisible


  ideal interreduce(const ideal& I)
  {
    const SparsePolyRing& P  =  RingOf(I);
    vector<RingElem> g  =  gens(I);
    vector<PPMonoidElem> p;
    sort(g.begin(), g.end(), BoolCmpLPP);
    for (vector<RingElem>::reverse_iterator it = g.rbegin(); it != g.rend() ; ++it)
      if ( !IsDivisibleL(LPP(*it), p) )
        back_inserter(p)  =  LPP(*it);
    g.clear();
    for (vector<PPMonoidElem>::const_iterator it = p.begin(); it != p.end() ; ++it)
      back_inserter(g)  =  monomial(P, 1, *it);
    return ideal(P, g);
  }//interreduce


// this can be done with remove_copy_if
  vector<PPMonoidElem> interreduce(const vector<PPMonoidElem>& p)
  {
    vector<PPMonoidElem> res;
    for (vector<PPMonoidElem>::const_iterator it = p.begin();it != p.end();++it)
      if (!IsStrictlyDivisible(*it, p))
        res.push_back(*it);
    sort(res.begin(), res.end());
    vector<PPMonoidElem>::iterator ToErase = unique(res.begin(),res.end());
    res.erase(ToErase,res.end());
    return res;
  }//interreduce


  PPMonoidElem SQFR(ConstRefPPMonoidElem T)
  {
    vector<long> expv;
    exponents(expv,T);
    replace_if(expv.begin(),expv.end(),bind2nd(not_equal_to<long>(),0),1);
    return PPMonoidElem(owner(T),expv);
  }//SQFR


  vector<PPMonoidElem> colon(const vector<PPMonoidElem>& L,ConstRefPPMonoidElem T)
  {
    vector<PPMonoidElem> result;
    for (vector<PPMonoidElem>::const_iterator i = L.begin(); i != L.end(); ++i)
      result.push_back(colon(*i,T));
    return result;
  }


/////////////////////////// REAL CODE

  ideal PPMonoidElemVector2Ideal(const SparsePolyRing& the_SPR,
                                  const vector<PPMonoidElem>& L)
  {  
    vector<RingElem> DummyGens;
    for (vector<PPMonoidElem>::const_iterator it = L.begin();
         it != L.end();++it)
      DummyGens.push_back(monomial(the_SPR,1,*it));
    return ideal(the_SPR,DummyGens);
  }//PPMonoidElemVector2Ideal


  long TmpDim(const vector<PPMonoidElem>& the_V)
  {  
    ring Z = RingZZ();
    ring ZMod2 = NewZZmod(2);
    SparsePolyRing PPRing = NewPolyRing(ZMod2,owner(the_V.front()));
    ideal J(PPMonoidElemVector2Ideal(PPRing,interreduce(the_V)));
    return DimQuot(J);
  }//dim


// Compares std lex two univariate polys
  bool IsLexLower(ConstRefRingElem the_f,ConstRefRingElem the_g)
  {
    SparsePolyIter it_f = BeginIter(the_f);
    SparsePolyIter it_g = BeginIter(the_g);
    for (; !IsEnded(it_f)&&!IsEnded(it_g); ++it_f,++it_g)
    {
      unsigned long df = StdDeg(PP(it_f));
      unsigned long dg = StdDeg(PP(it_g));
      if (df<dg) return true;
      if (df>dg) return false;
      if (coeff(it_f)<coeff(it_g)) return true;
      if (coeff(it_f)>coeff(it_g)) return false;
    }
    return false; 
  }//IsLexLower


  bool IsHPLexLower(const vector<PPMonoidElem>& the_V,// may be non const to save
                    ConstRefPPMonoidElem the_t1,
                    ConstRefPPMonoidElem the_t2,
                    const SparsePolyRing& the_PPRing)
  {
    CoCoA_ASSERT_ALWAYS(PPM(the_PPRing) == owner(the_t1));
    CoCoA_ASSERT_ALWAYS(owner(the_t2) == owner(the_t2));
    vector<PPMonoidElem> V = the_V;
  
    // Computing the H-poly of the_V + the_t1
    V.push_back(the_t1);
    ideal J(PPMonoidElemVector2Ideal(the_PPRing,interreduce(V)));
    RingElem HNum1(HilbertNumQuot(J));
    V.pop_back();
  
    // Computing the H-poly of the_V + the_t2
    V.push_back(the_t2);
    J = ideal(PPMonoidElemVector2Ideal(the_PPRing,interreduce(V)));
    RingElem HNum2(HilbertNumQuot(J));
    V.pop_back();
  
    clog<<"IsLexLower: HNum1 = "<<HNum1<<endl;
    clog<<"IsLexLower: HNum2 = "<<HNum2<<endl;
    clog<<"IsLexLower: HNum1<HNum2 = "<<IsLexLower(HNum1,HNum2)<<endl;

    return IsLexLower(HNum1,HNum2);
  }//IsHPLexLower


  bool IsHPEqual(const vector<PPMonoidElem>& the_V,// may be non const to save
                 ConstRefPPMonoidElem the_t1,        
                 ConstRefPPMonoidElem the_t2,
                 const SparsePolyRing& the_PPRing)
  {
    CoCoA_ASSERT_ALWAYS(PPM(the_PPRing) == owner(the_t1));
    CoCoA_ASSERT_ALWAYS(owner(the_t2) == owner(the_t2));
    vector<PPMonoidElem> V = the_V;
  
    // Computing the H-poly of the_V + the_t1
    V.push_back(the_t1);
    ideal J(PPMonoidElemVector2Ideal(the_PPRing,interreduce(V)));
    RingElem HNum1(HilbertNumQuot(J));
    V.pop_back();
  
    // Computing the H-poly of the_V + the_t2
    V.push_back(the_t2);
    J = ideal(PPMonoidElemVector2Ideal(the_PPRing,interreduce(V)));
    RingElem HNum2(HilbertNumQuot(J));
    V.pop_back();
  
    clog<<"IsHPEqual: HNum1 = "<<HNum1<<endl;
    clog<<"IsHPEqual: HNum2 = "<<HNum2<<endl;
    clog<<"IsHPEqual: HNum1 == HNum2 "<<(HNum1 == HNum2)<<endl;
  
    return HNum1 == HNum2;
  }//IsHPEqual


/*
  the_OldPP may be empty, the_CurrentPP can't. 
  The elements of the_CurrentPP can't be one (constants can't be leading terms)
  Kills all the elements in the_CurrentPP which are Hilbert lex less than the others
  optimizations: use more complex data types, filling in dim and HPpoly and then
  reasoning with them
  use lists to optimize the killing.
*/
  void HilbertPrune(const vector<PPMonoidElem>& the_OldPP,
                    vector<PPMonoidElem>& the_CurrentPP)
  {
    // pruning by degree
    clog<<"HilbertPrune: Candidates = "<<the_CurrentPP<<endl;
  
    CoCoA_ASSERT_ALWAYS(!the_CurrentPP.empty());  
    // computing minimum degree
    long CurrentDeg = 0;
    long MinDeg = StdDeg(the_CurrentPP.front());
    vector<PPMonoidElem>::const_iterator BeginPlus = the_CurrentPP.begin();
    ++BeginPlus;
    for (vector<PPMonoidElem>::const_iterator i = BeginPlus;
         i != the_CurrentPP.end(); ++i)
    {  
      CurrentDeg = StdDeg(*i);
      if (CurrentDeg<MinDeg)
        MinDeg = CurrentDeg;
    }

    // pruning by degree
    vector<PPMonoidElem> tmp;
    for (vector<PPMonoidElem>::const_iterator it = the_CurrentPP.begin();
         it != the_CurrentPP.end();
         ++it)
      if (StdDeg(*it) == MinDeg)
        tmp.push_back(*it);
    swap(tmp,the_CurrentPP);
    if (the_CurrentPP.size() == 1)
      return;
  
//clog<<"HilbertPrune: after deg pruning = "<<the_CurrentPP<<endl;

    // pruning by dim
    //  ring Z  =  RingZZ();
    ring QQ  =  RingQQ();
    //    SparsePolyRing HPRing  =  NewPolyRing(QQ, 1);// The ring for the HPoly
    SparsePolyRing PPRing  =  NewPolyRing(NewZZmod(2),owner(the_CurrentPP.front()));//The utility poly ring associated to the vector<PPMonoidElem> 
    //    RingElem BinPower = power(1-indet(HPRing,0),NumIndets(PPRing));

  
    // computing MinDim
    vector<PPMonoidElem> PPVListSQFR;
    transform(the_OldPP.begin(),the_OldPP.end(),back_inserter(PPVListSQFR),SQFR);
    PPVListSQFR.push_back(SQFR(the_CurrentPP.front()));
    long MinDim = TmpDim(PPVListSQFR);
    PPVListSQFR.pop_back();
    clog<<"HilbertPrune: StartingDim = "<<MinDim<<" "<<SQFR(the_CurrentPP.front())<<endl;

    long CurrentDim = 0;
    BeginPlus = the_CurrentPP.begin();++BeginPlus;
    for (vector<PPMonoidElem>::const_iterator it = BeginPlus;
         it != the_CurrentPP.end(); 
         ++it)
    { 
      PPVListSQFR.push_back(SQFR(*it));
      CurrentDim = TmpDim(PPVListSQFR);
      PPVListSQFR.pop_back();
      clog<<"HilbertPrune:  CurrentDim = "<<CurrentDim<<" "<<SQFR(*it)<<endl;
      if (CurrentDim<MinDim)
        MinDim = CurrentDim;
    }
   

// dim pruning
    tmp.clear();
    for (vector<PPMonoidElem>::const_iterator it = the_CurrentPP.begin();
         it != the_CurrentPP.end();
         ++it)
    {  
      PPVListSQFR.push_back(SQFR(*it));
      if (TmpDim(PPVListSQFR) == MinDim)
        tmp.push_back(*it);
      PPVListSQFR.pop_back();
    }
    swap(tmp,the_CurrentPP);
    if (the_CurrentPP.size() == 1)
      return;

    clog<<"HilbertPrune: after dim pruning = "<<the_CurrentPP<<endl;
    
// HP pruning

    // determining the position of the lex lower PP
    vector<PPMonoidElem>::const_iterator LexMin = the_CurrentPP.begin();
    BeginPlus = the_CurrentPP.begin();BeginPlus++;
    for (vector<PPMonoidElem>::const_iterator it = BeginPlus;it != the_CurrentPP.end();++it)
      if (IsHPLexLower(the_OldPP,*it,*LexMin,PPRing))
        LexMin = it;
  
    // pruning by lex lower
    tmp.clear();
    for (vector<PPMonoidElem>::iterator it = the_CurrentPP.begin();
         it != the_CurrentPP.end();
         ++it)
      if (IsHPEqual(the_OldPP,*it,*LexMin,PPRing))
        tmp.push_back(*it);
    swap(tmp,the_CurrentPP);
  
  
//clog<<"HilbertPrune: after HP pruning = "<<the_CurrentPP<<endl;
  }//HilbertPrune

// given the PP t1, t2, returns the vector representi the inequality for the
// var weights necessary to have t1>t2 [ie -t1+t2+1 <= 0]
  vector<long> CreateInquality(ConstRefPPMonoidElem t1,ConstRefPPMonoidElem t2)
  {
    vector<long> expv1, expv2;
    exponents(expv1,t1);
    exponents(expv2,t2);
    vector<long> v(NumIndets(owner(t1))+1);
    for (unsigned int i = 0;i != expv1.size();++i)
      v[i] = expv2[i]-expv1[i];
    v[v.size()-1] = 1;// to have the inequality strict
    return v;
  }//CreateInquality


// Q has to be the rationals
// t is one element of f = <t1,..,tn>, which hence is not empty
// result is the matrix whose rows give the inequalities t>t1,..,t>tn
// To guarantee M is without useless inequalities, kill the elements in f that
// divide some other element of f
// optimization: don't add the rows whose real var vector is the multiple of some othe
// real var vector
  matrix CreateSimplexMatrix(const ring& QQ,
                             ConstRefPPMonoidElem t,
                             const vector<PPMonoidElem>& f)
  {
    vector<vector<long> > result;
    for (vector<PPMonoidElem>::const_iterator it = f.begin();it != f.end();++it)
      if (*it != t)
        result.push_back(CreateInquality(t,*it));

    matrix M = NewDenseMat(QQ,f.size(),NumIndets(owner(t))+2);
 
    // setting the inequalities
    for (unsigned int i = 0;i != result.size();++i)
      for (unsigned int j = 0;j != result[i].size();++j)
        SetEntry(M,i,j,result[i][j]);
      

    // setting the marginal labels row, one for each vars
    for (unsigned int j = 0;j != NumIndets(owner(t));++j)
      SetEntry(M,NumRows(M)-1,j,j);
   
    // setting the marginal labels col, one for each row execpt the last
    for (unsigned int i = 0;i != NumRows(M)-1;++i)
      SetEntry(M,i,NumCols(M)-1,NumIndets(owner(t))+i);

    return M;
  }//CreateSimplexMatrix


  void pivot(matrix& M,const unsigned int i, const unsigned int j)
  {
    RingElem p = M(i,j);
    unsigned int MarginRow = NumRows(M)-1;
    unsigned int MarginCol = NumCols(M)-1;
  
    for (unsigned int l = 0;l<i;++l)
      for (unsigned int k = 0;k<j;++k)
        SetEntry(M,l,k,M(l,k)-M(i,k)*M(l,j)/p);

    for (unsigned int l = i+1;l<MarginRow;++l)
      for (unsigned int k = 0;k<j;++k)
        SetEntry(M,l,k,M(l,k)-M(i,k)*M(l,j)/p);

    for (unsigned int l = 0;l<i;++l)
      for (unsigned int k = j+1;k<MarginCol;++k)
        SetEntry(M,l,k,M(l,k)-M(i,k)*M(l,j)/p);
 
    for (unsigned int l = i+1;l<MarginRow;++l)
      for (unsigned int k = j+1;k<MarginCol;++k)
        SetEntry(M,l,k,M(l,k)-M(i,k)*M(l,j)/p);
   
    // doing the pivot row and col
    M->myRowMul(i,1/p);
    M->myColMul(j,-1/p);
    SetEntry(M,i,j,-M(i,j));// resetting the pivot
    SetEntry(M,i,MarginCol,M(i,MarginCol)*p);// resetting the margin
    SetEntry(M,MarginRow,j,-M(MarginRow,j)*p);// resetting the margin

    // swapping the marginal labels
  
//clog << "\n swapping i,j  = " <<i<<","<<j<< "\n"<<endl;

    const RingElem dummy = M(i,MarginCol);
    SetEntry(M,i,MarginCol,M(MarginRow,j));
//clog << "\n M  = " <<M<<endl;
    SetEntry(M,MarginRow,j,dummy);
//clog << "\n M  = " <<M<<endl;

    return;
  }//pivot


  enum PivotPossibleResult {SolutionExist,SolutionDoesNotExist,StillWorking};


  ostream& operator<<(ostream& out, PivotPossibleResult res)
  {
    switch (res)
    {
    case SolutionExist:        out<<"SolutionExist";        break;
    case SolutionDoesNotExist: out<<"SolutionDoesNotExist"; break;
    case StillWorking:         out<<"StillWorking";         break;
    default:            out << "-- unknown value --"; break;
    }
    return out;
  }//operator<<


// Look for pivot
  PivotPossibleResult LookPivot(const matrix& the_M,unsigned int& the_i,unsigned int& the_j)
  {
//clog << "LookPivot: begin"<<endl;
//clog << "LookPivot: the_M"<<the_M<<endl;
    ring P = RingOf(the_M);
    unsigned int LastRealRow = NumRows(the_M)-2;
    unsigned int LastRealCol = NumCols(the_M)-3;
    unsigned int CoeffCol = NumCols(the_M)-2;
    // looking for the positive coeff
    unsigned int l = LastRealRow;
    while (l != 0&&the_M(l,CoeffCol) <= 0)
      --l;
    if (l == 0&&the_M(l,CoeffCol) <= 0) 
      return SolutionExist;
  
    //looking for pivot col
    unsigned int j = 0;
    while (j<LastRealCol&&the_M(l,j) >= 0)
      ++j;

    if (j == LastRealCol&&the_M(l,j) >= 0) 
      return SolutionDoesNotExist;
    the_j = j;//pivot col found

    //looking for pivot row
    unsigned int CandidatePivotRow = l;
    for (unsigned int i = l;i <= LastRealRow;++i)
    {
      if (the_M(i,the_j)>zero(P)
          &&
          the_M(l,CoeffCol)/the_M(i,the_j)<the_M(l,CoeffCol)/the_M(CandidatePivotRow,the_j))
        CandidatePivotRow = i;
    }
    the_i = CandidatePivotRow;//pivot row found

    return StillWorking;
  }//LookForPivot


// v != emptyset
  RingElem lcm(const vector<RingElem>& v)
  {
 
    CoCoA_ASSERT_ALWAYS(!v.empty());
    CoCoA_ASSERT_ALWAYS(IsTrueGCDDomain(owner(v.front())));
    RingElem d = v.front();
    vector<RingElem>::const_iterator it = v.begin();++it;
    for (;it != v.end();++it)
      d = (*it)*d/gcd(d,*it);
    return d;
  }//gcd

// multiplies by the gcd of the denoms
  vector<RingElem> normalize(const vector<RingElem>& v)
  {
    if (v.empty())
      return v;

    const FractionField& P = owner(v.front());
    CoCoA_ASSERT_ALWAYS(IsFractionFieldOfGCDDomain(P));
  
    vector<RingElem> tmp;
    for (vector<RingElem>::const_iterator it = v.begin();it != v.end();++it)
      if (!IsZero(*it))
        tmp.push_back(den(*it));
    RingElem d(BaseRing(P));// should be ZZ
    RingHom phi = EmbeddingHom(P);
    vector<RingElem> res;
    if (tmp.size() == 0) 
      d = 1;
    else
      d = lcm(tmp);
    for (vector<RingElem>::const_iterator it = v.begin();it != v.end();++it)
      res.push_back(num((*it)*phi(d)));
    return res;
  }//normalize


  vector<RingElem> DividesByGCD(const vector<RingElem>& the_v)
  {
    if (the_v.empty())
      return the_v;
    ring P = owner(the_v.front());
    CoCoA_ASSERT_ALWAYS(IsTrueGCDDomain(P));
    RingElem d = the_v.front();
    vector<RingElem>::const_iterator it = the_v.begin();++it;
    for (;it != the_v.end();++it)
      d = gcd(d,*it);
    vector<RingElem> v;
    for (vector<RingElem>::const_iterator it = the_v.begin();it != the_v.end();++it)
      v.push_back((*it)*d);
    return v;
  }//DividesByGCD


// Short term, while I understand better the canonical homomorphisms
// crashes if x is not an unsigned int, and this is right
  unsigned int Convert2UInt(RingElem x)
  {
    BigInt N;
    unsigned int n = 0;
    bool b =  IsInteger(N,x);
    if (!b)
      CoCoA_THROW_ERROR("convert2UInt","ring element was NOT an unsigned int");
    IsConvertible(n,N);
    return n;
  }//convert2UInt


// the_v has to be over QQ, but really over ZZ. This converts it to an unsigned int
// vector
  vector<unsigned int> ConvertV2UIntV(const vector<RingElem> the_v)
  {
    vector<unsigned int> v;
    for (vector<RingElem>::const_iterator it = the_v.begin();it != the_v.end();++it)
      v.push_back(Convert2UInt(*it));
    return v;
  }//QQv2UIntv


  vector<RingElem> WVector(const matrix& M)
  { 
//clog << "WVector: begin" <<endl;
    unsigned int NumIndets = NumCols(M)-2;
    vector<RingElem> v(NumIndets,zero(RingOf(M)));
    unsigned int CoeffCol = NumCols(M)-2;
    unsigned int MarginCol = NumCols(M)-1;
    unsigned int MarginRow = NumRows(M)-1;
    BigInt N;
    unsigned int num;

    for (unsigned int i = 0;i != MarginRow;++i)
      if (M(i,MarginCol)<NumIndets)
      {
        IsInteger(N,M(i,MarginCol));// if not, horrible error!
        IsConvertible(num,N);// if not possible, horrible error!
        v[num] = -M(i,CoeffCol);
      }
    return normalize(v);
  }//WVector


//Inequalities system solving
// result is a solution of the inequalities system if not empty.
// if the result is empty, no solution exists.
// MODIFIES the matrix M!!!!!
  vector<RingElem> ISSolving(matrix& M)
  {  
    PivotPossibleResult status;
    unsigned int pivot_i,pivot_j;
    do 
    {
      status = LookPivot(M,pivot_i,pivot_j);

      if (status == StillWorking)
        pivot(M,pivot_i,pivot_j);
    }
    while (status == StillWorking);
    vector<RingElem> sol;
    if (status == SolutionDoesNotExist) return sol;
    return WVector(M);
  }//ISSolving


// Short term, while I understand better the canonical homomorphisms
// crashes if x is not an unsigned int, and this is right
  unsigned int convert2UInt(RingElem x)
  {
    BigInt N;
    unsigned int n = 0;
    bool b = IsInteger(N,x);
    if (!b)
      CoCoA_THROW_ERROR("convert2UInt","ring element was NOT an unsigned int");
    IsConvertible(n,N);
    return n;
  }//convert2UInt


// check if the tableau is Standard ie no real var has been exchanged for a
// slack var
  bool IsTableauStandard(const matrix& M)
  {
    for (unsigned int i = 0;i != NumCols(M)-2;++i)
      if (Convert2UInt(M(NumRows(M)-1,i)) != i)
        return false;
    return true;
  }//CheckTableauCompatibility


// the_a = the_b in top of the_c. In the context of the special matrix
// representing symplextableaus only!!!! 
// DOUBLE WARN works only if both tableaus are in standard form 
  matrix combine(const matrix& the_b,const matrix& the_c)
  {
    CoCoA_ASSERT_ALWAYS(NumCols(the_b) == NumCols(the_c));
    //CoCoA_ASSERT_ALWAYS(IsTableauStandard(the_b)||IsTableauStandard(the_c));
    // checking the tableaus compatibility
    if (!IsTableauStandard(the_b)||!IsTableauStandard(the_c))
      CoCoA_THROW_ERROR("combine","Dynamic GB:combine wrong tableau matrices");
    
    matrix a(NewDenseMat(RingOf(the_b),NumRows(the_b)+NumRows(the_c)-1,NumCols(the_c)));
    for (unsigned int i = 0;i != NumRows(the_b);++i)
      for (unsigned int j = 0;j != NumCols(the_b);++j)
        SetEntry(a,i,j,the_b(i,j));
    for (unsigned int i = 0;i != NumRows(the_c);++i)
      for (unsigned int j = 0;j != NumCols(the_c)-1;++j)
        SetEntry(a,i+NumRows(the_b)-1,j,the_b(i,j));
    
    RingElem LastSlackVarOfTheB(the_b(NumRows(the_b)-2,NumCols(the_b)-1));
    for (unsigned int i = 0;i != NumRows(the_c)-1;++i)
      SetEntry(a,i+NumRows(the_b)-1,NumCols(the_c)-1,LastSlackVarOfTheB+i+1);
    return a;   
  }//combine


// produces a matrix with two more cols for const and margin, and one more row for margin
  matrix CreateNewEmptyTableau(const ring& QQ,const unsigned int n,const unsigned int m)
  {
    return NewDenseMat(QQ,n+1,m+2);
  }


// V is NOT empty. V is already a matrix. Check if and when.
  matrix VectorVector2Matrix(const ring& P,const vector<vector<RingElem> >& V)
  { 
    unsigned int n = V.size();
    unsigned int m = V.front().size();
    matrix M = NewDenseMat(P,n,m);
    for (unsigned int i = 0;i != n;++i)
      for (unsigned int j = 0;j != m;++j)
        SetEntry(M,i,j,V[i][j]);
    return M;
  }//


  vector<vector<RingElem> > Matrix2VectorVector(const matrix& M)
  { 
    vector<vector<RingElem> > V;
    vector<RingElem> v;
    for (unsigned int i = 0;i != NumRows(M);++i)
    { 
      for (unsigned int j = 0;j != NumCols(M);++j)
        v.push_back(M(i,j));
      V.push_back(v);
      v.clear();
    }
    return V;
  }//


  void ImpossibleLPPPrune(const matrix& the_M,
                          vector<PPMonoidElem>& the_candidates,
                          vector<vector<RingElem> >& the_solutions)
  {
//cout << "ImpossibleLPPPrune: begin "<<endl;
    CoCoA_ASSERT_ALWAYS(IsQQ(RingOf(the_M)));
    vector<vector<RingElem> > solutions;
    if (the_candidates.empty())
    {
      the_solutions.clear();
      return;
    }
    const FractionField& QQ = RingOf(the_M);
    vector<PPMonoidElem> pruned;
    PPMonoid PPM1 = owner(the_candidates.front());
    vector<RingElem> sol;
//cout << "ImpossibleLPPPrune: before tableau creation "<<endl;
  
    // used for the newly created mat
    matrix tmp1 = CreateNewEmptyTableau(QQ,the_candidates.size()-1,NumIndets(PPM1));
    // this has space enough for the tableau union of the_M and tmp1
    matrix tmp2 = CreateNewEmptyTableau(QQ,the_candidates.size()+NumRows(the_M)-2,NumIndets(PPM1));// used for the newly created mat

    for (vector<PPMonoidElem>::const_iterator it = the_candidates.begin();
         it != the_candidates.end(); ++it)
    { 
      AssignMatrix(tmp1,CreateSimplexMatrix(QQ,*it,the_candidates));
      sol = ISSolving(tmp1);
      if (!sol.empty())
      {
        pruned.push_back(*it);
        solutions.push_back(sol);
      }
    }
    swap(the_candidates,pruned);
    swap(the_solutions,solutions);
  }//ImpossibleLPPPrune


// write a class for pp+sol+matrix and use that
  void SelectOptimalSolutions(const vector<PPMonoidElem>& the_candidates,
                              const vector<PPMonoidElem>& the_PossibleCandidates,
                              vector<vector<RingElem> >& the_PossibleSolutions,
                              vector<vector<RingElem> >& the_OptimalSolutions)
  {
    vector<vector<RingElem> > OptimalSolutions;
    for (vector<PPMonoidElem>::const_iterator it = the_candidates.begin();
         it != the_candidates.end();
         ++it)
    {
      unsigned int i = 0;
      while (i != the_PossibleCandidates.size()&&the_PossibleCandidates[i] != *it)
        ++i;
      if (i != the_PossibleCandidates.size())
        OptimalSolutions.push_back(the_PossibleSolutions[i]);
    }
    swap(the_OptimalSolutions,OptimalSolutions);			  
  }//SelectOptimalSolutions


  void CreateLPPCandidates(const vector<PPMonoidElem> the_OldPP,
                           const matrix& the_OldPPMat,
                           ConstRefRingElem the_f,
                           vector<PPMonoidElem>& the_OptimalCandidates,
                           vector<PPMonoidElem>& the_PossibleCandidates,
                           vector<vector<RingElem> >& the_OptimalSolutions,
                           vector<vector<RingElem> >& the_PossibleSolutions)
  {
    CoCoA_ASSERT_ALWAYS(IsSparsePolyRing(owner(the_f)));
    vector<PPMonoidElem> candidates = interreduce(supp(the_f));
    vector<vector<RingElem> > PossibleSolutions;
    vector<vector<RingElem> > OptimalSolutions;
    ImpossibleLPPPrune(the_OldPPMat,candidates,PossibleSolutions);
    vector<PPMonoidElem> PossibleCandidates = candidates;
    HilbertPrune(the_OldPP,candidates);
    SelectOptimalSolutions(candidates,PossibleCandidates,PossibleSolutions,OptimalSolutions);
    swap(candidates,the_OptimalCandidates);
    swap(PossibleCandidates,the_PossibleCandidates);
    swap(PossibleSolutions,the_PossibleSolutions);
    swap(OptimalSolutions,the_OptimalSolutions);
  }//PossibleLPP


  void ChooseBestLPPAndUpdate(vector<PPMonoidElem> the_OldPP,
                              vector<vector<RingElem> >& the_OldPPRelations,
                              ConstRefRingElem the_f,
                              RingElem the_best,
                              vector<RingElem>& the_solution)
  {
    vector<PPMonoidElem> OptimalCandidates;
    vector<PPMonoidElem> PossibleCandidates;
    vector<vector<RingElem> > solutions;
    vector<vector<RingElem> > OptimalSolutions;
    FractionField QQ  =  RingQQ();
    CreateLPPCandidates(the_OldPP,VectorVector2Matrix(QQ,the_OldPPRelations),the_f,
                        OptimalCandidates,PossibleCandidates,solutions,OptimalSolutions);// those are filled up
    const SparsePolyRing& P  =  owner(the_best);
    RingElem best = monomial(P,1,OptimalCandidates.front());
    vector<RingElem> solution = solutions.front();

// Updating therelations
    matrix NewConditions(CreateSimplexMatrix(QQ,LPP(best),PossibleCandidates));
    matrix tmp = VectorVector2Matrix(QQ,the_OldPPRelations);
    the_OldPPRelations = Matrix2VectorVector(combine(tmp,NewConditions));

    the_OldPP.push_back(LPP(best));
    swap(best,the_best);
    swap(solution,the_solution);
  }//ChooseBestLPPAndUpdate


  void ChooseBestLPP(const vector<PPMonoidElem> the_OldPP,
                     ConstRefRingElem the_f,
                     RingElem the_best,
                     vector<RingElem>& the_solution)
  {
    vector<PPMonoidElem> OptimalCandidates;
    vector<PPMonoidElem> PossibleCandidates;
    vector<vector<RingElem> > OldPPRelations_dummy;
    vector<vector<RingElem> > solutions;
    vector<vector<RingElem> > OptimalSolutions;
    FractionField QQ  =  RingQQ();
    CreateLPPCandidates(the_OldPP,VectorVector2Matrix(QQ,OldPPRelations_dummy),the_f,
                        OptimalCandidates,PossibleCandidates,OptimalSolutions,solutions);// those are filled up
    const SparsePolyRing& P  =  owner(the_best);
    RingElem best = monomial(P,1,OptimalCandidates.front());
    vector<RingElem> solution = solutions.front();
    cout << "ChooseBestLPP: PossibleCandidates  "<<PossibleCandidates<<endl;
    cout << "ChooseBestLPP: OptimalCandidates  "<<OptimalCandidates<<endl;
    cout << "ChooseBestLPP: solutions  "<<solutions<<endl;
    cout << "ChooseBestLPP: best  "<<best<<endl;
    cout << "ChooseBestLPP: solution  "<<solution<<endl;

    swap(best,the_best);
    swap(solution,the_solution);
  }//ChooseBestLPPAndUpdate


  SparsePolyRing NewPolyRingA(const SparsePolyRing& P,
                                    const vector<RingElem>& the_OrdVector)
  {
    CoCoA_ASSERT_ALWAYS(the_OrdVector.size() != NumIndets(P));
    ring EntryRing = owner(the_OrdVector.front());
    // the_OrdVector+1-vector
    matrix M = NewDenseMat(EntryRing,2,the_OrdVector.size());
    for (unsigned int i = 0;i != NumCols(M);++i)
    {
      SetEntry(M,0,i,the_OrdVector[i]);
      SetEntry(M,1,i,1);
    }
    matrix OrdM = MakeTermOrdMat(M);
    clog << "NewPolyRingA:OrdMatrix = " <<OrdM<<endl;
    const PPOrdering ord  =  NewMatrixOrdering(OrdM, 0);
    vector<symbol> IndetNames = symbols(PPM(P));
    return NewPolyRing(CoeffRing(P),IndetNames,ord);
  }//NewPolyRing


  RingHom OrderingChangeRingHom(const SparsePolyRing& the_OldP,
                                const SparsePolyRing& the_NewP)
  {
    CoCoA_ASSERT_ALWAYS(NumIndets(the_OldP) == NumIndets(the_NewP));
  
    vector<RingElem> images;
    for (unsigned int i = 0;i != NumIndets(the_OldP);++i)
      images.push_back(indet(the_NewP,i));    
    return PolyRingHom(the_OldP,
                       the_NewP,
                       CoeffEmbeddingHom(the_NewP),
                       images);
  }//OrderingChangeRingHom


  vector<RingElem> ChangeOrdering(const SparsePolyRing& the_NewSPR,
                                  const vector<RingElem>& the_V)
  {
    vector<RingElem> V;
    if (the_V.empty())
      return V;
    SparsePolyRing OldSPR = owner(the_V.front());
    RingHom phi = OrderingChangeRingHom(OldSPR,the_NewSPR);
    for (vector<RingElem>::const_iterator it = the_V.begin();it != the_V.end();++it)
      V.push_back(phi(*it));
    return V;
  }//ChangeOrdering

// DYNAMIC ALG: reduces until a non zero SPoly has a LPP with is different from the best possibile LPP w.r.t. HPoly.
// In that case, settles the WrongLPP filed in the reductor. For each LPP found, updates the LPP list.
// when this procedure stops, if myWrongLPPFoundValue is true, we have a GB
// otherwise, we have to resume computations from the candidate basis plus the
// special pairs                           
  void GReductor::myReduceUntilWrongLPPFound(RefPPMonoidElem theBestLPP,
                                             vector<RingElem>& theBestSol) 
  {
    vector<PPMonoidElem> CandidateBasisLPP; 
    ring QQ  =  RingQQ();

    matrix OldConstraints(NewDenseMat(QQ,0,0));// empty at the moment
    myPrepareGBasis(); 
    myStampaPairs(clog);   
    while (!myPairs.empty())
    {
      myReduceCurrentSPoly();
         
      if (!IsZero(mySPoly))
      {
        vector<PPMonoidElem> PossibleCandidates;
        vector<PPMonoidElem> OptimalCandidates;
        vector<vector<RingElem> > solutions;
        vector<vector<RingElem> > OptimalSolutions;
        CreateLPPCandidates(CandidateBasisLPP,OldConstraints,poly(mySPoly),
                            OptimalCandidates,PossibleCandidates,
                            OptimalSolutions,solutions);
        clog << "myReduceUntilWrongLPPFound:CandidateBasisLPP " <<CandidateBasisLPP<<endl;
        clog << "myReduceUntilWrongLPPFound:mySPoly " <<supp(poly(mySPoly))<<endl;
        clog << "myReduceUntilWrongLPPFound:Actual LPP " <<LPPForOrd(mySPoly)<<endl;
        clog << "myReduceUntilWrongLPPFound:OptimalCandidates " <<OptimalCandidates<<endl;
        clog << "myReduceUntilWrongLPPFound:Actual among Optimal? "
             <<(OptimalCandidates.end()
                != 
                find(OptimalCandidates.begin(),OptimalCandidates.end(),LPPForOrd(mySPoly)))
             <<endl;
        if (OptimalCandidates.end()
            == 
            find(OptimalCandidates.begin(),OptimalCandidates.end(),LPPForOrd(mySPoly)))
        {       
          clog << "myReduceUntilWrongLPPFound:WRONG LPP FOUND "<<endl;
          myWrongLPPFoundValue = true;// we have the wrong ordering, resume
          // computations from the candidate basis plus the special pairs
          swap(theBestLPP,OptimalCandidates.front());
          swap(theBestSol,OptimalSolutions.front());
          myUpdateBasisAndPairs();
          return;
        }
        CandidateBasisLPP.push_back(LPPForOrd(mySPoly));
        myUpdateBasisAndPairs();
      }
    }//while
    myFinalizeGBasis();
  }//myReduceUntilWrongLPPFound


  // to GOperations
  void ComputeDynamicGBasis(PolyList& theGB,const PolyList& thePL, int level = 0)
  {
    clog << "ComputeDynamicGBasis:level =  " <<level<<endl;
    PolyList GB;
    if (thePL.empty())
    {
      swap(theGB,GB);
      return;
    }
    PolyList NewCandidateBasis;
    {
      const SparsePolyRing& CurrentSPR(owner(thePL));
      PPMonoidElem BestLPP(PPM(CurrentSPR));
      vector<RingElem> BestSolution;
      bool ArePolyHomogeneous = IsHomog(thePL);
      bool IsSatAlg = false;
      GRingInfo GRI(CurrentSPR,ArePolyHomogeneous,IsSatAlg,NewDivMaskEvenPowers());
      clog <<"ComputeDynamicGBasis:ArePolyHomogeneous = "<<ArePolyHomogeneous<<endl;
      GReductor GBR(GRI, thePL);
      GBR.myReduceUntilWrongLPPFound(BestLPP,BestSolution);// dynamic algorithm
      if (!GBR.WrongLPPFound())
      {
        GBR.myGBasis(GB);
        swap(theGB,GB);
        return;
      }
      clog << "ComputeDynamicGBasis:WrongLPPFound " <<endl;
      SparsePolyRing TmpRing = NewPolyRingA(CurrentSPR,BestSolution);
      RingHom phi = OrderingChangeRingHom(CurrentSPR,TmpRing);
      PolyList CandidateBasis;
      clog << "ComputeDynamicGBasis:Pairs "<<endl;
      GBR.myStampaPairs(clog);
      GBR.GetCandidateGBasis(CandidateBasis);
      clog << "ComputeDynamicGBasis:CandidateBasis " <<CandidateBasis<<endl;
      NewCandidateBasis = ChangeOrdering(TmpRing,CandidateBasis);
      clog << "ComputeDynamicGBasis:BestLPP " <<BestLPP<<endl;
      clog << "ComputeDynamicGBasis:BestSolution " <<BestSolution<<endl;
      clog << "ComputeDynamicGBasis:TmpRing " <<TmpRing<<endl;
      clog << "ComputeDynamicGBasis:phi " <<phi<<endl;
      clog << "ComputeDynamicGBasis:NewCandidateBasis " <<NewCandidateBasis<<endl;
      GB.clear();
    }
    ComputeDynamicGBasis(GB,NewCandidateBasis,level+1);
    swap(theGB,GB);
    clog << "ComputeDynamicGBasis: theGB =  " <<theGB<<endl;
  }//ComputeGBasis
		  

  void program()
  {
    // This is a test for Dynamic Buchberger's Algorithm:
  
/*  OLD TESTING  
    ring Z  =  RingZZ();
    ring ZMod  =  NewZZmod(Z, 101);
    ring Q  =  NewFractionField(Z);
    vector<symbol> X  =  range(symbol("x",0), symbol("x",4));
    SparsePolyRing P(NewPolyRing_DMPI(QQ, X));
    RingElem x  =  indet(P,0),  y  =  indet(P,1);
    RingElem z  =  indet(P,2), t = indet(P,3);
    RingElem h  =  indet(P,4);
    vector<PPMonoidElem> OldLPP;
    OldLPP.push_back(LPP(x));
    matrix OldM(NewDenseMatrix(QQ,0,0));
    RingElem f = x*y+power(x,2)+power(y,2)+power(z,2)+power(t,2)+t*z;
    clog << "testing PossibleLPP poly =  " << f << endl;
    vector<PPMonoidElem> PossibleCandidates;
    vector<PPMonoidElem> OptimalCandidates;
    vector<vector<RingElem > > solutions;
    vector<vector<RingElem > > OptimalSolutions;
    CreateLPPCandidates(OldLPP,OldM,f,OptimalCandidates,PossibleCandidates,OptimalSolutions,solutions);
    clog << "PossibleCandidates+ =  " << PossibleCandidates << endl;
    clog << "OptimalCandidates =  " << OptimalCandidates << endl;
    clog << "solutions =  " << solutions << endl;
    clog << "OptimalSolutions =  " << OptimalSolutions << endl;

    vector<RingElem> OrdVector;
    OrdVector.push_back(one(Z));
    OrdVector.push_back(2*one(Z));
    OrdVector.push_back(3*one(Z));
    OrdVector.push_back(zero(Z));
    OrdVector.push_back(20*one(Z));
    SparsePolyRing NewP =  NewPolyRing(P,OrdVector);

    clog << "OrdVector =  " << OrdVector << endl;
    clog << "NewP =  " << NewP << endl;
    clog << "P =  " << P << endl;

    RingHom phi  =  OrderingChangeRingHom(P,NewP);
    RingElem g = x*y+power(x,2)+power(y,2)+power(z,2)+power(t,2)+t*z+h;
    clog << "g =  " << g << endl;
    clog << "phi(g) =  " << phi(g) << endl;
   
    vector<RingElem> VV;
    VV.push_back(x*y+power(x,2)+power(y,2)+power(z,2)+power(t,2)+t*z+h);
    VV.push_back(x*y+power(x,3)+power(y,3)+power(z,2)+power(t,2)+t*z+power(h,3));
    clog << "VV =  " << VV << endl;
    clog << "ChangeOrdering(VV) =  " << ChangeOrdering(NewP,VV) << endl;
*/
    GlobalManager CoCoAGlobals;
    ring QQ  =  RingQQ();
    vector<symbol> X  =  SymbolRange(symbol("x",0), symbol("x",4));
    SparsePolyRing P(NewPolyRing_DMPI(QQ, X));

    cout << "Executing test C4_h "<<endl;
    RingElem x  =  indet(P,0),  y  =  indet(P,1),  z  =  indet(P,2), t = indet(P,3);
    RingElem h  =  indet(P,4);
    PolyList InputPolys,GB,GB1,outMinGens;
    InputPolys.push_back(x+y+z+t);
    InputPolys.push_back(x*y+y*z+z*t+t*x);
    InputPolys.push_back(x*y*z+y*z*t+z*t*x+t*x*y);
    InputPolys.push_back(x*y*z*t-power(h,4));
    ComputeDynamicGBasis(GB, InputPolys);
    ComputeGBasis(GB1, outMinGens, InputPolys);
    monic(GB);
    cout << "C4_h  =  [" << GB1 << "]\n"<<endl;

/*
  RingElem f(P1);
  RingElem tmp1(power(x,3)*y*z);
  RingElem tmp2(x*y*t*power(h,2));

  clog << "\n c1  =  " << tmp1 << "\n"<<endl;
  clog << "\n c2  =  " << tmp2 << "\n"<<endl;
  clog << "\n inequality  =  " << CreateInquality(LPP(tmp1),LPP(tmp2))<< "\n"<<endl;
*/
  }
} // end of namespace CoCoA


// Use main() to handle any uncaught exceptions and warn the user about them.
int main()
{
  try
  {
    CoCoA::program();
    return 0;
  }
  catch (const CoCoA::ErrorInfo& err)
  {
    cerr << "***ERROR***  UNCAUGHT CoCoA Error";
    ANNOUNCE(cerr, err);
  }
  catch (const std::exception& exc)
  {
    cerr << "***ERROR***  UNCAUGHT std::exception: " << exc.what() << endl;
  }
  catch(...)
  {
    cerr << "***ERROR***  UNCAUGHT UNKNOWN EXCEPTION" << endl;
  }

  CoCoA::BuildInfo::PrintAll(cerr);
  return 1;
}
