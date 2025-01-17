//   Copyright (c)  2008  Anna Bigatti

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

#include "CoCoA/TmpPPVector.H"
#include "CoCoA/SparsePolyOps-RingElem.H"
//#include "CoCoA/SparsePolyRing.H" // from SparsePolyOps-RingElem.H
#include "CoCoA/VectorOps.H"
//#include "CoCoA/PPWithMask.H"  -- included by CoCoA/TmpPPVector.H
//#include "CoCoA/DivMask.H"  -- included by CoCoA/PPWithMask.H
//#include "CoCoA/PPMonoid.H"  -- included by CoCoA/PPWithMask.H
//#include "CoCoA/ring.H"  -- included by CoCoA/SparsePolyRing.H

#include <algorithm>
using std::sort;
//#include <vector>
using std::vector;

// ints MFI_Occurrences, MFI_Indets;
// int Ints120[121];
// MixedPPVector  MTL120[121];
// PPVector GlobalSplitterTList;

/**********************************************************/
namespace CoCoA
{
  
  void PPVector::myInterreduce()
  { myInterreduceSort(); } // AMB: faster than without sorting


  void PushBack(PPVector& PPs, ConstRefPPMonoidElem pp)
  {
    if (PPM(PPs) != owner(pp))
      PPs.myPushBack(PPMonoidElem(PPM(PPs), exponents(pp)));
    else
      PPs.myPushBack(pp);
  }


  void PushBack(PPVector& PPs, const PPWithMask& pp)
  {
    CoCoA_ASSERT(PPM(PPs) == owner(PP(pp)));
    PPs.myPushBack(pp);
  }


  void PushBackPopBack(PPVector& ToPPs, PPVector& FromPPs)
  { ToPPs.myPushBackPopBack(FromPPs); }


  bool IsDivisible(const PPWithMask& pp, const PPVector& ByL)
  {
    // CoCoA_ASSERT(compatible)
    return ByL.myDivides(pp);
  }


  bool IsDivisible(ConstRefPPMonoidElem pp, const PPVector& ByL)
  {
    // CoCoA_ASSERT(compatible)
    return ByL.myDivides(pp);
  }


  void swap(PPVector& PPs1, PPVector& PPs2)
  { 
    // CoCoA_ASSERT(compatible)
    //    PPs1.mySwap(PPs2);
    swap(PPs1.myVec, PPs2.myVec);
  }


  void interreduce(PPVector& PPs)
  { PPs.myInterreduce(); }


  void InterreduceSort(PPVector& PPs)
  { PPs.myInterreduceSort(); }


  void lcms(PPVector& PPs, const PPVector& PPs1, const PPVector& PPs2)
  {
    // CoCoA_ASSERT(compatible)
    PPs.myLcms(PPs1, PPs2);
  }


  void convert(std::vector<RingElem>& v, ring P, const PPVector& PPs)
  {
    // CoCoA_ASSERT(compatible)
    PPs.myConvert(v, P);
  }
  

  void convert(PPVector& PPs, const std::vector<RingElem>& v)
  {
    PPVector tmpPPs(PPM(PPs), DMR(PPs));
    CoCoA_ASSERT(AreMonomials(v));
    for (const auto& f: v)  PushBack(tmpPPs, LPP(f));
    swap(tmpPPs, PPs);
  }


//   const PPWithMask& PPVector::operator[](MachineInt i) const
//   {
//     // CoCoA_ASSERT(range) ??? anna
//     //    return myAt(AsUnsignedLong(i)); ???
//     return myVec[AsUnsignedLong(i)];
//   }


  std::ostream& operator<<(std::ostream& out, PPVector PPs)
  {
    if (!out) return out;  // short-cut for bad ostreams
    out << PPs.myVec;
    return out;
  }


//   class PPVector
//   {
//   public:
//     PPVector(PPMonoid PPM, DivMaskRule DMR);
//     virtual ~PPVector() {};
//     long myNumIndetsEnv() const { return NumIndets(myPPM); }
//     long mySize() const { return len(myVec); }
//     void myReserve(long n);
//     bool IamEmpty() const;
//     PPVectorElem myDefaultElem() const;
//     const PPVectorElem myComp(unsigned long i) const;
//     void myPushBack(ConstRefPPMonoidElem pp);
//     void myPushBack(const PPWithMask& pp);
//     void myPushBackPopBack(PPVector& FromPPs);
//     void mySort();
//     bool myDivides(const PPWithMask& pp) const;
//     bool myDivides(ConstRefPPMonoidElem pp) const;
//     void myInterreduceSort();
//     void myLcms(const PPVector& PPs1, const PPVector& PPs2);
//     void myConvert(std::vector<RingElem>& v, SparsePolyRing P) const;
//   private:
//     const PPWithMask& myNth(int n) const;
//   private: // member fields
//     PPMonoid myPPM;
//     DivMaskRule myDMR;
//     vector<PPWithMask> myVec;
//   };
  

  //-------- auxiliary
  bool PPWithMaskLessThan(const PPWithMask& f, const PPWithMask& g)
  {
    return (owner(PP(f)))->myCmp(raw(PP(f)), raw(PP(g))) < 0;
  }
  
  //-------- auxiliary end

  
  PPVector::PPVector(PPMonoid PPM, DivMaskRule DMR):
    myPPM(PPM), myDMR(DMR)
  { }


  PPVector::PPVector(PPMonoid PPM, DivMaskRule DMR, const std::vector<RingElem>& v):
    myPPM(PPM), myDMR(DMR)
  {
    convert(*this, v);
  }


  PPVector::PPVector(const PPVector& PPs):
    myPPM(PPs.myPPM), myDMR(PPs.myDMR)
  {
    myVec = PPs.myVec;
  }


  PPVector& PPVector::operator=(const PPVector& rhs)
  {
    if (this == &rhs) return *this;
    PPVector copy(rhs);
    swap(*this, copy);
    return *this;
  }
  
  
//   PPVector* PPVector::myZeroClone() const
//   { return new PPVector(myPPM, myDMR); }
  
  
  void PPVector::myReserve(long n) // negative n check???
  { myVec.reserve(n); }
  
  
  void PPVector::myClear() // negative n check???
  { myVec.clear(); }
  
  
  bool PPVector::IamEmpty() const
  { return myVec.empty(); }
  
  
//   const PPWithMask PPVector::myAt(unsigned long i) const
//   {
//     // anna: make all checks
//     return myVec[i];
//   }
  
  
  void PPVector::myPushBack(ConstRefPPMonoidElem pp)
  { myVec.push_back(PPWithMask(pp, myDMR)); }
  
  
  void PPVector::myPushBack(const PPWithMask& pp)
  { myVec.push_back(pp); }
  
  
  void PPVector::myPushBackPopBack(PPVector& FromPPs)
  {
    myPushBack((FromPPs.myVec).back());
    FromPPs.myVec.pop_back();
  }

  // NOT exception clean
  void PPVector::mySort()
  {
    sort(myVec.begin(), myVec.end(), PPWithMaskLessThan);
  }
  
  
  bool PPVector::myDivides(const PPWithMask& pp) const
  {
//     const long n=myLen();
//     for (long i=0; i < n; ++i)
//       if ( IsDivisibleFast(pp, myVec[i]) ) return true;
    for (const auto& pp_i: myVec)
      if ( IsDivisibleFast(pp, pp_i) ) return true;
    return false;
  }
  
  
  bool PPVector::myDivides(ConstRefPPMonoidElem pp) const
  {
    return myDivides(PPWithMask(pp, myDMR));
  }
  
  
  // NOT exception clean
  void PPVector::myInterreduceSort()
  {
    if (IamEmpty())  return;
    mySort(); // ascending
    vector<PPWithMask> v;
    v.reserve(myLen());
    swap(myVec, v); //// ANNA: swap afterwards for exception safety?
    for (const auto& pp: v)
      if (!myDivides(pp)) myPushBack(pp);
  }


  // anna: careful with push_back and multiple copies of lcms
  void PPVector::myLcms(const PPVector& PPs1, const PPVector& PPs2)
  {
    // anna: this calls the destructor too many times (2*mylcm calls)
    // CoCoA_ASSERT(compatible)
    if (IsEmpty(PPs1) || IsEmpty(PPs2))
    {
      myVec.clear();
      return;
    }
    PPMonoidElem pp(myPPM);
    vector<PPWithMask> res;
    res.reserve(len(PPs1)*len(PPs2));
    for (const auto& pp1: PPs1.myVec)
      for (const auto& pp2: PPs2.myVec)
      {
        myPPM->myLcm(raw(pp), raw(PP(pp1)), raw(PP(pp2)));
        res.push_back(PPWithMask(pp, myDMR));
      }
    swap(myVec, res);
  }
  
  
  void PPVector::myAlexanderDual(const PPVector& PPs)
  {
    CoCoA_ASSERT(PPM(PPs) == myPPM);
    CoCoA_ASSERT(DMR(PPs) == myDMR);
    PPVector res(PPM(PPs), DMR(PPs));
    PPVector tmp(PPM(PPs), DMR(PPs));
    PushBack(res, one(PPM(PPs)));
    //    for (long i=0; i<len(PPs); ++i)
    for (const auto& ppwm: PPs.myVec)
    {
      tmp.mySupport(PP(ppwm));
      res.myLcms(res, tmp);
      res.myInterreduceSort();
    }
    swap(myVec, res.myVec);
  }


  void PPVector::mySupport(ConstRefPPMonoidElem pp)
  {
    CoCoA_ASSERT(owner(pp) == myPPM);
    const PPMonoid& PPM = myPPM;
    vector<PPWithMask> res;
    for (long i=0 ; i < NumIndets(PPM) ; ++i )
      if ( exponent(pp, i) != 0 )
        res.push_back(PPWithMask(indet(PPM, i), myDMR));
    swap(myVec, res);
  }


  void PPVector::myConvert(std::vector<RingElem>& v, SparsePolyRing P) const
  {
    vector<RingElem> tmp;
    if (myPPM == PPM(P))
      for (const auto& ppwm: myVec)
        tmp.push_back(monomial(P, PP(ppwm)));
    else
    {
      vector<long> expv(NumIndets(P));
      for (const auto& ppwm: myVec)
        tmp.push_back(monomial(P, exponents(expv,PP(ppwm))));
    }
    swap(v, tmp);
  }


/********************    list of terms    ********************/

/*************************************************************/

//#define Random(len)( 1 )

// void InsInTList(PPVector TL, eterm t, int *NewLen)
// {
//   ints OccInd = Indets(t);
//   if ( IntsGetLen(OccInd) == 1 )
//   {
//     InsInSPList (OccInd[1], eterm_get_nth(t,OccInd[1]), SPList(TL) );
//     eterm_free(t);
//   }
//   else
//     MTLPutLast (MTList(TL), *NewLen, t);
// }



/********    Simple Power List    ********/

// void InsInSPList(size_t Index, unsigned Exp, eterm SPL)
// {
//   size_t CurrentExp =SPL[Index];
  
//   if (CurrentExp == 0)
//   {  
//     eterm_put_non0_nth(SPL, Index, Exp);
//     IntsPutLast(Indets(SPL), Index);
//   }
//   else if (CurrentExp > Exp)
//     eterm_put_non0_nth(SPL, Index, Exp);
// }

// #define SPLDividesTerm(SPL,t)  sp_SmallMult(t,SPL) 

// bool TListSimplePowersDivide(PPVector PPs, eterm t)
// {
//   eterm SPL =SPList(PPs);
//   int   i;
  
//   for ( i=eterm_get_indetsNo(SPL) ; i>0 ; --i )
//     if ( SPL[i]!=0 && SPL[i]<=t[i] )   return TRUE;
//   return FALSE;
// }

// /************************    SPLIT    ************************/
 
// void DelLinearSimplePowers(PPVector PPs, int *LPSNo)
// {
//   eterm SPL = SPList(PPs);
//   ints  OccInd = Indets(SPL);
//   size_t   n = IntsGetLen(OccInd), sum = 0, index;
  
//   for ( ; n>0 ; --n )
//     if ( SPL[(index=OccInd[n])] == 1 )
//     {
//       ++sum;
//       eterm_put0_nth(SPL, index);
//       IntsMoveLastToNth(OccInd, n); 
//     }
//   *LPSNo = sum;  
// }

// eterm GetCoprimeSPList(PPVector PPs)
// {
//   int    MTLLen =MTListLen(PPs), *OccInd;
//   MixedPPVector   MTL =MTList(PPs);
//   size_t    i, n =MTLLen, index;
//   eterm res, theSPL; 
//   shortbitset CoprimeSupp;

//   if ( IntsGetLen(Indets(SPList(PPs))) == 0 )   return nullptr;
//   CoprimeSupp = SqFr(SPList(PPs));  
//   for ( ; n > 0 ; --n)
//   {
//     CoprimeSupp &=(~(SqFr(MTL[n]) ) );
//     if ( CoprimeSupp == 0 )   return nullptr;
//   }
//   theSPL = SPList(PPs);
//   OccInd = Indets(theSPL);  
//   res = eterm_init(eterm_get_indetsNo(theSPL));  
//   for ( i=IntsGetLen(OccInd) ; i>0 ; --i )
//   {
//     if ( bs_get_n(CoprimeSupp, (size_t)OccInd[i]) )
//     {
//       index = OccInd[i];
//       eterm_put_non0_nth(res, index, theSPL[index]);
//       IntsPutLast(Indets(res), index);
//       eterm_put0_nth(theSPL, index);      
//       IntsMoveLastToNth(OccInd, i);
//     }
//   }  
//   return res;
// }


// void MoveNotCoprimeSP(eterm FromSPList, eterm ToSPList, eterm theTerm)
// {
//   ints OccInd = Indets(FromSPList);
//   size_t   i = IntsGetLen(OccInd), index;
  
//   for ( ; i>0 ; --i )
//     if ( theTerm[(index = OccInd[i])] != 0 )
//     {
//       eterm_put_non0_nth(ToSPList, index, FromSPList[index]);
//       IntsPutLast(Indets(ToSPList), index);
//       eterm_put0_nth(FromSPList, index);
//       IntsMoveLastToNth(OccInd, i);
//     }
// }

// void MoveNotCoprime(PPVector FromTList, PPVector ToTList, eterm theTerm)
// {
//   int       FromMTLLen =MTListLen(FromTList),  ToMTLLen=0;
//   MixedPPVector   FromMTL =MTList(FromTList),  ToMTL =MTList(ToTList);
//   int     n =FromMTLLen;
  
//   for ( ; n > 0 ; --n)
//     if ( !eterm_coprime(FromMTL[n], theTerm) ) 
//     { 
//       ToMTL.push_back(FromMTL[n]);
//       FromMTL.myEraseNth(n);
//     }
  
//   MoveNotCoprimeSP(SPList(FromTList), SPList(ToTList), theTerm);
// }

// PPVector SplitIndets(PPVector PPs)
// {
//   PPVector      res =GlobalSplitterTList;
//   int           resLen = 0;
//   MixedPPVector MTL = MTList(PPs),  resMTL =MTList(res);
//   int   n = MTListLen(PPs),  m;
//   eterm tn, tm;

//   for ( ; n > 0 ; --n)
//   {
//     tn = MTL[n];
//     for ( m=resLen ; m>0 ; --m)
//     {
//       tm = resMTL[m];
//       if ( !eterm_coprime(tn, tm)  )
//       {
//         eterm_union_and_assign(tm, tn);
//         if ( tn != MTL[n] ) eterm_free(tn);
//         tn = tm;
//         resMTL.myEraseNth(m);
//       }  
//     }
//     if ( tn != MTL[n] )    MTLPutLast(resMTL, resLen, tn);
//     else MTLPutLast(resMTL, resLen, eterm_dup(tn));
//   }
//   if ( resLen == 1 )
//   {
//     eterm_free(resMTL[1]);
//     SetMTListLen(res, 0);
//     return nullptr;
//   }
//   SetMTListLen(res, resLen);
  
//   return res;
// }

// /************************    INTERREDUCE    ************************/

// void MTLOrderByDecrDegree(MixedPPVector theMTList, int Len, int MaxDeg)
// {
//   eterm t;
//   MixedPPVector *MTLDeg=(MixedPPVector*)mycalloc
//    ((MaxDeg+1),sizeof(MixedPPVector),"*MixedPPVector");
//   int *MTLLenDeg=(int*)mycalloc((MaxDeg+1), sizeof(int),"*int"), auxLen =0;
//   int i, d;

//   for ( d=MaxDeg ; d>=2 ; --d)    MTLDeg[d] = malloc_MTList(Len);

//   for ( i=Len ; i>0 ; --i)
//   {
//     d = eterm_degree((t = theMTList[i]) ); 
//     MTLPutLast( MTLDeg[d], MTLLenDeg[d], t);
//   }
//   for ( d=MaxDeg ; d>1 ; --d)
//   {
//     if (MTLLenDeg[d] != 0)
//       for ( i=MTLLenDeg[d]; i>0; --i)
//         MTLPutLast(theMTList, auxLen, (MTLDeg[d])[i]);
//     free_MTList(Len, MTLDeg[d]);
//   }
//   myfree((MaxDeg+1)*sizeof(MixedPPVector), MTLDeg, "*MixedPPVector");  
//   myfree((MaxDeg+1)*sizeof(int), MTLLenDeg, "*int");

// }

// bool MTLDividesTerm(MixedPPVector MTL, int MTLLen, eterm T)
// {
//   int i = MTLLen;

//   for ( ; i>0 ; --i ) if ( eterm_divides(MTL[i], T) )  return TRUE;

//   return FALSE;
// }

// void InterreduceTList(PPVector PPs)
// {
//   MixedPPVector   MTL =MTList(PPs);
//   int      MTLLen =MTListLen(PPs);
//   int      n =MTLLen, j, MaxDeg =0;
//   eterm    auxT;

//   for ( ; n> 0 ; --n)
//     if ( TListSimplePowersDivide(PPs,(auxT=MTL[n]) )  )
//       MTL.myEraseNth(n);
//     else
//       MaxDeg = MAX(MaxDeg, eterm_degree(auxT));
//   if (MTLLen > 1)
//   {
//     MTLOrderByDecrDegree(MTL, MTLLen, MaxDeg);
//     for ( n=MTLLen-1 ; n>0; --n)
//     {
//       auxT =  MTL[n];
//       for ( j = MTLLen ; j>n ; --j)
//         if ( eterm_divides(MTL[j], auxT)  )
//         {
//           MTL.myEraseNth(n);
//           break;
//         }
//     }
//   }
//   SetMTListLen(PPs, MTLLen);
// }

// /************************  PIVOT  ************************/

// int MostFrequentIndet(PPVector PPs)
// {
//   int IndNo = TListIndetsNo(PPs), exp, MFIndNo = 1, i, j;
//   ints OccInd;
  
//   for ( j=IndNo ; j>0 ; --j )    MFI_Occurrences[j] = 0;
//   for ( i=len(PPs) ; i>0 ; --i)
//   {
//     OccInd = Indets(MTL[i]);
//     for ( j=IntsGetLen(OccInd) ; j>0 ; --j )   MFI_Occurrences[OccInd[j]]++;
//   }
//   MFI_Indets[1] = (j=IndNo);
//   i = MFI_Occurrences[j];
//   for ( --j ; j>0 ; --j )
//     if ((exp=MFI_Occurrences[j]) >= i)
//     {
//       if (exp > i)  { MFIndNo = 0;  i = exp; }
//       MFI_Indets[++MFIndNo] = j;
//     }  
 
//   if ( i == 1 )       return -1;
//   if ( MFIndNo == 1 ) return MFI_Indets[1];
//   else                return MFI_Indets[MFIndNo/2];
// }


// void PPVector::myPivotSetIndet()
// {
//   int MFIndet = MostFrequentIndet(PPs);

//   AssignOne(myPivot());
//   if (MFIndet == -1)  return;

//   myPivot() = indet(myPivot().myPPM(), MFIndet);
// }

// void PPVector::myPivotSetGCD3()
// {
//   int MFIndet = MostFrequentIndet(PPs);

//   AssignOne(myPivot());
//   if (MFIndet == -1)  return;
 
//   // reverse order for compatibility with old C code
//   for (int i=len(PPs), int count=0; --i>=0 && count<3; )
//     if ( PPs[i][MFIndet] != 0 )
//     {
//       myPivot() = gcd(myPivot(), PPs[i]);
//       count++;
//     }
// }

// void PPVector::myPivotSetSimplePower()
// {  
//   int MFIndet = MostFrequentIndet(PPs);

//   AssignOne(myPivot());
//   if (MFIndet == -1)  return;
 
//   int r = -1;
//   while (myPPs[++r][MFIndet] == 0);
//   int e1 = myPPs[r][MFIndet];
//   r = len(myPPs);
//   while (myPPs[--r][MFIndet] == 0);
//   int e2 = myPPs[r][MFIndet];

//   myPivot() = IndetPower(myPivot().myPPM(), MFIndet, min(e1, e2));
// }

/*
#define PivotOf(PPs,MTLLen)\
    ((MTLLen<3)?GCD3PivotOf(PPs):BigPivotOf(PPs))
*/

// #define PivotOf(PPs,MTLLen) BigPivotOf(PPs)

// eterm eterm_colon_SP(eterm theTerm, size_t index, unsigned exp)
// {
//   unsigned OldExp;
 
//   if ((OldExp =theTerm[index]) != 0 )
//   {
//     if ( OldExp <= exp ) printf("OldExp <= exp");
//     eterm_put_non0_nth(theTerm, index, OldExp-exp );
//   }
//   return theTerm;
// }

// void OccIndDel(eterm T, int index)
// {
//   ints    OccInd = Indets(T);
//   int     n = IntsGetLen(OccInd);

//   while ( OccInd[n] != index )  --n;
//   IntsMoveLastToNth(OccInd, n);
// }

// void ReduceAndDivideBySimplePower( PPVector PPs, PPVector *DivTList,
//                                    size_t PIndex, unsigned PExp)
// {
//   int MTLLen = MTListLen(PPs), auxLen = MAX(MTListLen(PPs),199);
//   int DivMTLLen, MTLLenEe, i;
//   unsigned TExp, e, index;
//   ints DivMTLLenExp;
//   eterm  T;
//   MixedPPVector   *DivMTLExp, MTL =MTList(PPs), DivMTL, MTLEe;
//   eterm   DivSPL;
  
//   if ( PExp>119 )
//   {    
//     DivMTLLenExp = (int*)mymalloc((PExp+2)*sizeof(int), "*int");
//     DivMTLExp = (MixedPPVector*)mymalloc((PExp+2)*sizeof(MixedPPVector),
//                                          "*MixedPPVector");
//   }
//   else
//   {  DivMTLLenExp = Ints120;     DivMTLExp = MTL120;   }
//   for ( i=PExp+1 ; i>=0 ; --i )
//   {  DivMTLExp[i] = malloc_MTList(auxLen);    DivMTLLenExp[i] = 0;  }
  
//   /* DivTList */
//   DivSPL = eterm_colon_SP(eterm_dup(SPList(PPs)), PIndex, PExp);
//   /* PPs */
//   InsInSPList(PIndex, PExp, SPList(PPs));
 
//   for ( i=MTLLen ; i > 0 ; --i)
//   {
//     TExp = (T=MTL[i])[PIndex];
//     if ( TExp > PExp )
//     {
//       PPsSwapLastAndNth(MTL, MTLLen, i);
//       eterm_put_non0_nth(T, PIndex, TExp-PExp);
//       MTLPutLast(DivMTLExp[PExp+1], DivMTLLenExp[PExp+1], T);
//     }
//     else if ( TExp == PExp )
//     {
//       PPsSwapLastAndNth(MTL, MTLLen, i);
//       if ( eterm_get_OccIndNo(T) == 2 )
//       {
// 	if ((index = (Indets(T))[1]) == PIndex )  index = (Indets(T))[2];
//         InsInSPList(index, T[index], DivSPL);
// 	eterm_free(T);
//       }
//       else
//       {
// 	eterm_put0_nth(T, PIndex);
// 	OccIndDel(T, PIndex);
//         MTLPutLast(DivMTLExp[TExp], DivMTLLenExp[TExp], T);
//       }
//     }
//     else if ( TExp != 0 )
//     {
//       if ( eterm_get_OccIndNo(T) == 2 )
//       {
// 	if ((index = (Indets(T))[1]) == PIndex )  index = (Indets(T))[2];
//         InsInSPList(index, T[index], DivSPL);
//       }
//       else
//       {
// 	eterm_put0_nth(T, PIndex);
// 	OccIndDel(T, PIndex);
//         MTLPutLast(DivMTLExp[TExp], DivMTLLenExp[TExp], T);
//       }
//     }
//     else /* TExp == 0 */ MTLPutLast(DivMTLExp[TExp], DivMTLLenExp[TExp], T);
//   }
 
//   /* DivTList(Interreduction) */
//   DivMTL    = DivMTLExp[PExp];
//   DivMTLLen = DivMTLLenExp[PExp];
//   for ( e=PExp-1 ; e>0 ; --e )
//   {
//     MTLEe    = DivMTLExp[e];
//     MTLLenEe = DivMTLLenExp[e];
//     for ( i=MTLLenEe ; i>0 ; --i )
//     {
//       T = MTLEe[i];
//       if ( SPLDividesTerm(DivSPL, T) )  PPsSwapLastAndNth(MTLEe, MTLLenEe, i);
//       else if ( MTLDividesTerm(DivMTL, DivMTLLen, T) )
//         PPsSwapLastAndNth(MTLEe, MTLLenEe, i);
//       else
//         MTLPutNth(MTLEe, i, eterm_dup(T));
//       eterm_put_non0_nth(T, PIndex, e);
//       IntsPutLast(Indets(T), PIndex);
//     }
//     if ( MTLLenEe!=0 )    MTLAppend(DivMTL, &DivMTLLen, MTLEe, MTLLenEe);
//     free_MTList(auxLen, MTLEe);
//   }
//   MTLEe    = DivMTLExp[0];
//   MTLLenEe = DivMTLLenExp[0];
//   for ( i=MTLLenEe ; i>0 ; --i )
//   {
//     T = MTLEe[i];
//     if ( SPLDividesTerm(DivSPL, T) )   PPsSwapLastAndNth(MTLEe, MTLLenEe, i);
//     else if ( MTLDividesTerm(DivMTL, DivMTLLen, T) )
//       PPsSwapLastAndNth(MTLEe, MTLLenEe, i);
//     else
//       MTLPutNth(MTLEe, i, eterm_dup(T));
//   }
//   if ( MTLLenEe!=0 )    MTLAppend(DivMTL, &DivMTLLen, MTLEe, MTLLenEe);
//   free_MTList(auxLen, MTLEe);
  
//   MTLAppend(DivMTL, &DivMTLLen, DivMTLExp[PExp+1], DivMTLLenExp[PExp+1]);
//   free_MTList(auxLen, DivMTLExp[PExp+1]);
//   *DivTList = TListMake(DivSPL, DivMTL, auxLen, DivMTLLen);
 
//   /* PPs */
//   TListReduceSize(PPs, MTLLen);
//   if ( PExp>119 )
//   {
//     myfree((PExp+2)*sizeof(int), DivMTLLenExp, "*int");
//     myfree((PExp+2)*sizeof(MixedPPVector), DivMTLExp, "*MixedPPVector");
//   }
// }

// void ReduceAndDivideByMixedTerm( PPVector PPs, PPVector *DivTList,
//                                  eterm Pivot)
// {
//   int  MTLLen =MTListLen(PPs), OldMTLLen =MTLLen,
//        BMLen =0, CMTLLen =0, DivMTLLen =0;
//   int index, TDeg;
//   eterm     DivT,  t, auxSPL;
//   int     i = MTLLen;
//   MixedPPVector   MTL =MTList(PPs),  BigMultMTL, CoprimeMTL, DivMTL;
//   eterm   DivSPL;
 
//   *DivTList = NewTList(MTLLen, TListIndetsNo(PPs));
//   DivMTL = MTList(*DivTList);
//   DivSPL = eterm_colon(eterm_dup(SPList(PPs)), Pivot);
//   auxSPL =(SPList(*DivTList));  /* eterm_init() */

//   BigMultMTL = malloc_MTList(MTLLen);
//   CoprimeMTL = malloc_MTList(MTLLen);
 
//   for ( ; i > 0 ; --i )
//   {
//     t=MTL[i];
//     TDeg = eterm_degree(t);

//     if ( eterm_divides(Pivot,t) )
//     {
//       /* PPs */
//       PPsSwapLastAndNth(MTL, MTLLen, i);
//       /* DivTList */
//       if (sp_BigMult(t,Pivot))
// 	MTLPutLast(BigMultMTL, BMLen, eterm_colon(t,Pivot));
//       else
//       {
// 	DivT = eterm_colon(t, Pivot);
// 	if ( IntsGetLen(Indets(DivT)) == 1)
// 	{
// 	  index =  (Indets(DivT))[ 1];
// 	  InsInSPList(index, DivT[index], DivSPL);
// 	  InsInSPList(index, DivT[index] + Pivot[index], auxSPL);
// 	  eterm_free(DivT);
// 	}
// 	else
// 	  MTLPutLast(DivMTL, DivMTLLen, DivT);
//       }
//     }
//     else if ( SPLDividesTerm(auxSPL, t) )
//       ;
//     else if (TDeg==eterm_degree((DivT=eterm_colon(eterm_dup(t),Pivot) )) )
//       MTLPutLast(CoprimeMTL, CMTLLen, DivT);
//     else if ( IntsGetLen(Indets(DivT)) == 1 )
//     {
//       index =  (Indets(DivT))[ 1];
//       InsInSPList(index, DivT[index], DivSPL);
//       InsInSPList(index, DivT[index] + Pivot[index], auxSPL);
//       eterm_free(DivT);
//     }
//     else
//       MTLPutLast(DivMTL, DivMTLLen, DivT);
//   }
//   /* PPs */
//   MTLPutLast(MTL, MTLLen, Pivot);
//   TListReduceSize(PPs, MTLLen);
//   /* DivTList */
//   eterm_colon(auxSPL, Pivot);
//   SetMTListLen(*DivTList, DivMTLLen);
//   InterreduceTList(*DivTList);
//   DivMTLLen = MTListLen(*DivTList);

//   if ( DivMTLLen != 0 )

//     for ( i=CMTLLen ; i>0 ; --i )
//     {
//       if ( SPLDividesTerm(DivSPL,(t=(CoprimeMTL)[i]) ) )
//       {
// 	eterm_free(t);
// 	PPsSwapLastAndNth(CoprimeMTL, CMTLLen, i);
//       }
//       else if ( MTLDividesTerm(DivMTL, DivMTLLen, t) )
//       {
// 	eterm_free(t);
// 	PPsSwapLastAndNth(CoprimeMTL, CMTLLen, i);
//       }
//     }
//   else
//     for ( i=CMTLLen ; i>0 ; --i )
//       if ( SPLDividesTerm(DivSPL,(t=(CoprimeMTL)[i]) ) )
// 	PPs.myEraseNth(i);
//   SPList(*DivTList) = DivSPL;
//   if (CMTLLen!=0)  MTLAppend(DivMTL, &DivMTLLen, CoprimeMTL, CMTLLen);
//   if (BMLen!=0)    MTLAppend(DivMTL, &DivMTLLen, BigMultMTL, BMLen);
//   free_MTList(OldMTLLen, CoprimeMTL);
//   free_MTList(OldMTLLen, BigMultMTL);
//   eterm_free(auxSPL);
//   TListReduceSize(*DivTList, DivMTLLen);
// }

// void ReduceAndDivideByPivot(PPVector PPs, PPVector *DivTList, eterm Pivot)
// {
//   if ( eterm_get_OccIndNo(Pivot) == 1)
//   {
//     int index = (Indets(Pivot))[1], exp = eterm_degree(Pivot);
//     eterm_free(Pivot);
//     ReduceAndDivideBySimplePower(PPs, DivTList, index, exp);
//   }
//   else
//     ReduceAndDivideByMixedTerm(PPs, DivTList, Pivot);
// }


  //----------------------------------------------------------------------
  
} // end of namespace CoCoA
