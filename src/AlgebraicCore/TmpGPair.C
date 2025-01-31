//   Copyright (c)  2005  John Abbott,  Anna M. Bigatti
//   Original author: 2005  Massimo Caboara

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

#include "CoCoA/TmpGPair.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/SparsePolyOps-RingElem.H"

#include <algorithm>
using std::find;
using std::find_if;
using std::max;
#include <iostream>
using std::ostream;
using std::endl;

namespace CoCoA
{

  // -- ANNA: temporary
  SugarDegree NewSugar(const GPair& gp)
  {
    PPMonoidElem cofactor1 = LPP(OrdPoly(gp))/LPPForOrd(gp.myFirstGPoly());
    //PP(myLCMwMask)/LPPForDiv(myFirstGPoly());
    PPMonoidElem cofactor2 = LPP(OrdPoly(gp))/LPPForOrd(gp.mySecondGPoly());
    SugarDegree s(sugar(gp.myFirstGPoly()));
    s->myMul(cofactor1);
    s->myUpdate(cofactor2, gp.mySecondGPoly());
    return s;
  }



// GPair ////////////////////////////////////////////////////////////////////
 // Special pair: it represents an input polynomial 
  GPair::GPair(const GPoly& the_p):
    myLCMwMask(the_p.myGRingInfo().myPPM(), the_p.myGRingInfo().myDivMaskRule()),
    myOrdPoly(monomial(owner(the_p), LPPForOrd(the_p))),
    myWDeg(wdeg(myOrdPoly)),
    mySugar(sugar(the_p))
  {
    myFirstGPolyPtr = &the_p;
    mySecondGPolyPtr = nullptr;	
    myLCMwMask = LPPForDivwMask(the_p);
    IamCoprimeFlag = false;
    myComponent = component(the_p);
    //    myStdDeg=0;
 }//GPair::GPair


  GPair::GPair(const GPoly& p, const GPoly& q):
    myLCMwMask(p.myGRingInfo().myPPM(), p.myGRingInfo().myDivMaskRule()),
    myOrdPoly(one(owner(p))),// Fake, filled by myComplete
    myWDeg(GradingDim(owner(p))),// Fake, filled by myComplete
    mySugar(uninitialized)// Fake, filled by myComplete
  {
    myFirstGPolyPtr = &p;
    mySecondGPolyPtr = &q;
    myLCMwMask = lcm(LPPForDiv(p), LPPForDiv(q));
    IamCoprimeFlag = IsCoprime(LPPForDiv(p), LPPForDiv(q));
    myComponent = component(p);//had to be : p.Component()==q.Component()
    //    myStdDeg=0;// Fake, filled by myComplete
  }

  // ANNA: not called
//   GPair::GPair(const GPoly* p,const GPoly* q,  unsigned int i, unsigned int j):
//     myLCMwMask(p->myGRingInfo().myPPM(), p->myGRingInfo().myDivMaskRule()),
//     myOrdPoly(one(owner(*p))),// Fake, filled by myComplete
//     myWDeg(GradingDim(owner(*p))),// Fake, filled by myComplete
//     mySugar(uninitialized)// Fake, filled by myComplete
//   {
//     myFirstGPolyPtr=p;
//     mySecondGPolyPtr=q;
//     myLCMwMask = lcm(LPPForDiv(*p), LPPForDiv(*q));
//     //    IamCoprimeFlag = (PP(myLCMwMask) == LPP(*p)*LPP(*q) );// MOD Aggiustare
//     IamCoprimeFlag = IsCoprime(LPPForDiv(*p), LPPForDiv(*q));
//     std::cout << "GPair*: i = " << i << " " << Age(*p) << " = Age" << std::endl;
//     std::cout << "GPair*: j = " << j << " " << Age(*q) << " = Age" << std::endl;
//     myFirstIndexValue=i;
//     mySecondIndexValue=j;	
//     myComponent=Component(*p);//had to be : p->Component()==q->Component()
//     //    myStdDeg=0;// Fake, filled by myComplete
//  }


  void GPair::myComplete()
  {
    const SparsePolyRing P = owner(myOrdPoly);
    myOrdPoly = monomial(P, 1, exponents(PP(myLCMwMask)));
    myWDeg = wdeg(PP(myLCMwMask));
    mySugar = NewSugar(*this);
    //  myStdDeg=GRI.TmpStdDeg(PP(myLCMwMask));  // ANNA is this used??
  }//myComplete


bool GPair::operator==(const GPair& P)const
{
//   return (this->myFirstIndexValue==P.myFirstIndexValue
//           &&
//           this->mySecondIndexValue==P.mySecondIndexValue);
// check on pointers
   return (myFirstGPolyPtr == P.myFirstGPolyPtr
           &&
           mySecondGPolyPtr == P.mySecondGPolyPtr);
}//operator==


ostream& operator<<(ostream& out, const GPair& P)
{
  if (!out) return out;  // short-cut for bad ostreams

  out << "<" << age(*(P.myFirstGPolyPtr));//P.myFirstIndexValue;
  if (P.IsInputPoly()) out <<",InputPoly";
  else out << "," << age(*(P.mySecondGPolyPtr));//P.mySecondIndexValue;
  if (P.myComponent!=0) out<<", Comp="<<P.myComponent;
  out << ", " << PP(P.myLCMwMask);
  if (P.IamCoprimeFlag)   out << ", coprime";
  out << ", Deg="<< P.myWDeg;
  out << ", Sugar="<< P.mySugar;
  //out<<", FirstPoly="<< poly(*(P.myFirstGPolyPtr));
  //out<<", FirstPoly Deg="<< wdeg(*(P.myFirstGPolyPtr));
  //out<<", OPType "<<GetGRingInfo(P).myInputAndGrading();
  out << ">";
  return out;
}//operator<<


void Ordered_Insert(GPairList& L,const GPair& P)
{
  GPairList::iterator it = find_if(L.begin(),L.end(),[&](const GPair& arg) { return (P < arg); });
  L.insert(it,P);
}//Ordered_Insert


void RemoveFromGPList(GPairList& L,GPair& P)
{
  GPairList::iterator it = find(L.begin(),L.end(),P);
  if (it != L.end())
    L.erase(it);
}//RemoveFromGPList


const ring& CoeffRing(const GPair& P)
{
  return CoeffRing(*(P.myFirstGPolyPtr));
}//CoeffRing

const SparsePolyRing& owner(const GPair& P)
{
  return owner(*(P.myFirstGPolyPtr));
}//owner


// ***************************************************  ModuleGPairsList

  namespace // anonymous
  {

    // *********** These 3 functions are strictly for use in the procedures below.
    // JAA: 2022-02-18  switched to "lambda fns"
    bool GPairDividesLCM(const GPair& P1, const GPair& P2)
    {
      return IsDivisibleFast(LCMwMask(P2), LCMwMask(P1)) && GPairComponent(P1)==GPairComponent(P2);
    }


    // Currently unused!  (comment out to avoid compiler warning)
    // bool GPairDividesLCMProperly(const GPair& P1, const GPair& P2)
    // {
    //   return IsDivisibleFast(LCMwMask(P2), LCMwMask(P1))
    //     && (LCMwMask(P1)!=LCMwMask(P2)) &&
    //     GPairComponent(P1)==GPairComponent(P2);
    // }


    bool GPairEqualLCM(const GPair& P1, const GPair& P2)
    {
      return LCMwMask(P1)==LCMwMask(P2) && GPairComponent(P1)==GPairComponent(P2);
    }

  } // end of namespace anonymous

// ***********

ModuleGPairList::ModuleGPairList()
{
  myMGPList.resize(10000);
}

void ModuleGPairList::Insert(GPair& P)
{
  Ordered_Insert(myMGPList[P.mySecondIndex()],P);
  //make a copy of en element that will be destroyed - think about using splice
}

GPairList::iterator ModuleGPairList::FindDivides(const GPair& P, bool& found)
{
  GPairList::iterator it;
  it = find_if(myMGPList[P.mySecondIndex()].begin(),
               myMGPList[P.mySecondIndex()].end(),
               [&](const GPair& arg) { return GPairDividesLCM(arg,P); });
  found = (it!=myMGPList[P.mySecondIndex()].end());
  return it;
}//FindDivides



bool ModuleGPairList::IsIn(const GPair& P)
{
  return find(myMGPList[P.mySecondIndex()].begin(),
              myMGPList[P.mySecondIndex()].end(),
              P)
    !=
    myMGPList[P.mySecondIndex()].end();
}//IsIn


GPairList::iterator ModuleGPairList::FindSameLCMAndSecondInd(const GPair& P,bool& found)
{
  GPairList::iterator it;
  it=find_if(myMGPList[P.mySecondIndex()].begin(),
             myMGPList[P.mySecondIndex()].end(),
             [&](const GPair& arg) { return GPairEqualLCM(arg,P);});
  found = (it!=myMGPList[P.mySecondIndex()].end());
  return it;
}//FindSameLCMAndSecondInd


long ModuleGPairList::size() const
{
  long S=0;
  for (int i=0; i!=len(myMGPList);++i)
    S += len(myMGPList[i]);
  return S;
}



ostream& operator<<(ostream& out,const ModuleGPairList& MGPL)
{
  if (!out) return out;  // short-cut for bad ostreams

  //out<<endl<<"SIZE"<<len(MGPL.myMGPList)<<endl;
  for (const auto& mgp: MGPL.myMGPList)
  {
    //out<<".";
    if (mgp.empty()) continue;
    GPairList::const_iterator it1=mgp.begin();
    out<<endl<<"Index="<<it1->mySecondIndex()<<endl;
    for (;it1!=mgp.end();++it1) {out<<*it1;};
  }
  //  MGPairList::const_iterator it;
  // for (it=MGPL.myMGPList.begin();it!=MGPL.myMGPList.end();++it)
  // {
  //   //out<<".";
  //   if (!it->empty())
  //   {
  //     GPairList::const_iterator it1=it->begin();
  //     out<<endl<<"Index="<<it1->mySecondIndex()<<endl;
  //     for (;it1!=it->end();++it1) {out<<*it1;};
  //   }
  // }
  return out;
}//operator<<


/********** End Module GPLists ********************************************/

// ex inline //

bool GPair::BCriterion_OK(const PPWithMask& NewPP)const
{
  if (IsInputPoly()) return true;	
  return  !(IsDivisibleFast(myLCMwMask, NewPP)
            &&(lcm(LPPForDiv(*mySecondGPolyPtr),PP(NewPP)) != PP(myLCMwMask))
            &&(lcm(LPPForDiv(*myFirstGPolyPtr),PP(NewPP)) != PP(myLCMwMask)));
}


///*
//  You have the choice of DegRevLex (Default) and the ring ordering
//  GPair ordering is SparsePolyRing (owner(myOrdPoly)) ordering
bool GPair::operator<(const GPair& the_gp)const
{
  int CMP;
  if ( (CMP=cmp(mySugar, the_gp.mySugar))!=0 )   return CMP<0;
  //  if (! IsConstWSugar(mySugar) )  // Anna: avoid useless computation in homog case
  if ( (CMP=FastCmp(myWDeg, the_gp.myWDeg))!=0 ) return CMP<0;
  if (IsInputPoly()!=the_gp.IsInputPoly())
    return the_gp.IsInputPoly();
  return SparsePolyRingPtr(owner(myOrdPoly))->myCmpLPP(raw(myOrdPoly), raw(the_gp.myOrdPoly))<0;
}

}// end namespace cocoa
