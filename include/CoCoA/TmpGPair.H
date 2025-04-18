#ifndef CoCoA_GPair_H
#define CoCoA_GPair_H

//   Copyright (c)  2005, 2013  Massimo Caboara, Anna Bigatti

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


#include "CoCoA/SugarDegree.H"
#include "CoCoA/TmpGPoly.H"
#include "CoCoA/TmpGTypes.H"

#include <iosfwd>
//using std::ostream;
#include <list>
//using std::list;  // for GPairList
#include <vector>
//using std::vector;// for MGPairList
#include <cstddef>
//using std::size_t;


namespace CoCoA
{

  class GPair
  {
  public:
    //explicit GPair(const ring&, const PPMonoid&);

    // WARN i<j
    GPair(const GPoly& p);// special ctor for input polys
    GPair(const GPoly& p,  const GPoly& q);// usual ctor
    //    GPair(const GPoly* p,  const GPoly* q);// usualctor
    ~GPair(){};

    void myComplete(void); ///< fills to the GPair the field myOrdPoly and myWDeg (=deg(myOrdPoly))
    bool operator<(const GPair& gp) const;  ///< based on (w/std)deg + IsInputPoly + PPOrdering

    inline const GPoly& myFirstGPoly() const{return *myFirstGPolyPtr;};
    inline const GPoly& mySecondGPoly() const{return *mySecondGPolyPtr;};
    //  inline const GPoly* FirstGPolyPtr() const{return myFirstGPoly;};
    //  inline const GPoly* SecondGPolyPtr() const{return mySecondGPoly;};
    inline const GRingInfo& myGRingInfo() const {return myFirstGPoly().myGRingInfo();}

    inline int myFirstIndex() const{return age(myFirstGPoly());} ///< for stats
    inline int mySecondIndex() const{return age(mySecondGPoly());} ///< for stats
    inline bool IsInputPoly() const{return mySecondGPolyPtr==nullptr;}
    inline bool IamCoprime() const{return IamCoprimeFlag;}
    //    inline bool IsInputPoly() const{return mySecondIndex()==-1;}
    bool BCriterion_OK(const PPWithMask&)const;
    bool operator==(const GPair&)const;// uses the ordering on lcm_fg

    friend ConstRefPPMonoidElem LCM(const GPair&);
    friend const PPWithMask& LCMwMask(const GPair&);
    friend ConstRefRingElem OrdPoly(const GPair&);
    friend long GPairComponent(const GPair& GP);
    friend const degree& wdeg(const GPair& GP);
    //    friend long StdDeg(const GPair& GP);
    friend const SugarDegree& sugar(const GPair& GP);
    friend std::ostream& operator<<(std::ostream&, const GPair&);// the indices are printed ++
    friend const ring& CoeffRing(const GPair& P);
    friend const SparsePolyRing& owner(const GPair& P);
    friend const GRingInfo& GetGRingInfo(const GPair& P);

  private:
    const GPoly* myFirstGPolyPtr; ///< f
    const GPoly* mySecondGPolyPtr;///< g
    PPWithMask myLCMwMask; ///< lcm(lpp(f),lpp(g)), stored with DivMask
    RingElem myOrdPoly; ///< pairs will be ordered according to this polynomial (now = lcm)
    bool IamCoprimeFlag; ///< are lpp(f) and lpp(g) coprime?
    degree myWDeg; ///< degree of myLCM, hence of the pair, it being homog.
    long myComponent;///< LPP(f) and LPP(g) have the same component, this one
    //    std::size_t myStdDeg; ///< standard degree of myLCM.  Used only if GradingDim ==0 // ANNA: unused!
    SugarDegree mySugar;  ///< sugar of the pair
  };// end GPairs class

  inline long GPairComponent(const GPair& GP) {return GP.myComponent;}
  inline const degree& wdeg(const GPair& GP) {return GP.myWDeg;}
  //  inline long StdDeg(const GPair& GP) {return GP.myStdDeg;}  // ANNA: unused
  inline const SugarDegree& sugar(const GPair& GP) {return GP.mySugar;}
  inline const GRingInfo& GetGRingInfo(const GPair& P) {return P.myGRingInfo();}

  typedef std::list<GPair> GPairList;

//WARN the pairs are at the moment in ORD ordering
  void Ordered_Insert(GPairList&,const GPair&);
  void RemoveFromGPList(GPairList&,GPair&);


// This class is not used. See if it is still useful
  typedef std::vector<GPairList> MGPairList;
  class ModuleGPairList
  {

  public:
    ModuleGPairList();
    ~ModuleGPairList(){};
    void Insert(GPair&);
    friend std::ostream& operator<<(std::ostream&, const ModuleGPairList&);
    // if found, there is an element in this which divides P;, if !found, not.
    GPairList::iterator FindDivides(const GPair& P,bool& found);
    bool IsIn(const GPair& P);
    GPairList::iterator FindSameLCMAndSecondInd(const GPair& P,bool& found);
    //splice the part of degree D
    //   void GetPartInDegD(ModuleGPairList&, degree D);
    //add MGPL to this. suppose MGPL is all in a degree lesser than any in this
    void PushPartOfLesserDeg(ModuleGPairList&MGPL);
    long size() const;
  private:
    MGPairList myMGPList;
  };

  // -------------------- inline functions -------------------- //

  inline ConstRefPPMonoidElem LCM(const GPair& gp)
  {
    return PP(gp.myLCMwMask);
  }


  inline const PPWithMask& LCMwMask(const GPair& gp)
  {
    return gp.myLCMwMask;
  }


  inline ConstRefRingElem OrdPoly(const GPair& gp)
  {
    return gp.myOrdPoly;
  }

}// end namespace cocoa

#endif
