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
    friend std::ostream& operator<<(std::ostream&, const GPair&);// the indexes are printed ++
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



// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/TmpGPair.H,v 1.16 2022/02/18 15:37:53 abbott Exp $
// $Log: TmpGPair.H,v $
// Revision 1.16  2022/02/18 15:37:53  abbott
// Summary: Changed args to const
//
// Revision 1.15  2022/02/18 14:11:50  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.14  2019/03/18 17:56:49  abbott
// Summary: Changed 0 into nullptr where appropriate
//
// Revision 1.13  2017/04/18 09:29:31  bigatti
// --removed first/second index (the information is in age(poly))
//
// Revision 1.12  2014/04/30 16:02:23  abbott
// Summary: Replaced unsigned int by long (rtn type of mem fn size)
// Author: JAA
//
// Revision 1.11  2014/01/28 11:03:18  bigatti
// -- changed names myFirstGPolyValue --> myFirstGPolyPtr
//
// Revision 1.10  2013/10/28 13:15:04  bigatti
// -- IsSpecial --> IsInputPoly
//
// Revision 1.9  2011/05/26 16:31:53  bigatti
// -- changed unsigned int --> long
//
// Revision 1.8  2009/11/20 15:51:07  bigatti
// -- removed unused myStdDeg
//
// Revision 1.7  2009/10/27 17:15:14  bigatti
// -- fixed: using sugar(g)->myWSugar() insted of wsugar(g)
//
// Revision 1.6  2009/02/05 15:21:22  bigatti
// -- split inline/friend function definitions to please the new compiler
//
// Revision 1.5  2008/09/19 13:33:42  bigatti
// -- added: Sat algorithm (M.Caboara)
//
// Revision 1.4  2008/09/16 15:03:43  bigatti
// -- added LPPForDiv
// -- changed LPP into LPPForOrd
//
// Revision 1.3  2007/12/05 11:06:24  bigatti
// -- changed "size_t StdDeg/myStdDeg(f)" into "long"  (and related functions)
// -- changed "log/myLog(f, i)" into "MaxExponent/myMaxExponent(f, i)"
// -- fixed bug in "IsOne(ideal)" in SparsePolyRing.C
//
// Revision 1.2  2007/10/30 17:14:11  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1  2007/03/09 18:56:56  bigatti
// -- added Tmp prefix to Groebner related files
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.15  2007/03/08 18:42:06  cocoa
// Cleaned up whitespace.
//
// Revision 1.14  2007/03/07 22:24:03  bigatti
// -- reintroduced TmpGTypes.H (waiting for a better solution)
//
// Revision 1.13  2007/03/07 17:04:31  cocoa
// -- several changes by M.Caboara: more operations on ideals,
//    exception cleaner, coding conventions, WSugar, dynamic
//
// Revision 1.12  2006/12/06 17:11:16  cocoa
// -- removed #include "config.H"
//
// Revision 1.11  2006/12/04 13:55:54  cocoa
// -- added: sugar for GradingDim > 0  (called wsugar)
//
// Revision 1.10  2006/11/24 17:17:58  cocoa
// -- reorganized includes of header files
//
// Revision 1.9  2006/11/22 14:43:32  cocoa
// -- minor cleaning (indicated by Intel compiler)
//
// Revision 1.8  2006/11/20 14:57:17  cocoa
// -- added: (standard) sugar for modules
// -- fixed: non-homogeneous sysygies
// -- minor fixes     [M.Caboara]
//
// Revision 1.7  2006/10/11 13:33:06  cocoa
// -- rearranged code for sugar in reduce.C
// -- activated sugar in GPair.C
// -- removed operator<< for GPairList (template in io)
//
// Revision 1.6  2006/10/06 14:04:16  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.5  2006/08/17 09:24:17  cocoa
// -- added: sugar
// -- changed: coding conventions
//
// Revision 1.4  2006/06/20 17:19:45  cocoa
// -- moved  GPair::operator<  into .C file
//
// Revision 1.3  2006/06/09 16:21:47  cocoa
// -- small adjustments on StdDeg computations
//
// Revision 1.2  2006/06/09 16:06:04  cocoa
// -- myStdDeg computed only if GradingDim==0
//
// Revision 1.1.1.1  2006/05/30 11:39:36  cocoa
// Imported files
//
// Revision 1.12  2006/04/27 14:02:53  cocoa
// -- tidied up include files
//
// Revision 1.11  2006/04/21 16:12:51  cocoa
// -- just commented code
//
// Revision 1.10  2006/04/11 15:04:44  cocoa
// -- added: StdDeg(gp)
//
// Revision 1.9  2006/04/11 14:54:21  cocoa
// -- added: myStdDeg
//
// Revision 1.8  2006/04/10 17:02:47  cocoa
// -- changed: BCriterion_OK now uses PPWithMask instead of ConstRefPPMonoidElem
//
// Revision 1.7  2006/03/17 18:28:01  cocoa
// -- changed: ordering on pairs (depending on current ordering)
//
// Revision 1.6  2006/03/12 21:28:34  cocoa
// Major check in after many changes
//
// Revision 1.5  2006/03/02 13:44:39  cocoa
// -- just comments
//
// Revision 1.4  2006/02/13 13:41:22  cocoa
// -- fixed: "const RingElem&" --> "ConstRefRingElem"
//
// Revision 1.3  2006/01/20 15:43:30  cocoa
// -- fixed: use of RefPPMonoidElem and ConstRefPPMonoidElem
//
// Revision 1.2  2006/01/17 15:44:56  cocoa
// -- chamges by Max for operations with modules
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.3  2005/09/29 15:49:58  cocoa
// -- fixed bugs in "operator<"
//
// Revision 1.2  2005/07/01 16:08:16  cocoa
// Friday check-in.  Major change to structure under PolyRing:
// now SparsePolyRing and DUPolyRing are separated (in preparation
// for implementing iterators).
//
// A number of other relatively minor changes had to be chased through
// (e.g. IndetPower).
//
// Revision 1.1.1.1  2005/05/03 15:47:30  cocoa
// Imported files
//
// Revision 1.2  2005/04/19 14:06:04  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.9  2004/11/19 16:15:51  cocoa
// (a) Removed unused error message about degree of zero;
//     replaced it by a more general message about needing a
//     non-zero polynomial (for various ops such as LC, LPP).
// (b) Added some missing arg checking in LC, LPP and deg
//     (for elems of a PolyRing).
// (c) Updated some commented out code in GPair and GPoly.
//
// Revision 1.8  2004/11/02 15:55:18  cocoa
// -- change: "deg" now returns a const reference
//
// Revision 1.7  2004/10/29 15:51:37  cocoa
// -- changed myLCM into myLCMwMask (PPMonoidElem --> PPWithMask)
// -- function IsDivisible had wrong semantics --> swapped arguments everywhere
//
// Revision 1.6  2004/06/16 16:13:41  cocoa
// Improved I/O facilities with knock-on changes
//
// Revision 1.5  2003/11/14 14:10:16  cocoa
// Max - New Pair ordering
//
// Revision 1.4  2003/11/12 16:41:06  cocoa
// CVS Test by Max
//
// Revision 1.3  2003/10/09 12:48:17  cocoa
// New coding convention for rings.
//
// Revision 1.2  2003/10/01 10:35:32  cocoa
// - applied "my" coding convention to PPMonoid and PPOrdering
//
// Revision 1.1.1.1  2003/09/24 12:55:43  cocoa
// Imported files
//
// Revision 1.10  2003/09/22 17:20:47  bigatti
// - new field myOrdPoly to order the GPairs
//
// Revision 1.9  2003/06/23 17:10:34  abbott
// Minor cleaning prior to public release.
// Improved the include directives,
//
// Revision 1.8  2003/05/29 16:12:33  bigatti
// - fix: unsigned int --> int  for RingSpecialIndex in Component(..)
//
// Revision 1.7  2003/05/28 14:20:07  bigatti
// - new code for modules
//
// Revision 1.6  2003/05/14 16:37:39  bigatti
// - myDeg is now of type degree
// - new ring syntax
//
// Revision 1.5  2002/09/19 17:21:33  bigatti
// - Cleaner code based on PolyRing
//
// Revision 1.4  2002/04/15 17:14:58  bigatti
// - Max's new code
//
// Revision 1.3  2002/04/09 14:08:15  bigatti
// - SPoly now takes a GPair as argument
//
// Revision 1.2  2001/12/12 18:22:53  bigatti
// - new structure of reduction
//
// Revision 1.1  2001/12/05 12:57:13  bigatti
// Initial revision
//

#endif
