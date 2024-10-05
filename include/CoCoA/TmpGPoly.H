#ifndef CoCoA_GPoly_H
#define CoCoA_GPoly_H

//   Copyright (c)  2005-2017 John Abbott, Anna M. Bigatti
//   Author: 2005-2010  Massimo Caboara, 2010-2017 Anna M. Bigatti

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


#include "CoCoA/FreeModule.H"
#include "CoCoA/GBEnv.H"
#include "CoCoA/PPWithMask.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/SugarDegree.H"
#include "CoCoA/TmpGTypes.H"
#include "CoCoA/utils.H"

#include<list>
// using std::list; // for GPolyList; GPolyPtrList;
#include<vector>
// using std::vector; // for GPolyPtrVector;


namespace CoCoA
{

  class GPair;     // fwd decl
  class Reductors; // fwd decl (here)


enum ClearMarker {clear};

class GPoly;


typedef std::list<GPoly> GPolyList;


// This class is only a wrapper around the polys.
// Some other data useful for GBasis computation is added.

/*-----------------------------------------------------------------*/

  class GPoly
  {
  public:
    enum ReductionFlag { Full, OnlyTail };

  public:
    GPoly(ConstRefRingElem,
          const GRingInfo&,
          long age=0);
    GPoly(ConstRefRingElem,
          ConstRefPPMonoidElem theLPP,//          const PPWithMask&,  ANNA
          ConstRefRingElem,
          const GRingInfo&,
          long age=0);///< this ctor takes LPP and LC from input for Dynamic Algorithm
    GPoly(RingElem&,
          const GRingInfo&,
          const ClearMarker,
          long age=0);  ///< this ctor resets RingElem&
    GPoly(RingElem&,
          ConstRefPPMonoidElem theLPP,//          const PPWithMask&,  ANNA
          ConstRefRingElem,
          const GRingInfo&,
          const ClearMarker,
          long age=0);///< this ctor resets RingElem& and takes LPP and LC from input for Dynamic Algorithm
    explicit  GPoly(const GRingInfo&);
    GPoly(const GPoly&);
    ~GPoly(){};
    inline bool operator<(const GPoly& f) const {return myAge<f.myAge;};
    GPoly& operator=(const GPoly& rhs);
    void AssignClear(GPoly&);
    bool operator==(const GPoly&)const;
    void myAppendClear(RingElem&);  ///< Dangerous
    void myAppendClear(GPoly&); ///< Dangerous
    const GRingInfo& myGRingInfo()const{return myGRingInfoValue;};
    CoeffEncoding::type myCoeffRingType()const{return myGRingInfoValue.myCoeffRingType();};
    bool operator!=(const GPoly&)const;
    inline void Deactivate(){IamActive=false;};
    inline void mySetMinimalGen() {myMinimalGenLevel=0;};
    void myInitializeSugar(const SugarDegree& s);
    void myAssignSPoly(const GPair&, const long age);
    ConstRefPPMonoidElem myLPPForDiv() const { return PP(myLPPForDivwMask); }
    friend long component(const GPoly& g) { return g.myComponent; }
    void MultiplyByPP(ConstRefPPMonoidElem);
    friend bool IsZero(const GPoly&);
    friend bool IsMinimalGen(const GPoly& g) {return g.myMinimalGenLevel==0;}
    friend std::ostream& operator<<(std::ostream&, const GPoly&);

    inline friend ConstRefPPMonoidElem LPPForOrd(const GPoly& f) {return f.myLPPForOrd;}
    inline friend ConstRefRingElem LC(const GPoly& f) {return f.myLCValue;}
    inline friend const PPWithMask& LPPForDivwMask(const GPoly& f) {return f.myLPPForDivwMask;}

    inline friend bool IsActive(const GPoly& f) {return f.IamActive;}
    inline friend long NumTerms(const GPoly& gp) {return gp.myNumTerms;}
    inline friend const degree& wdeg(const GPoly& gp) {return gp.myWDeg;}
    inline friend const SugarDegree& sugar(const GPoly& gp) {return gp.mySugar;}
    inline friend ConstRefRingElem poly(const GPoly& gp) {return gp.myPolyValue;}
    inline friend long age(const GPoly& gp) {return gp.myAge;}
    friend long MakeKey(const GPoly&);

    friend const SparsePolyRing& owner(const GPoly& gp);
    friend const PPMonoid& PPM(const GPoly& gp);
    friend const ring& CoeffRing(const GPoly& gp);

    /**@name Pseudo-saturation functions */
    //@{
    // real body should be in SparsePolyRing
    void smart_dehomog_DRL(long DHVar);
    void smart_dehomog(long DHVar);
    friend long max_common_wdeg(GPoly&,long Var);
    friend void GPolyList2PolyListClear(GPolyList&,PolyList&);
    friend void GPolyList2PolyList(const GPolyList&,PolyList&);

    //@}

    /**@name Public reduction functions */
    //@{
    //  void reduce(const std::vector<GPoly*>& G); // is this useful??
    void myReduce(const Reductors& theReductors);
    void myReduceTail(const GPoly& g);
    void myReduceTail(const Reductors& theReductors);
    //@}

    /**@name member fields of GPoly  */
    //@{
  private:
    friend class GReductor;

    PPWithMask myLPPForDivwMask; ///< precomputed LPP of the polynomial, stored with DivMask, for divisibility tests
    PPMonoidElem myLPPForOrd; ///< precomputed LPP of the polynomial, for sorting
    RingElem myLCValue;// use iterator, can be big
    RingElem myPolyValue; ///< the polynomial proper
    //const SparsePolyRing& myPolyRing; ///< owner(myPoly) as SparsePolyRing
    const GRingInfo& myGRingInfoValue;
    bool IamActive; ///< mutable, the rest const
    //    bool IamMinimalGen; ///< mark for minimal generators
    long myMinimalGenLevel; ///< mark for minimal generators
    long myNumTerms; ///< precomputed length of the polynomial: for geobuckets
    long myComponent; /// <the module component of the LPP
    degree myWDeg; ///< precomputed degree of the polynomial
    SugarDegree mySugar; ///< "sugar" of the polynomial: a sort of degree
    long myAge; ///< age of the polynomial: for ordering reducers
    //@}

    //----------  GPoly private member functions  ----------
    void myPolySetSPoly(const GPoly& f, const GPoly& g); ///< set only s-polynomial in myPoly
    inline ConstRefRingElem myPoly() const {return myPolyValue;}; // DANGER DANGER DANGER: the LPP, LC, wdeg, StdDeg, sugar, wsugar of this polynomial DIFFER form the corresponding fields for the
                                                                   // GPoly if we are using the Dynamic Algorithm. NEVER use it if you can avoid it

    void myUpdateLenLPPLCDegComp();
    friend void monic(std::list<GPoly>&);
    //    friend void ChooseReduction(char r='g');   ///< 'g' = geobucket; 'b' = basic;

    /**@name Private functions defined in reduce.C*/
    //@{
    friend GPoly* FindReducer(const GPoly&, const Reductors& theReductors);
    //@}

    // these functions need to access myPoly(). See if they can be pared down/put here
    // they are located in GReductor.C
    friend void DeEmbedPoly(ModuleElem&,const GPoly&,long);
    friend void IntEmbedGPolyList(GPolyList&, GPolyList&);
    friend void DeEmbedPoly(RingElem&, const GPoly&);
  };//class GPoly

  // These are only for internal use
  typedef std::list<GPoly*> GPolyPtrList;
  typedef std::vector<GPoly*> GPolyPtrVector;


  GPoly* FindReducer(const PPWithMask& pm,
                     const long component,
                     const Reductors& theReductors);
  GPoly* FindReducer(const GPoly&, const Reductors& theReductors);

  //  std::vector<GPoly*>::const_iterator FindReducer(const geobucket&, const std::vector<GPoly*>& G);


   void monic(PolyList&);
   void monic(GPolyList&);
   void power(RingElem& theResult,
               const std::vector<RingElem>& theV,
               const degree& the_d);

  // -------------------- inline functions -------------------- //

  inline ConstRefPPMonoidElem LPPForDiv(const GPoly& f)
  { return f.myLPPForDiv(); }

  inline const SparsePolyRing& owner(const GPoly& gp)
  { return gp.myGRingInfoValue.myNewSPR(); }


  inline const PPMonoid& PPM(const GPoly& gp)
  { return PPM(owner(gp)); }


  inline const ring& CoeffRing(const GPoly& gp)
  { return CoeffRing(owner(gp)); }


  // -------------------- ReductorData -------------------- //
  class ReductorData

  {
  public:
    ReductorData(GPoly* p, long p_component, long count=0);
    ReductorData(const ReductorData&);
    ~ReductorData(){};
    bool operator<(const ReductorData& RD) const {return myKey<RD.myKey;};
    friend std::ostream& operator<<(std::ostream&, const ReductorData&);
    //inline friend long Deg(const ReductorData& RD) {return Deg(*myGPoly);};
    friend const SparsePolyRing& owner(const Reductors&);
    friend const PPMonoid& PPM(const Reductors&);
    friend long component(const ReductorData& RD) { return RD.myComponent; }
    bool IamNotToBeUsed() const { return myIamNotToBeUsedValue; }
    void mySetIamNotToBeUsed(const bool b) { myIamNotToBeUsedValue=b; }
    GPoly* myGetGPolyPtr() const { return myGPolyPtr; }



  private:
    friend class Reductors;
    GPoly* myGPolyPtr;
    long myKey; ///< integer "key" used to sort the list of reducers
    long myComponent;
    mutable PPWithMask myLPPForDivwMask;
    mutable long myCount; ///< how many times it has been used
    mutable bool IamBorelUpdated; ///< true if its degree is the current degree (Homogeneous algorithm)
    bool myIamNotToBeUsedValue;///<true if not to be used for reducing (used in final interreduction)

    friend GPoly* FindReducer(const PPWithMask& pp, const long pp_component,
                              const Reductors& theReductors);
    bool myBorelUpdate(ConstRefPPMonoidElem pp, const Reductors& theReductors);
  };


  class Reductors{

  public:
    enum UseBorelFlag { UseBorel, DontUseBorel };

    Reductors(const GRingInfo& P);
    Reductors(const GRingInfo& P, const UseBorelFlag);
    ~Reductors(){};
    void Insert(GPoly*, const long count=0);
    UseBorelFlag IhaveBorelReductors()const
    {if(IhaveBorelReductorsFlag) return UseBorel; else return DontUseBorel;};

    long myLen() const {return len(myReductors);};
    void myStampaReductors(std::ostream&) const;
    const GRingInfo& myGRingInfo()const{return myGRingInfoValue;};
    //    GPoly* LastGPolyPtr(){return myReductors.back().myGPolyPtr;};//copia ptr  // never used
    //    PPMonoidElem LastGPolyLPP(){return myReductors.back().myLPP;};//copia PP  // never used
    void myReduceTails(const GPoly&);
    void OrderedInterreduce(const GPoly&);
    void SuperInterreduce(const GPoly&);

    friend const SparsePolyRing& owner(const Reductors&);
    friend const PPMonoid& PPM(const Reductors&);
    friend void interreduce(Reductors&, GPoly);
    friend void OrderedInterreduce(Reductors&, GPoly);
    friend void reduce(const Reductors&);
    friend void ReduceTail(const Reductors&);
    friend GPoly* FindReducer(const PPWithMask& pp,
                              const long pp_component,
                              const Reductors& theReductors);
    // HERE TO BE NOTICED friend const SparsePolyRing& owner(const Reductors&);
    // HERE TO BE NOTICED friend const PPMonoid& PPM(const Reductors&);
			
    // This function takes params that are part of GReductor class.
    // In the future, a new interface class will be built that have this funct. as a member
    // This new class would insulate Reductors from end user, which will see only
    // this new interface class which would rule all reductors functions.
    // data in a GReductor will be pruned accordingly.
    void myBorelReductorsUpdateInNextDegree();
    void myClear();
    std::vector<ReductorData>::iterator find(GPoly*);// Find the (unique) Reductor data which ptr is equal to GPolyPtr. Returns end() if nothing is found.

    //    friend void swap(Reductors&,Reductors&); // ANNA: no longer needed?

    //@name member fields for Reductors
    //@{
  private:
    /// the poly used for reductions NB polys are ptrs to myPolys
    std::vector<ReductorData> myReductors;
    mutable std::vector<ReductorData> myBorelReductors;
    GPolyList myBorelGPolys;
    bool IhaveBorelReductorsFlag;
    const GRingInfo& myGRingInfoValue;
    //@}

  };// End class Reductors


  SparsePolyRing owner(const PolyList& thePL);
  SparsePolyRing owner(const std::vector<RingElem>& thePV);
  SparsePolyRing owner(const GPolyList& theGPL);

  void PolyList2PolyVectorClear(PolyList&,std::vector<RingElem>&);
  void PolyVector2PolyListClear(std::vector<RingElem>&,PolyList&);

}// end namespace cocoa



// RCS header/log on the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/TmpGPoly.H,v 1.28 2024/09/04 15:25:31 bigatti Exp $
// $Log: TmpGPoly.H,v $
// Revision 1.28  2024/09/04 15:25:31  bigatti
// Summary: just spacing
//
// Revision 1.27  2024/09/04 15:18:04  bigatti
// Summary: renamed Component into component
//
// Revision 1.26  2024/05/30 15:07:03  bigatti
// Summary: renamed interreduce into myReduceTails
//
// Revision 1.25  2022/02/18 14:11:50  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.24  2017/05/02 12:10:44  bigatti
// -- commented dould lines for friend LC and LPPForDivwMask
//
// Revision 1.23  2017/04/26 12:53:12  bigatti
// -- some cleaning: len --> NumTerms, f.IsActive() --> IsActive(f), ..
//
// Revision 1.22  2017/04/18 09:33:32  bigatti
// -- "age" now lower-case
//
// Revision 1.21  2013/10/28 13:15:35  bigatti
// -- IamMinimalGen --> myMinimalGenLevel
//
// Revision 1.20  2013/06/12 08:50:03  bigatti
// -- added IamMinimalGen
//
// Revision 1.19  2012/10/16 09:51:53  abbott
// Replaced  RefRingElem  by  RingElem&
//
// Revision 1.18  2011/12/05 16:54:40  bigatti
// -- changed: MachineInteger --> MachineInt (just in comment)
//
// Revision 1.17  2011/03/11 17:39:54  bigatti
// -- changed  unsigned int --> long
//
// Revision 1.16  2011/03/11 15:32:31  bigatti
// -- changed size_t into long
//
// Revision 1.15  2011/03/11 14:50:39  abbott
// Changed some size_t into long.
//
// Revision 1.14  2010/07/16 09:29:58  bigatti
// -- minor cleaning and coding conventions
//
// Revision 1.13  2010/05/14 09:53:09  bigatti
// -- removed empty ctor for SugarDegree
// -- added marker for SugarDegree(uninitialized)
// -- SugarDegree for GBasis input is initialized by myPrepareGBasis
//
// Revision 1.12  2010/03/23 14:43:07  bigatti
// -- class GRingInfo estracted from TmpGPoly
//
// Revision 1.11  2009/10/27 17:15:14  bigatti
// -- fixed: using sugar(g)->myWSugar() insted of wsugar(g)
//
// Revision 1.10  2009/04/27 13:12:44  bigatti
// -- added some temporary functions (DegLess...)
//
// Revision 1.9  2009/04/27 12:23:36  bigatti
// -- added SaturatingAlgNoDRL  to enum
//
// Revision 1.8  2008/09/19 13:33:42  bigatti
// -- added: Sat algorithm (M.Caboara)
//
// Revision 1.7  2008/09/16 15:03:42  bigatti
// -- added LPPForDiv
// -- changed LPP into LPPForOrd
//
// Revision 1.6  2008/07/09 16:08:04  abbott
// Cosmetic tidying.
//
// Revision 1.5  2008/04/21 12:32:54  abbott
// Corrected size_t into std::size_t in several header files; in some cases,
// replaced size_t with MachineInt (with consequent changes to impl files).
//
// Revision 1.4  2007/12/05 11:06:24  bigatti
// -- changed "size_t StdDeg/myStdDeg(f)" into "long"  (and related functions)
// -- changed "log/myLog(f, i)" into "MaxExponent/myMaxExponent(f, i)"
// -- fixed bug in "IsOne(ideal)" in SparsePolyRing.C
//
// Revision 1.3  2007/12/04 14:27:07  bigatti
// -- changed "log(pp, i)" into "exponent(pp, i)"
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
// Revision 1.13  2007/03/08 18:42:06  cocoa
// Cleaned up whitespace.
//
// Revision 1.12  2007/03/07 22:24:03  bigatti
// -- reintroduced TmpGTypes.H (waiting for a better solution)
//
// Revision 1.11  2007/03/07 17:04:31  cocoa
// -- several changes by M.Caboara: more operations on ideals,
//    exception cleaner, coding conventions, WSugar, dynamic
//
// Revision 1.10  2006/12/22 11:37:40  cocoa
// Removed operator<< for PolyList and GPolyList.  A few consequential changes.
//
// Revision 1.9  2006/12/06 17:11:16  cocoa
// -- removed #include "config.H"
//
// Revision 1.8  2006/12/04 13:55:54  cocoa
// -- added: sugar for GradingDim > 0  (called wsugar)
//
// Revision 1.7  2006/11/24 17:17:58  cocoa
// -- reorganized includes of header files
//
// Revision 1.6  2006/11/20 14:57:17  cocoa
// -- added: (standard) sugar for modules
// -- fixed: non-homogeneous sysygies
// -- minor fixes     [M.Caboara]
//
// Revision 1.5  2006/11/02 13:25:44  cocoa
// Simplification of header files: the OpenMath classes have been renamed.
// Many minor consequential changes.
//
// Revision 1.4  2006/10/06 16:32:06  cocoa
// -- changed: GPoly::SPoly --> GPoly::myAssignSPoly
// -- changed: Len(const GPoly&) --> len(const GPoly&)
// -- added:   poly(const GPoly&)
// -- added:   GPoly::myUpdateLenLPPDegComp()
// -- in reduce.C added functions for computing sugar during reduction
//
// Revision 1.3  2006/10/06 14:04:16  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.2  2006/08/17 09:28:56  cocoa
// -- added: sugar
// -- added: controls on homogeneity
//
// Revision 1.1.1.1  2006/05/30 11:39:36  cocoa
// Imported files

#endif
