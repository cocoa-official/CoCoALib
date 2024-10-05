#ifndef CoCoA_IsTree_H
#define CoCoA_IsTree_H

//   Copyright (c)  2006-2010  John Abbott,  Anna M. Bigatti
//   Original author: 2006-2010  Massimo Caboara

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


#include "CoCoA/DynamicBitset.H"
#include "CoCoA/TmpGTypes.H"
#include "CoCoA/time.H"
#include "CoCoA/utils.H"


#include <cstddef>
//using std::size_t;
#include <vector>
//using std::vector;
#include <list>
//using std::list;
#include <bitset>
//using std::bitset;
#include <utility>
//using std::pair;


namespace CoCoA
{

  class PolyRing; // forward declaration -- defined in PolyRing.H
  class RingElem; // forward declaration -- defined in ring.H

  typedef DynamicBitset facet;

// WARNING TODO: facets are supposed to have all the same len.
// Add a static field and some check.



// Add the PolyRing (at least Monomial) to the FacetComplex class.
// This helps the I/O. Or, change the I/O in such a way to get the
// PR from the PL, special case list []
// functions: PolyList2FacetList, FacetComplex(const PolyRing&, const PolyList&);
// a FacetComplex is, more or less, a set of facets with lex ordering.

  typedef std::list<facet>::const_iterator FacetComplexConstIter;
  typedef std::list<facet>::iterator FacetComplexIter;
  typedef std::vector<std::pair<FacetComplexConstIter,std::vector<FacetComplexConstIter> > >::const_iterator conn_block_const_itr;
  typedef std::vector<std::pair<FacetComplexConstIter,std::vector<FacetComplexConstIter> > >::iterator ConnBlockIter;


  class FacetComplex
  {
  public:
    friend class ConnectionBlock;

    FacetComplex(const PolyRing&, const PolyList&);
    FacetComplex(const FacetComplex& theFacetComplex){myElems=theFacetComplex.myElems;};
    FacetComplex(){};
    FacetComplex(const std::list<facet>&);
    ~FacetComplex(){};
    unsigned int myNumIndets()const;
    unsigned int myGetIndex(const facet&)const;
    friend std::ostream& operator<<(std::ostream&, const FacetComplex&);
    FacetComplex& operator=(const FacetComplex&);
    void myInsert(const facet& f){myElems.push_back(f);};
    void myErase(const facet& f){myElems.remove(f);};
    unsigned int mySize()const{return len(myElems);};
    bool myAmIEmpty()const{return myElems.empty();};
    std::list<facet> myFacetComplex2FacetList()const;
    FacetComplex mySetDifference(const facet&)const;
    FacetComplex delta(const facet&,const facet&,const facet&)const;// TREE SPECIFIC the delta operation in the extended abstract
    FacetComplex delta_new(const facet&,const facet&,const facet&)const;// TREE SPECIFIC the delta operation in the paper

    bool AreConnected(const facet&,const facet&)const;// old algorithm
    bool AreConnected_new(const facet&,const facet&)const;// new O(nl) algorithm
    bool IsTripleCondition(const facet&,const facet&,const facet&)const;// TREE SPECIFIC triple condition
    std::list<facet> myIsTreeNoOpt();// TREE SPECIFIC
    std::list<facet> myIsTreeOpt();// TREE SPECIFIC
    std::list<facet> myIsTreeCBNoOpt();// TREE SPECIFIC
    std::list<facet> myIsTreeCBOpt();// TREE SPECIFIC
    void myClear(){myElems.clear();};

  private:
    // they are just used in delta_new(...)
    void myMakeXj(std::list<unsigned int>&,const unsigned int j)const;
    void myMakeG(std::vector<unsigned int>&,const std::vector<unsigned int>& P,
                 const std::list<unsigned int>& xj)const;
    std::list<facet> myElems;
  };// end class FacetComplex


  std::vector<RingElem> FacetList2PolyList(const SparsePolyRing&, const std::list<facet>&);


// Graph sparse description of the direct connectedeness of the facets of a FacetComplex
// A connection block is a std::vector of <itr,std::vector<itr>>, where the itrs are really ptrs to some facet in the FacetComplex
// the meaning: if p is in the vector associated to q, p and q are directly connected
  class ConnectionBlock
  {
  public:
    friend class FacetComplex;
    ConnectionBlock(const FacetComplex&);
    ConnectionBlock(){};
    friend std::ostream& operator<<(std::ostream&, const ConnectionBlock&);
    ConnBlockIter erase(ConnBlockIter&);
  private:
    std::vector<std::pair<FacetComplexConstIter,std::vector<FacetComplexConstIter> > > my_array;
  };// end class ConnectionBlock


///////////// inline facet functions ///////////////////////////////////////////

  // The IsFLesser operation from the paper
  inline bool IsFLesser(const facet& f,const facet& g1,const facet& g2)
  {
    facet F1(FacetIntersection(f,g1));
    facet F2(FacetIntersection(f,g2));
    return contains(F1,F2);
  }//IsFLesser

}// end namespace cocoa#endif


// RCS header/log on the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/TmpIsTree.H,v 1.12 2022/02/18 14:11:50 abbott Exp $
// $Log: TmpIsTree.H,v $
// Revision 1.12  2022/02/18 14:11:50  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.11  2022/02/14 14:09:36  bigatti
// Summary: added date for original author
//
// Revision 1.10  2022/02/11 09:49:34  abbott
// Summary: Updated copyright notices (redmine 855)
//
// Revision 1.9  2014/07/31 14:45:16  abbott
// Summary: Merged io.H and UtilsTemplate.H into new header VectorOperations.H
// Author: JAA
//
// Revision 1.8  2010/05/28 15:50:42  bigatti
// -- cleaning
// -- commented out ctor taking a list of long
// -- moved some "facet" functions to TmpIsTree
//
// Revision 1.7  2010/03/30 15:36:39  bigatti
// -- using DynamicBitset.HC (with former code for facet)
//
// Revision 1.6  2010/03/16 16:38:29  bigatti
// -- added comments for TREE SPECIFIC operations (in the way for
//    separating general square free operations)
//
// Revision 1.5  2008/09/17 14:38:07  bigatti
// -- changed WordLen
// -- removed comment
//
// Revision 1.4  2007/10/30 17:14:11  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.3  2007/09/24 14:17:23  abbott
// Changed post-increment to pre-increment several times.
// Minor changes to keep prerelease gcc-4.3 quiet.
// Minimal code cleaning.
//
// Revision 1.2  2007/03/27 15:27:05  bigatti
// -- minor update for TmpIsTree + test
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.7  2007/03/07 22:24:03  bigatti
// -- reintroduced TmpGTypes.H (waiting for a better solution)
//
// Revision 1.6  2007/03/07 17:39:24  bigatti
// -- fix for PolyList and VectorList
//
// Revision 1.5  2007/03/07 14:59:02  cocoa
// -- renamed complex --> FacetComplex
//
// Revision 1.4  2006/11/24 17:35:22  cocoa
// -- reorganized includes of header files
//
// Revision 1.3  2006/10/06 14:04:15  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.2  2006/08/07 21:23:25  cocoa
// Removed almost all publicly visible references to SmallExponent_t;
// changed to long in all PPMonoid functions and SparsePolyRing functions.
// DivMask remains to sorted out.
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.1  2006/05/16 09:03:11  cocoa
// -- first import
//

#endif
