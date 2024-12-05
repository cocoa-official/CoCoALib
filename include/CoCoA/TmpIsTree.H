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

#endif
