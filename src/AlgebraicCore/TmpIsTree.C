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


#include "CoCoA/TmpIsTree.H"
#include "CoCoA/DynamicBitset.H"
#include "CoCoA/SparsePolyOps-RingElem.H"
//#include "CoCoA/SparsePolyRing.H"

#include <algorithm>
using std::find;
//#include <bitset>
using std::bitset;
#include <functional>
//??
#include <iostream>
using std::ostream;
using std::endl;
using std::flush;
//#include <list>
using std::list;
#include <utility>
using std::pair;
//#include <vector>
using std::vector;

//static const bool MAX_DEBUG = false;

namespace CoCoA
{



/*
--------------------------------------------------------------------------------
------ FacetComplex class functions -----------------------------------------------
--------------------------------------------------------------------------------
*/
  FacetComplex::FacetComplex(const PolyRing& /*P*/, const PolyList& PL)
   {
     for (const RingElem& f :PL)
      {
        myInsert(facet(LPP(f)));
      }
    }//FacetComplex


    FacetComplex::FacetComplex(const list<facet>& L)
    {
      for (const facet& f: L)
        myInsert(f);
    }//FacetComplex


  std::ostream& operator<<(std::ostream& out, const FacetComplex& the_FacetComplex)
  {
    if (!out) return out;  // short-cut for bad ostreams
    out<<"Length is "<<the_FacetComplex.mySize()<<endl;
    for (const facet& elem: the_FacetComplex.myElems)
    {
       out << elem << endl;
    }
    return out;
  }//operator<<


  unsigned int FacetComplex::myNumIndets()const
  {
    if (myAmIEmpty())
      return 0;
    else
      return len(myElems.front());
    // Anna: was return myElems.front().myVecSize()*DynamicBitset::ourNumBitsInBlock;

  }//FacetComplex::myNumIndets


  // f1=f2
  FacetComplex& FacetComplex::operator=(const FacetComplex& the_FacetComplex)
  {
    myElems=the_FacetComplex.myElems;
    return *this;
  }

  list<facet> FacetComplex::myFacetComplex2FacetList() const
  {
    list<facet> L;
    for (const facet& elem: myElems)
    {
	 L.push_back(elem);
     }
     return L;
  } //FacetComplex::myFacetComplex2FacetList


   FacetComplex FacetComplex::mySetDifference(const facet& f) const
   {
     FacetComplex c;
     facet f_tmp;
     for (const facet& elem: myElems)
     {
       f_tmp = elem - f;// possible OPT: operation and check at the same time
       if (f_tmp.IamNotEmpty()) {c.myInsert(f_tmp);}
     }
     return c;
   } //FacetComplex::mySetDifference

   // The delta operation from the extended abstract
   FacetComplex FacetComplex::delta(const facet& f,const facet& g1,const facet& g2) const
   {
      FacetComplex L;
      facet f_tmp;
      for (const facet& elem: myElems)
      {
        if (  (elem != g1 && (IsFLesser(f,g1,elem)))//OPT: this three ops can be combined in one?
	                 ||
	      (elem != g2 && (IsFLesser(f,g2,elem)))
	    )
	 {
           f_tmp = elem-f;
           if (f_tmp.IamNotEmpty())
           {
             L.myInsert(f_tmp);
           }
          }// Big check if
      }// For
      return L;
   }//FacetComplex::delta

   // The delta operation from the paper
  FacetComplex FacetComplex::delta_new(const facet& f,
                             const facet& g1,
			                       const facet& g2)const
   {
      FacetComplex L;
      facet f_tmp;
      facet g1_int_g2=FacetIntersection(g1,g2);
      for (const facet& elem: myElems)
      {
        if (FacetIntersection(f,elem) == g1_int_g2)
	 {
           f_tmp = elem-f;
           if (f_tmp.IamNotEmpty())
             L.myInsert(f_tmp);
         }// Big check if
      }// For
      L.myInsert(g1-f);
      L.myInsert(g2-f);
      return L;
   }//delta_new


   // Are g1 and g2 connected in the FacetComplex?
   bool FacetComplex::AreConnected(const facet& g1,const facet& g2)const
   {
     if (g1.IamEmpty() || g2.IamEmpty()) return false;
     if (g1==g2 || AreDirectlyConnected(g1,g2)) return true;
     FacetComplex d(*this);
     FacetComplexIter it;
     facet g1_tmp=g1;
     d.myErase(g1);
     d.myErase(g2);
     bool DoneSomething;
     do {
       DoneSomething = false;
       it = d.myElems.begin();
       while (it!=d.myElems.end())
       {
         if (AreDirectlyConnected(g1_tmp,*it))
	 {
	   g1_tmp+=*it;
	   it=d.myElems.erase(it);
	   if (AreDirectlyConnected(g1_tmp,g2))
	     return true;
	   DoneSomething=true;
	 }
	 else
	   ++it;
       }// end while do
     } while (DoneSomething);
    return false;
   }// AreConnected


  unsigned int FacetComplex::myGetIndex(const facet& f) const
  {
    unsigned int i=0;
    for (const facet& elem: myElems)
    {
      if (elem == f)
         return i;
      ++i;
    }
    return 0;
  }


 // For use in re_connected_new
  void FacetComplex::myMakeXj(list<unsigned int>& theL,const unsigned int j)const
 {
   theL.clear();
   unsigned int i=0;
   for (const facet& elem: myElems)
   {
     if (elem.myIsEntryThere(j))
       theL.push_back(i);
     ++i;
   }
   return;
 } //FacetComplex::MakeXj



 // For use in re_connected_new
  void FacetComplex::myMakeG(vector<unsigned int>& theL,
                        const vector<unsigned int>& P,
                        const list<unsigned>& xj)const
 {

   theL.clear();
   bool found_x=false;
   if (xj.empty())
   {
     for (unsigned int i=0; i!=mySize(); ++i)
       theL.push_back(i);
     return;
   }//if
   const unsigned int x=xj.front();
   for (unsigned int i: P)
   {
     found_x=false;
     for (unsigned int k: xj)
     {
       if (P[k] == i)
       { found_x=true; break; }
     }//for
     if (found_x) theL.push_back(x); else theL.push_back(i);
   }// for
   return;
 }//FacetComplex::MakeG


   // OPT: don't create a new FacetComplex, use a non modifying alg on the old one
   bool FacetComplex::AreConnected_new(const facet& g1,const facet& g2)const
   {
     if (g1.IamEmpty() || g2.IamEmpty()) return false;
     if (g1==g2 || AreDirectlyConnected(g1,g2)) return true;
     if (this->mySize()==2) return false;// they are 2, they are not directly connected, they are not connected
     unsigned int g1_index=myGetIndex(g1);
     unsigned int g2_index=myGetIndex(g2);
     vector<unsigned int> P,P1;
     for (unsigned int i=0;i!=mySize();++i)
       P.push_back(i);
     list<unsigned int> xj;
     vector<unsigned int> g;
     for (unsigned int j=0; j!=myNumIndets(); ++j)
     {
       myMakeXj(xj,j);
       if (len(xj)>1)
       {
         myMakeG(g,P,xj);
         P1.clear();
         for (unsigned int i: P)
           P1.push_back(g[P[i]]);
        P1.swap(P);
       }//if
     }//for
     return (P[g1_index]==P[g2_index]);
    }// AreConnected_new


   bool FacetComplex::IsTripleCondition(const facet& f,const facet& g1,const facet& g2)const
   {
     if (IsFLesser(f,g1,g2) || IsFLesser(f,g2,g1)) return false;
     return (this->delta(f,g1,g2)).AreConnected(g1-f,g2-f);
   }


   // The paper algorithm without any optimization
  list<facet> FacetComplex::myIsTreeNoOpt()
  { //long triples_checked=0;long conn_check_no=0;//stats
    //    unsigned int row_being_worked_on=0;//stats
    //    std::cout<<endl<<"Standard Algorithm without optimization"<<endl;
    FacetComplexConstIter it1stop=--myElems.end();
    FacetComplexConstIter it2start;
    list<facet> cycle; // void if this is a tree, contains a 3-cycle otherwise
    for (FacetComplexConstIter it=myElems.begin();
         it!=myElems.end();++it)
    {// std::cout<<"working with row "<<row_being_worked_on<<endl;
      //      ++row_being_worked_on;
      for (FacetComplexConstIter it1=myElems.begin();
          it1!=it1stop;++it1)
      {  if (it!=it1)
         { it2start=it1;++it2start;
           for (FacetComplexConstIter it2=it2start;it2!=myElems.end();++it2)
	         { if (it!=it2)
                   { //++triples_checked; // stats
	            if (!(IsFLesser(*it,*it1,*it2) || IsFLesser(*it,*it2,*it1)))
	            {  // ++conn_check_no; // stats
		             if   (this->delta_new(*it,*it1,*it2).AreConnected(*it1-*it,*it2-*it))
	              { cycle.push_back(*it);
	                cycle.push_back(*it1);
	                cycle.push_back(*it2);
                        //	                std::cout<<endl<<"Standard Algorithm without optimization"<<endl;
                        //                  std::cout <<"triples_checked= "<<triples_checked<<endl;
                        //	                std::cout <<"conn_check_no  = "<<conn_check_no<<endl;
	                return cycle;
		            }// if (this->
	            }// if if_f_lesser
	          }//if it!=it2
	        }// third for
        }// if it!=it1
      }// second for
     // No opt here if (is_f_useless){myElems.erase(it);}
   }// first for
    //   std::cout<<endl<<"Standard Algorithm without optimization"<<endl;
    //   std::cout <<"triples_checked= "<<triples_checked<<endl;
    //   std::cout <<"conn_check_no  = "<<conn_check_no<<endl;
   return cycle;
   } //myIsTreeNoOpt


   // The paper algorithm with the useless facet optimization
  list<facet> FacetComplex::myIsTreeOpt()
  {
    // long triples_checked=0;long conn_check_no=0;//stats
    //    unsigned int row_being_worked_on=0;//stats
    bool is_f_useless=true;
    //    std::cout<<endl<<"Standard Algorithm with optimization"<<endl;
    FacetComplexConstIter it1stop=--myElems.end();
    FacetComplexConstIter it2start;
    list<facet> cycle; // void if this is a tree, contains a 3-cycle otherwise
    FacetComplexIter it=myElems.begin();
    while (it!=myElems.end())
    { is_f_useless=true;
      //      std::cout<<"working with row "<<row_being_worked_on<<endl;
      //      ++row_being_worked_on;
      for (FacetComplexConstIter it1=myElems.begin();
          it1!=it1stop;++it1)
      {  if (it!=it1)
         { it2start=it1;++it2start;
           for (FacetComplexConstIter it2=it2start;
                it2!=myElems.end();++it2)
	  { if (it!=it2)
	    {
              // ++triples_checked; // stats
	      if (!(IsFLesser(*it,*it1,*it2) ||
	            IsFLesser(*it,*it2,*it1)))
	      { is_f_useless=false;
		// ++conn_check_no; // stats
		if   (this->delta_new(*it,*it1,*it2).AreConnected(*it1-*it,*it2-*it))
	        { cycle.push_back(*it);
	          cycle.push_back(*it1);
	          cycle.push_back(*it2);
                  //                  std::cout<<endl<<"Standard Algorithm with optimization"<<endl;
                  //                  std::cout <<"triples_checked= "<<triples_checked<<endl;
                  //                  std::cout <<"conn_check_no  = "<<conn_check_no<<endl;
	          return cycle;
		}// if (this->
	      }// if if_f_lesser
	    }//if it!=it2
	  }// third for
        }// ifit!=it1
      }// second for
     if (is_f_useless)
     {//  std::cout<<"@";
        it=myElems.erase(it);}
     else
       ++it;
   }// first for
    //   std::cout<<endl<<"Standard Algorithm with optimization"<<endl;
    //   std::cout <<"triples_checked= "<<triples_checked<<endl;
    //   std::cout <<"conn_check_no  = "<<conn_check_no<<endl;
    return cycle;
  }//myIsTreeOpt


  // This creates a connection block and walks all the possible triples there
 list<facet> FacetComplex::myIsTreeCBNoOpt()
 {  double T;
   // long triples_checked=0;long conn_check_no=0;//stats
   //    unsigned int row_being_worked_on=0;
    //    std::cout<<endl<<"Incidence matrix Algorithm without optimization"<<endl;
    T=CpuTime();
    ConnectionBlock cb(*this);
    T = CpuTime() - T;
    //    std::cout<<endl<< "Connection Block time expended is " << T << endl;
    list<facet> cycle; // void if this is a tree, contains a 3-cycle otherwise
    ConnBlockIter it=cb.my_array.begin();
    while (it!=cb.my_array.end())
    {// std::cout<<"working with row "<<row_being_worked_on<<endl;
      //      ++row_being_worked_on;
      if (len(it->second) > 1)
      { for (vector<FacetComplexConstIter>::const_iterator it1=(*it).second.begin();
	    it1!=(*it).second.end();++it1)
	{ for (vector<FacetComplexConstIter>::const_iterator it2=it1;
	      it2!=(*it).second.end();++it2)
	  { if (it1!=it2)
	    { // ++triples_checked; // stats
	      if (!(IsFLesser(*(*it).first,*(*it1),*(*it2)) ||
	            IsFLesser(*(*it).first,*(*it2),*(*it1))))
	      { // ++conn_check_no; // stats
		if (this->delta_new(*(*it).first,*(*it1),*(*it2)).AreConnected(*(*it1)-*(*it).first,*(*it2)-*(*it).first))
	        { cycle.push_back(*(*it).first);
	          cycle.push_back(*(*it1));
	          cycle.push_back(*(*it2));
                  //	          std::cout<<endl<<"Incidence matrix Algorithm without optimization"<<endl;
                  //                  std::cout <<"triples_checked= "<<triples_checked<<endl;
                  //	          std::cout <<"conn_check_no  = "<<conn_check_no<<endl;
	          return cycle;
		}// if (this->
	      }// if if_f_lesser
	    }// if it1!=it2
	  }// third for : it2
	}// second for: it1
      }// if len(it->second)>1)
      ++it;
    }// main while loop : it
    //    std::cout<<endl<<"Incidence matrix Algorithm without optimization"<<endl;
    //    std::cout <<"triples_checked= "<<triples_checked<<endl;
    //    std::cout <<"conn_check_no  = "<<conn_check_no<<endl;
    return cycle;
  }//myIsTreeCBNoOpt

   // This creates a connection block and walks all the possible triplets there
  list<facet> FacetComplex::myIsTreeCBOpt()
  { double T;
    // long triples_checked=0;long conn_check_no=0;//stats
    //    unsigned int row_being_worked_on=0;//stats
    bool is_f_useless=true;
    //    std::cout<<endl<<"Connection Block Algorithm with optimization"<<endl;
    T=CpuTime();
    ConnectionBlock cb(*this);
    T = CpuTime() - T;
    //    std::cout<<endl<< "Connection Block time expended is " << T << endl;
    list<facet> cycle; // void if this is a tree, contains a 3-cycle otherwise
    ConnBlockIter it=cb.my_array.begin();
    while (it!=cb.my_array.end())
    {// std::cout<<"working with row "<<row_being_worked_on<<endl;
      //      ++row_being_worked_on;
      is_f_useless=true;
      if (len(it->second) > 1)
      {  for(vector<FacetComplexConstIter>::const_iterator it1=(*it).second.begin();
	    it1!=(*it).second.end();++it1)
	{ for(vector<FacetComplexConstIter>::const_iterator it2=it1;
	      it2!=(*it).second.end();++it2)
	  { if (it1!=it2)
	    { // ++triples_checked; // stats
	      if (!(IsFLesser(*(*it).first,*(*it1),*(*it2))||IsFLesser(*(*it).first,*(*it2),*(*it1))))
	      { is_f_useless=false;
		// ++conn_check_no; // stats
		if (this->delta_new(*(*it).first,*(*it1),*(*it2)).AreConnected(*(*it1)-*(*it).first,*(*it2)-*(*it).first))
	        { cycle.push_back(*(*it).first);
	          cycle.push_back(*(*it1));
	          cycle.push_back(*(*it2));
                  //	          std::cout<<endl<<"Incidence matrix Algorithm with optimization"<<endl;
                  //                  std::cout <<"triples_checked= "<<triples_checked<<endl;
                  //	          std::cout <<"conn_check_no  = "<<conn_check_no<<endl;
	          return cycle;
		}// if (this->
	      }// if if_f_lesser
	    }// if it1!=it2
	  }// third for : it2
	}// second for: it1
      }// if len(it->second) > 1
      if (is_f_useless)
      {// std::cout<<"@";
	it=cb.erase(it);// WARNING: after this operation the cb is correct only from the current row downward, excluded
      }
      else
      {++it;}
    }// main while loop : it
    //    std::cout<<endl<<"Incidence matrix Algorithm with optimization"<<endl;
    //    std::cout <<"triples_checked= "<<triples_checked<<endl;
    //    std::cout <<"conn_check_no  = "<<conn_check_no<<endl;
    return cycle;
  }//myIsTreeCBOpt


/* end class FacetComplex */


/*
--------------------------------------------------------------------------------
------ ConnectionBlock class functions -----------------------------------------------
--------------------------------------------------------------------------------
*/

  ConnectionBlock::ConnectionBlock(const FacetComplex& the_FacetComplex)
  {
     my_array.reserve(the_FacetComplex.mySize());

     for (FacetComplexConstIter it1=the_FacetComplex.myElems.begin();
          it1!=the_FacetComplex.myElems.end();++it1)
     {
       vector<FacetComplexConstIter> v;
       for (FacetComplexConstIter it2=the_FacetComplex.myElems.begin();
          it2!=the_FacetComplex.myElems.end();++it2)
       {
         if (it1!=it2 && AreDirectlyConnected(*it1,*it2))
	       {
	          v.push_back(it2);
	       }
       }// end nested for
       my_array.push_back(make_pair(it1,v));
     }// end first for
  }//ConnectionBlock

 std::ostream& operator<<(std::ostream& out, const ConnectionBlock& the_ConnectionBlock)
  {
    if (!out) return out;  // short-cut for bad ostreams
    out<<"[";
    for (conn_block_const_itr it=the_ConnectionBlock.my_array.begin();
         it!=the_ConnectionBlock.my_array.end();++it)
    {
       out<<*((*it).first)<<"| ";
       for (const auto& FC: it->second)
         out<<" "<<(*FC)<<" ";
       out<<"*";
    }//for
    out<<"]";
    return out;
  }//operator<<


// This version is safe, erases all the occurrences
    ConnBlockIter ConnectionBlock::erase(ConnBlockIter& row_to_erase)
    {
      vector<FacetComplexConstIter>::iterator to_be_killed;
      FacetComplexConstIter elem_to_del=((*row_to_erase).first);
//      for (ConnBlockIter it=my_array.begin();it!=my_array.end();++it)
      for (auto& iter: my_array)
      {
        to_be_killed=find(iter.second.begin(),iter.second.end(),elem_to_del);
	      if (to_be_killed!=iter.second.end())
	        iter.second.erase(to_be_killed);
      }
      return my_array.erase(row_to_erase);
    }//ConnectionBlock::erase


  //----------------------------------------------------------------------
  // functions of class facet

  // Not efficient, efficiency not needed at the moment
  RingElem Facet2RingElem(const SparsePolyRing& theP,const DynamicBitset& b)
  {
    return monomial(theP, 1, NewPP(PPM(theP), b));
  }//Facet2RingElem


  // Not efficient, efficiency not needed at the moment
  std::vector<RingElem> FacetList2PolyList(const SparsePolyRing& theP,const std::list<DynamicBitset>& theFL)
  {
    vector<RingElem> PL;
    for (const DynamicBitset& bitset: theFL)
      PL.push_back(Facet2RingElem(theP,bitset));
    return PL;
  }//FacetList2PolyList

 // end functions of class facet


}// end namespace cocoa





/*

Some future optimization:

for IsFLesser: proceed in the test word by word:
compute the first word of g1-f, the first word of g2-f, the check

ConnectionBlock: use ptr and not iterators.

sparse representation for facets? Use vector, (sort?).
Only reasonable if density is much much lower than #VARS!

*/

// RCS header/log on the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/TmpIsTree.C,v 1.24 2024/07/03 15:38:54 bigatti Exp $
// $Log: TmpIsTree.C,v $
// Revision 1.24  2024/07/03 15:38:54  bigatti
// Summary: commented out unused variable  row_being_worked_on
//
// Revision 1.23  2024/03/03 21:15:01  abbott
// Summary: Commented out unused vars (for stats)
//
// Revision 1.22  2023/03/13 21:06:40  abbott
// Summary: New iterator "for" syntax; some cruft removed where the "for" has already been changed.  Redmine 1346
//
// Revision 1.21  2022/02/18 14:12:00  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.20  2022/02/14 14:32:48  bigatti
// Summary: added date for original author
//
// Revision 1.19  2022/02/11 09:49:35  abbott
// Summary: Updated copyright notices (redmine 855)
//
// Revision 1.18  2021/07/19 11:24:01  abbott
// Summary: Updated "for" loops (redmine 1346)
//
// Revision 1.17  2020/02/18 11:32:16  abbott
// Summary: redmine 1346: new for loop syntax
//
// Revision 1.16  2018/05/18 16:38:51  bigatti
// -- added include SparsePolyOps-RingElem.H
//
// Revision 1.15  2016/11/11 14:15:34  abbott
// Summary: Added short-cut to operator<< when ostream is in bad state
//
// Revision 1.14  2016/09/08 14:14:30  bigatti
// -- commented out unused variable
//
// Revision 1.13  2014/06/17 10:15:28  abbott
// Summary: Commented out unused param in FacetComplex
// Author: JAA
//
// Revision 1.12  2014/04/30 16:17:23  abbott
// Summary: Commented out some useless code
// Author: JAA
//
// Revision 1.11  2011/03/11 12:36:55  bigatti
// -- changed size(DynamicBitset) into len(DynamicBitset)
//
// Revision 1.10  2010/05/28 15:50:42  bigatti
// -- cleaning
// -- commented out ctor taking a list of long
// -- moved some "facet" functions to TmpIsTree
//
// Revision 1.9  2010/04/15 16:01:39  bigatti
// -- DynamicBitset going towards final version
//
// Revision 1.8  2010/04/13 15:30:25  bigatti
// -- reorganized into (almost) final shape
//
// Revision 1.7  2010/03/31 13:53:42  bigatti
// -- naming convention: size --> mySize
//
// Revision 1.6  2010/03/30 15:36:39  bigatti
// -- using DynamicBitset.HC (with former code for facet)
//
// Revision 1.5  2007/10/30 17:14:07  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.4  2007/04/11 12:26:50  bigatti
// -- fix for Microsoft compiler
//
// Revision 1.3  2007/03/27 16:57:04  bigatti
// -- removed logging printouts
//
// Revision 1.2  2007/03/27 15:27:05  bigatti
// -- minor update for TmpIsTree + test
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.9  2007/03/08 18:22:28  cocoa
// Just whitespace cleaning.
//
// Revision 1.8  2007/03/07 14:59:02  cocoa
// -- renamed complex --> FacetComplex
//
// Revision 1.7  2007/02/10 18:44:03  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.6  2006/12/21 13:48:32  cocoa
// Made all increment/decrement calls prefix (except where the must be postfix).
//
// Revision 1.5  2006/11/24 17:12:05  cocoa
// -- reorganized includes of header files
//
// Revision 1.4  2006/11/22 14:43:32  cocoa
// -- minor cleaning (indicated by Intel compiler)
//
// Revision 1.3  2006/10/06 14:04:14  cocoa
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
