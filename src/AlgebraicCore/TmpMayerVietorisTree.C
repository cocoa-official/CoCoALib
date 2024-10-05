//   Copyright (c)  2008-2017  John Abbott and Anna Bigatti
//   Authors:  2008-2017  Anna Bigatti and Eduardo Saenz-de-Cabezon

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

#include "CoCoA/TmpMayerVietorisTree.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/SparsePolyOps-RingElem.H"
//#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/SparsePolyOps-ideal.H"
#include "CoCoA/TmpPPVector.H"
#include "CoCoA/VectorOps.H"   // for product
#include "CoCoA/ideal.H"
#include "CoCoA/ring.H"
#include "CoCoA/utils.H"
#include "CoCoA/verbose.H"

#include <iostream>
// #include <list>
using std::list;
#include <utility>
using std::pair;
// #include <vector> --- already included in TmpPPVector.H
using std::vector;

using std::endl;
/**********************************************************/
namespace CoCoA
{

  class MVTN1Node
  {
  public:
    MVTN1Node(PPMonoid PPM, DivMaskRule DMR);
    ~MVTN1Node() {};
  public:
    unsigned int pos;
    unsigned int dim;
    PPVector the_ideal;
    //varlist used_vars;
  };

  MVTN1Node::MVTN1Node(PPMonoid PPM, DivMaskRule DMR):
      the_ideal(PPM, DMR)
  { }
  

  typedef list<MVTN1Node> N1tree;



//******************************************************************//
//           TREES, NODES and BINODOS                               //
//******************************************************************//

//These types include everything needed to build a Mayer-Vietoris tree, in which each node is defined by a monomial ideal, a position and a dimension
//Binodos are the types in which the result of the computation of MAyer-Vietoris trees will be stored. A binodo has a multidegree and a list_of_dims

  class MVTNode
  {
  public:
    MVTNode(PPMonoid PPM, DivMaskRule DMR);
    ~MVTNode() {};
  public:
    position_t pos;
    unsigned int dim;
    PPVector the_ideal;
    //varlist used_vars;
  };

  MVTNode::MVTNode(PPMonoid PPM, DivMaskRule DMR):
      the_ideal(PPM, DMR)
  { }
  

  typedef list<MVTNode> tree;



  namespace // anonymous for functions local to this file/compilation unit.
  {

    void tilde(PPVector& g, PPVector& I)
    {
      PPVector m(PPM(I), DMR(I));
      PushBackPopBack(m, I);
      lcms(g,m,I);
      InterreduceSort(g);
    }
  

    bool HasHomology(const PPWithMask& pp, const PPVector& J)
    {
      PPMonoid M=owner(PP(pp));
      int n=NumIndets(M);
      PPMonoidElem ppi=PP(pp)/product(indets(M));
      if (IsDivisible(ppi, J))
        return false;
      else
        for (int i=0; i<n; ++i)
          if (!IsDivisible(indet(M,i)*ppi, J))
            return false;
      return true;
    }

// Compiler says this is never used
    // bool HasHomology_gen(const PPWithMask& pp, const PPVector& J, const PPMonoidElem& m)
    // {
    //   PPMonoid M=owner(PP(pp));
    //   int n=NumIndets(M);
    //   PPMonoidElem ppi=(PP(pp)*m)/product(indets(M));
    //   if (IsDivisible(ppi, J))
    //     return false;
    //   else
    //     for (int i=0; i<n; ++i)
    //       if(exponent(m,i)==0)
    //       {
    //         if (!IsDivisible(indet(M,i)*ppi, J))
    //         {return false;}
    //       }
    //   return true;
    // }
    /////////////// new code 4-july-2008
//This file contains the basic structures that will be used in the computation of Mayer-Vietoris trees and the basic functions involving them.


//*******************************************************//
//*******************************************************//
//          F U N C T I O N S                            //
//*******************************************************//
//*******************************************************//
//*******************************************************//





// FUNCTIONS FOR THE CONTAINER OF MULTIDEGREES, DIMENSIONS AND POSITIONS:   //

    MultidegreeMap::iterator special_insert(const PPMonoidElem& , const int& , const position_t& , MultidegreeMap& /* ,MultidegreeMap::iterator*/);

// void put(MultidegreeMap&);

// void put(ListOfDims&);

// AUXILIARY FUNCTIONS FOR INTEGERS AND LISTS OF INTEGERS //

    int number_zeros(list<int> );

    void binary_expression(position_t , list<int>& );


//const int N=100;//NUMBER OF VARIABLES





// FUNCTIONS FOR THE CONTAINER OF MULTIDEGREES, DIMENSIONS AND POSITIONS:   //



    MultidegreeMap::iterator special_insert(const PPMonoidElem& multidegree, const int& dimension, const position_t& position, MultidegreeMap& the_list /*, MultidegreeMap::iterator insertion_hint*/)
    {

      list<position_t> thispos(1,position);

      ListOfDims to_insert;
      to_insert.insert(std::make_pair(dimension,thispos));

      std::pair<MultidegreeMap::iterator, bool> insertion=the_list.insert(/*insertion_hint,*/ std::make_pair(multidegree,to_insert));//WHAT ABOUT THE INSERTION HINT???

      MultidegreeMap::iterator insert_point=insertion.first;

      if (!insertion.second)
      {
	pair<ListOfDims::iterator, bool> inner_insertion=(insert_point->second).insert(std::make_pair(dimension,thispos)); //repeated, use a variable!
	if (!inner_insertion.second)
        {
          (inner_insertion.first->second).push_front(position);
        }	
		
      }
      return insert_point;
    }



    void output(std::ostream& out, const ListOfDims& the_list)
    {
      if (!out) return;  // short-cut for bad ostreams

      ListOfDims::const_iterator pos;
      for (pos= the_list.begin(); pos!=the_list.end();++pos)
      {
	out<<pos->first<<":"<<pos->second<<endl;
      }
    }


    void output(std::ostream& out, const MultidegreeMap& the_list)
    {
      if (!out) return;  // short-cut for bad ostreams

      MultidegreeMap::const_iterator pos;
      for (pos= the_list.begin(); pos!=the_list.end();++pos)
      {
	out<<pos->first<<endl;
	output(out, pos->second);
	out<<"size of the list:"<<len(pos->second)<<endl;
      }
    }


// AUXILIARY FUNCTIONS FOR INTEGERS AND LISTS OF INTEGERS //


    int number_zeros(list<int> L)
    {
      int n=0;
      list<int>::iterator p = L.begin();
      while (p!=L.end())
      {
	if(*p==0)
        {
          ++n;
        }
	++p;
      }
      return n;
    }

    void binary_expression(position_t n, list<int>& L)
    {
      if (n<2)
      { if (n == 0) L.push_front(0); else L.push_front(1);}
      else
      {
        if (n%2 == 0) L.push_front(0); else L.push_front(1);
        binary_expression(n/2,L);
      }
    }


//----------------------------------------------------------------------

    bool test_pos(position_t p1, position_t p2)
    {
      list<int> b1,b2;

      binary_expression(p1,b1);
      binary_expression(p2,b2);

      while (b1.front()==b2.front())
      {
	b1.pop_front();
	b2.pop_front();
      }
      if(b1.front()==0)
      {
	return(number_zeros(b1)-1==number_zeros(b2));
      }
      else
      {
        return(number_zeros(b2)-1==number_zeros(b1));
      }
    }


// 20140809 JAA commented out because never used.
// void remove_pairs(list<position_t>& l1, list<position_t>& l2)
// {
//  std::list<position_t>::iterator p=l1.begin();
//  while (p!=l1.end())
//  	{
//  	std::list<position_t>::iterator q=l2.begin();
//  	bool found= false;
// 	while (!found && q!=l2.end())
//  		{
// 		if(test_pos(*p,*q))
// 			{
// 			found=true;
// 			p=l1.erase(p);
// 			q=l2.erase(q);
// 			}
// 			else
// 				{
// 				++q;
// 				}
// 		}
// 	if(!found && p!=l1.end())
// 		{
// 		++p;
// 		}
// 	}

// }

    void reduce_list (ListOfDims& the_list, ListOfDims& undecided_list) //reduces the list according to the position criterion 
    //VERY inefficient version!!
    {
      ListOfDims decided;
//list<position_t> this_undecided, this_decided;
      long k=len(the_list);
      int l=0;
      int m=0;
//FIRST TEST 
//If the multidegree only appears twice in the tree, then it is a reduction pair iff the positions are couple
// if (k==2)
// 	{
// 	l=len((the_list.begin())->second);//cout<<"l="<<l;
// 	m=len((++the_list.begin())->second);//cout<<"m="<<m<<endl;
// 
// 	if (l==1 and m==1)
// 		{
// 		if (test_pos((the_list.begin())->second.front(),(++the_list.begin())->second.front()))
// 			{
// 			the_list.clear();
// 			}
// 		}
// 	}
//SECOND TEST
//Every multidegree with no couples cannot be part of a reduction pair
      if (len(the_list)>1)
      {
	ListOfDims::iterator p=the_list.begin();
	
	while (p!=the_list.end())
        {
          int i=p->first;
		
          std::list<position_t>::iterator r=p->second.begin();
          while (r!=p->second.end())
          {
            bool found=false;
            ListOfDims::iterator q=the_list.begin();		
            while (q!=the_list.end())
            {
              if(q->first==i+1 || q->first==i-1)
              {
                std::list<position_t>::iterator s=q->second.begin();
                while(s!=q->second.end() && !found)
                {
                  if (test_pos(*r,*s))
                  {
                    list<int> b1,b2;
                    binary_expression(*r,b1);
                    binary_expression(*s,b2);
                    //cout<<"test da 1"<<"bin"<<*r<<"="<<b1<<"::"<<*s<<"="<<b2<<endl;
                    found=true;
                    undecided_list[p->first].push_front(*r);
                  }
                  else
                  {
                    ++s;
                  }
                }
                ++q;
              }
              else
              {
                ++q;
              }
            }
            if(!found)
            {
              decided[p->first].push_front(*r);
            }
            ++r;
          }
          //swap(undecided_list[p->first],p->second);
// 		undecided_list.insert(std::make_pair(p->first,this_undecided));
          //decided.insert(std::make_pair(p->first,this_decided));
		
          ++p;
        }
	swap(the_list,decided);
      }

//THIRD TEST
//If the list has two dimension and one has more elements than the other, some of these will have no pair
//Warning! no position-based tests allowed after this
      k=len(undecided_list);
      if (k==2)
      {
        l=len((undecided_list.begin())->second);//cout<<"l und="<<l;
        m=len((++undecided_list.begin())->second);//cout<<"m und="<<m<<endl;
	if (l>m)
        {
          for (int i=0; i < l-m; ++i)
          {
            the_list[(undecided_list.begin())->first].push_front((undecided_list.begin())->second.front());
            (undecided_list.begin())->second.pop_front();
          }
        }
	if (m>l)
        {
          for (int i=0; i < m-l; ++i)
          {
            the_list[(++undecided_list.begin())->first].push_front((++undecided_list.begin())->second.front());
            (++undecided_list.begin())->second.pop_front();
          }
        }
      }
    }


//////////////////



  }  // end of anonymous namespace

  std::ostream& operator<<(std::ostream& out, const MultidegreeMap& L)
  {
    if (!out) return out;  // short-cut for bad ostreams
    output(out, L);
    return out;
  }



  void ReduceMVTree(MultidegreeMap& the_map, MultidegreeMap& undecided_map)
  {
    MultidegreeMap::iterator pos;
    MultidegreeMap::iterator undec_pos=undecided_map.begin();
    for (pos=the_map.begin(); pos!=the_map.end();)
    {
      ListOfDims undecided_list;
      //cout<<"lista inicial";put(pos->second);cout<<"::";
      reduce_list(pos->second, undecided_list);
	
      if (!undecided_list.empty())
      {
        //cout<<pos->first<<endl;
        //cout<<"reducida ";put(pos->second);cout<<":: undecided ";put(undecided_list);cout<<endl;
        undecided_map.insert(undec_pos,std::make_pair(pos->first,undecided_list));
      }
      if (pos->second.empty())
      {
        the_map.erase(pos++);
      }
      else{++pos;}
    }

  }
  

  void MayerVietorisTreeN1(PPVector& betti, const PPVector& I)
  {
    PPVector N1(PPM(I), DMR(I));
    N1tree T;
    unsigned int N=NumIndets(PPM(betti));
    PPVector J=I;
    InterreduceSort(J);
    MVTN1Node nodo(PPM(J), DMR(J));
    swap(nodo.the_ideal, J);
    nodo.pos=1; nodo.dim=0;
    T.push_back(nodo);


    if (len(J)==1) return;

    N1tree undone=T;
    MVTN1Node new_node(PPM(J), DMR(J));
    int mipos, midim;
    PPVector mi_ideal(PPM(J), DMR(J));
    
    while (!undone.empty())
    {
      swap(mi_ideal, undone.front().the_ideal);
      mipos=undone.front().pos;
      midim=undone.front().dim;
      undone.pop_front();
      
      tilde(new_node.the_ideal, mi_ideal);
      
      //left child
      new_node.pos=mipos*2;
      new_node.dim=midim+1;
      
      if (new_node.dim==N-1)
      {for (long i=0; i<len(new_node.the_ideal); ++i)
          if (HasHomology(new_node.the_ideal[i], I))
            N1.myPushBack(new_node.the_ideal[i]);
      }
      else
      {
        if (len(new_node.the_ideal)>1 && len(new_node.the_ideal)+new_node.dim-1>=N-1 /* && HAS FULL SUPPORT*/)
          undone.push_back(new_node);
      }
      //right child
      new_node.pos = 1+mipos*2;
      new_node.dim = midim;
      swap(new_node.the_ideal, mi_ideal);
      if (len(new_node.the_ideal)>1 && len(new_node.the_ideal)+new_node.dim-1>=N-1)
        undone.push_back(new_node);
    }
    swap(betti, N1);
  }


//   void MayerVietorisTreeN1_gen(PPVector& betti, const PPVector& I, const PPMonoidElem& m, int cont)
//   {
// 
//     PPVector N1(PPM(I), DMR(I));
//     N1tree T;
//     unsigned int N=NumIndets(PPM(betti));
//     PPVector J=I;
//     InterreduceSort(J);
//     MVTN1Node nodo(PPM(J), DMR(J));
//     swap(nodo.the_ideal, J);
//     nodo.pos=1; nodo.dim=0;
//     T.push_back(nodo);
// 
// 
//     if (J.mySize()==1) 
// 	{
// 	if (N-cont==1)
// 		{betti.myPushBack(J[0]);}
// 		else{return;}
// 	}
// 
// 
//     N1tree undone=T;
//     MVTN1Node new_node(PPM(J), DMR(J));
//     int mipos, midim;
//     PPVector mi_ideal(PPM(J), DMR(J));
//     
//     while (!undone.empty())
//     {
//       swap(mi_ideal, undone.front().the_ideal);
//       mipos=undone.front().pos;
//       midim=undone.front().dim;
//       undone.pop_front();
//       
//       tilde(new_node.the_ideal, mi_ideal);
//       
//       //left child
//       new_node.pos=mipos*2;
//       new_node.dim=midim+1;
//       
//       if (new_node.dim==N-cont-1)
//         {for (unsigned int i=0; i<new_node.the_ideal.mySize(); ++i)
//           if (HasHomology_gen(new_node.the_ideal[i], I,m))
//             N1.myPushBack(new_node.the_ideal[i]);
// 	}
// 	else
// 	{
//    	   if (new_node.the_ideal.mySize()>1 && new_node.the_ideal.mySize()+new_node.dim-1>=N-cont-1)
//      	   undone.push_back(new_node);
// 	}
// 
//       //right child
//       new_node.pos = 1+mipos*2;
//       new_node.dim = midim;
//       swap(new_node.the_ideal, mi_ideal);
//       if (new_node.the_ideal.mySize()>1 && new_node.the_ideal.mySize()+new_node.dim-1>=N-cont-1)
//         undone.push_back(new_node);
//     }
//     swap(betti, N1);
//   }



  void SimplicialMayerVietorisTree(vector<int>& ranks, const PPVector& I)
  {
    N1tree T;
    //    unsigned int N=NumIndets(PPM(I));
    PPVector J=I;
    InterreduceSort(J);
    MVTN1Node nodo(PPM(J), DMR(J));
    swap(nodo.the_ideal, J);
    nodo.pos=1; nodo.dim=0;
    T.push_back(nodo);

    if (len(J)==1) return;

    N1tree undone=T;
    MVTN1Node new_node(PPM(J), DMR(J));
    int mipos, midim;
    PPVector mi_ideal(PPM(J), DMR(J));
    PPWithMask prod_vars(product(indets(PPM(J))),DMR(J));

    if (len(J)==1) return;
    
    while (!undone.empty())
    {
      swap(mi_ideal, undone.front().the_ideal);
      mipos=undone.front().pos;
      midim=undone.front().dim;
      undone.pop_front();
      
      tilde(new_node.the_ideal, mi_ideal);
      
      //left child
      new_node.pos=mipos*2;
      new_node.dim=midim+1;

//TO MAKE THE PRODUCT AND TEST IF THERE IS FULL SUPPORT
	
      if (len(new_node.the_ideal)==1)
      {
        ranks[new_node.dim]++;
      }
      else	{
        undone.push_back(new_node);
      }
	
        
      //right child
      swap(new_node.the_ideal, mi_ideal);
      PPMonoidElem my_producto=PP(new_node.the_ideal[0]);
      for (long k=1; k<len(new_node.the_ideal); ++k)
        my_producto=my_producto*PP(new_node.the_ideal[k]);
// 	std::cout<<new_node.the_ideal<<endl;      
// 	std::cout<<my_producto<<endl;
// 	std::cout<<prod_vars<<endl;
// 	std::cout<<"::"<<IsDivisibleFast(PPWithMask(my_producto,DMR(J)), prod_vars)<<endl;
      new_node.pos = 1+mipos*2;
      new_node.dim = midim;

      
      if (len(new_node.the_ideal)>1 && IsDivisibleFast(PPWithMask(my_producto,DMR(J)), prod_vars))
        undone.push_back(new_node);

    }

  }
  void MayerVietorisTree(MultidegreeMap& output_list, const PPVector& I)
  {
    VerboseLog VERBOSE("MayerVietorisTree");
    VERBOSE(50) << " starting " << std::endl;
    tree undone;
    MVTNode nodo(PPM(I), DMR(I));
    MultidegreeMap::iterator insert_pos;
    PPVector J=I;
    InterreduceSort(J);
    VERBOSE(50) << " after InterreduceSort " << std::endl;
    swap(nodo.the_ideal, J);
    nodo.pos=1; nodo.dim=0;
    undone.push_back(nodo);
    for (long i=0; i < len(nodo.the_ideal); ++i)    //inserting all the relevant monomials in a binary tree
    {
      //    insert_pos=
      special_insert(PP(nodo.the_ideal[i]),nodo.dim,nodo.pos, output_list /* ,insert_pos*/);
    }
    if (len(nodo.the_ideal)>1) 
    {

      MVTNode minodo(PPM(J), DMR(J));
      MVTNode new_node(PPM(J), DMR(J));
      position_t mipos;
      int midim;
      PPVector mi_ideal(PPM(J), DMR(J));
      PPVector mi_tilde(PPM(J), DMR(J));
      PPVector mi_prime(PPM(J), DMR(J));

      while (!undone.empty())
      {
        swap(mi_ideal,undone.front().the_ideal);
        mipos=undone.front().pos; midim=undone.front().dim;
        undone.pop_front();


        tilde(new_node.the_ideal,mi_ideal);

        //left child
        new_node.pos=mipos*2;
        new_node.dim=midim+1;
        //first insert the first monomial of the node, and obtain the position in which it was inserted
        //      insert_pos=
        special_insert(PP(new_node.the_ideal[0]),new_node.dim,new_node.pos, output_list /*, output_list.begin()*/);
        // 		
        for (long i=1; i < len(new_node.the_ideal); ++i)    //inserting all the relevant monomials in a binary tree
        {	
          //        insert_pos=
          special_insert(PP(new_node.the_ideal[i]),new_node.dim,new_node.pos, output_list /* ,insert_pos*/);

        }; 
		
        if (len(new_node.the_ideal)>2) {undone.push_back(new_node);}
        else
        {
          if (len(new_node.the_ideal)==2)
          {
            //          insert_pos=
            special_insert(lcm(PP(new_node.the_ideal[0]),PP(new_node.the_ideal[1])),1+new_node.dim,2*new_node.pos, output_list /* ,insert_pos*/);
          }
        }


        //right child
        new_node.pos=1+mipos*2;
        new_node.dim=midim;
        swap(new_node.the_ideal,mi_ideal);
		
        if (len(new_node.the_ideal)>2) {undone.push_back(new_node);}
        else
        {
          if (len(new_node.the_ideal)==2)
          {
            //          insert_pos=
            special_insert(lcm(PP(new_node.the_ideal[0]),PP(new_node.the_ideal[1])),1+new_node.dim,2*new_node.pos, output_list /* ,insert_pos*/);
          }
        }
      
      }


    }
    VERBOSE(50) << " end " << std::endl;
  }




  void MayerVietorisTreeExtremal(MultidegreeMap& output_list, list<PPMonoidElem>& /*extremal*/, const PPVector& I) //for the CoCoA school 2009
  {

    tree undone;
    MVTNode nodo(PPM(I), DMR(I));

    PPVector J=I;
    InterreduceSort(J);
    swap(nodo.the_ideal, J);
    nodo.pos=1; nodo.dim=0;
    undone.push_back(nodo);

    if (len(nodo.the_ideal)>1)
    {

      MVTNode minodo(PPM(J), DMR(J));
      MVTNode new_node(PPM(J), DMR(J));
      position_t mipos;
      int midim;
      PPVector mi_ideal(PPM(J), DMR(J));
      PPVector mi_tilde(PPM(J), DMR(J));
      PPVector mi_prime(PPM(J), DMR(J));

      while (!undone.empty())
      {
        swap(mi_ideal,undone.front().the_ideal);
        mipos=undone.front().pos; midim=undone.front().dim;
        undone.pop_front();


        tilde(new_node.the_ideal,mi_ideal);

        //left child
        new_node.pos=mipos*2;
        new_node.dim=midim+1;
        //first insert the first monomial of the node, and obtain the position in which it was inserted
        MultidegreeMap::iterator insert_pos=special_insert(PP(new_node.the_ideal[0]),new_node.dim,new_node.pos, output_list /*, output_list.begin()*/);
        // 		
        for (long i=1; i < len(new_node.the_ideal); ++i)    //inserting all the relevant monomials in a binary tree
        {	
          insert_pos=special_insert(PP(new_node.the_ideal[i]),new_node.dim,new_node.pos, output_list /* ,insert_pos*/);

        }; 
		
        if (len(new_node.the_ideal)>1) {undone.push_back(new_node);}


        //right child
        new_node.pos=1+mipos*2;
        new_node.dim=midim;
        swap(new_node.the_ideal,mi_ideal);
		
        if (len(new_node.the_ideal)>1) {undone.push_back(new_node);}
		
      }


    }
  }







  void ReducedMayerVietorisTree(MultidegreeMap& output_list, const PPVector& I)
  {
    MayerVietorisTree(output_list, I);
    MultidegreeMap  undecided;
    ReduceMVTree(output_list, undecided);
  }


  void Bettis(vector<int>& bettis, const MultidegreeMap& the_list)
  {
    VerboseLog VERBOSE("Bettis");
    MultidegreeMap::const_iterator pos;
    for (pos= the_list.begin(); pos!=the_list.end();++pos)
    {
      ListOfDims::const_iterator inner_pos;
      for (inner_pos= (pos->second).begin(); inner_pos!=(pos->second).end();++inner_pos)
      {
        bettis[inner_pos->first] += len(inner_pos->second);
      }
    }
  }


  int ResSize(const vector<int>& R)
  {
    vector<int>::const_iterator pos;
    int my_size=0;
    for(pos=R.begin(); pos!=R.end();++pos)
    {
      my_size+=*pos;
    }
    return my_size;
  }




  void GradedBettis(const MultidegreeMap& decided,BettiPseudoDiagram& graded_bettis, int& reg, int& pd)
  {
    reg=-1;pd=-1;
    if (!decided.empty())
    {
      MultidegreeMap::const_iterator pos;
      pd=0;
      for (pos= decided.begin(); pos!=decided.end();++pos)
      {
	int deg=StdDeg(pos->first);
	
	ListOfDims::const_iterator inner_pos;
	for (inner_pos= (pos->second).begin(); inner_pos!=(pos->second).end();++inner_pos)
        {
          if (graded_bettis[deg-(inner_pos->first)].empty())
          {
            int N=NumIndets(owner(pos->first));
            graded_bettis[deg-(inner_pos->first)].resize(N);
          }
		
          graded_bettis[deg-(inner_pos->first)][inner_pos->first]+=len(inner_pos->second);
          if(inner_pos->first > pd)
          {
            pd=inner_pos->first;
          }
        }

      }
      reg=(--graded_bettis.end())->first;
    }
  }

  int regularity(const MultidegreeMap& decided)
  {
    int reg=0;
    if (!decided.empty())
    {
      MultidegreeMap::const_iterator pos;
      for (pos= decided.begin(); pos!=decided.end();++pos)
      {
	int deg=StdDeg(pos->first);
	ListOfDims::const_iterator inner_pos;
	for (inner_pos= (pos->second).begin(); inner_pos!=(pos->second).end();++inner_pos)
        {
          if(deg-(inner_pos->first)>reg)
            reg=deg-inner_pos->first;
        }

      }
    }
    return reg;
  }



// int RegularityFromGradedBetti( const BettiPseudoDiagram& diagram)
// {
// return (--diagram.end())->first;
// }

  void PrintBettiDiagram(const BettiPseudoDiagram& graded_bettis)
  {
    BettiPseudoDiagram::const_iterator pos;
    for (pos= graded_bettis.begin(); pos!=graded_bettis.end();++pos)
    {
      std::cout<<pos->first<<" : "<<pos->second<<endl;
    }
  }

  void PseudoBettiDiagram(const BettiPseudoDiagram& decided, const BettiPseudoDiagram& undecided)
  {
    BettiPseudoDiagram::const_iterator pos, posdec;
    posdec=decided.begin();
    pos= undecided.begin();
    while (posdec!=decided.end() && pos!=undecided.end())
    {
      if(posdec->first<pos->first)
      {
        std::cout<<posdec->first<<" : ";
        for (long i=0; i!=len(posdec->second); ++i)
          std::cout<<posdec->second[i]<<" , ";
        std::cout<<endl;
        ++posdec;
      }
      else 
      {
        if(pos->first<posdec->first)
        {
          std::cout<<pos->first<<" : ";
          for (long i=0; i!=len(pos->second); ++i)
          {
            std::cout<<"0-"<<pos->second[i]<<" , ";	
          }
          std::cout<<endl;
          ++pos;
        }
      }
      while (pos->first==posdec->first)
      {
        std::cout<<pos->first<<" : ";
        for (long i=0; i!=len(pos->second); ++i)
        {
          if(pos->second[i]>0)
          {
            std::cout<<posdec->second[i]<<"-"<<posdec->second[i]+pos->second[i]<<" , ";
          }
          else
          {
            std::cout<<posdec->second[i]<<" , ";
          }
        }
	std::cout<<endl;
	++posdec;
        ++pos;
      }
// std::cout<<(posdec==decided.end())<<endl;
// std::cout<<"pos: "<<pos->first<<"posdec: "<<posdec->first<<endl;
      //Si hemos llegado al final de decided, imprimimos undecided
      if(posdec==decided.end())
      {
        for (/*pos*/; pos!=undecided.end(); ++pos)
	{
          std::cout<<pos->first<<" : ";
          for (long i=0; i!=len(pos->second); ++i)
          {
            if(pos->second[i]>0)
            {
              std::cout<<"0-"<<pos->second[i]<<" , ";
            }
            else
            {
              std::cout<<"0 , ";
            }
          }
          std::cout<<endl;
	}
      }
      else //si hemos llegado al final de undecided, imprimimos decided
      {if(pos==undecided.end())
        {
          for (/*posdec*/; posdec!=decided.end(); ++posdec)
          {
            std::cout<<posdec->first<<" : ";
            const long n = len(posdec->second);
            for (long i=0; i!=n; ++i)
              std::cout<<posdec->second[i]<<" , ";	
            std::cout<<endl;
          }
        }
      }
    }

	
  }


  void PrintMultidegreeMap(const MultidegreeMap& myMap)
  {
    MultidegreeMap::const_iterator pos;
    for (pos= myMap.begin(); pos!=myMap.end();++pos)
    {
      std::cout<<"Multidegree: "<<pos->first<<" Degree: "<<StdDeg(pos->first); //<<endl;
      ListOfDims::const_iterator inner_pos;
      for (inner_pos= (pos->second).begin(); inner_pos!=(pos->second).end();++inner_pos)
      {
        std::cout<<"  dim="<<inner_pos->first<<", pos="<<inner_pos->second<<endl;
      }

    }
  }

  void PrintMultidegreeMapDim(const MultidegreeMap& myMap)
  {
    MultidegreeMap::const_iterator pos;
    vector<list<position_t> > thispositions (NumIndets(owner(myMap.begin()->first)));
    vector<list<PPMonoidElem> > thismultis(NumIndets(owner(myMap.begin()->first)));
    for (pos= myMap.begin(); pos!=myMap.end();++pos)
    {
      PPMonoidElem thismulti=pos->first;
      ListOfDims::const_iterator inner_pos;
      for (inner_pos= (pos->second).begin(); inner_pos!=(pos->second).end();++inner_pos)
      {
        list<position_t>::const_iterator dims_pos;
        for (dims_pos= (inner_pos->second).begin(); dims_pos!=(inner_pos->second).end();++dims_pos)
        {
          thispositions[inner_pos->first].push_front(*dims_pos);
          thismultis[inner_pos->first].push_front(thismulti);
        }
// 		std::cout<<"Dimension: "<<inner_pos->first<<", Positions: "<<inner_pos->second<<endl;
      }

    }
    int cont=0;
    while (!thispositions[cont].empty())
    {
      std::cout<<"Dimension: "<<cont<<" ( "<<len(thispositions[cont])<< " )"<<endl;
      list<position_t>::const_iterator pos;
      list<PPMonoidElem>::const_iterator pos_mult=thismultis[cont].begin();
      for (pos=thispositions[cont].begin(); pos!=thispositions[cont].end();++pos)
      {
        std::cout<<(*pos)<<"--"<<(*pos_mult)<<endl;
        ++pos_mult;
      }
      ++cont;
    }
  }

  void irreducible(vector<PPVector>& f, const PPVector& components, const PPMonoidElem& lambda)
  {
    const PPMonoid P=PPM(components);

    const int N=NumIndets(P);

    for (long i=0; i<len(components); ++i)
    {
      PPVector g(PPM(components), DMR(components));
      for (int j=0; j<N; ++j)
      {
        if(exponent(PP(components[i]),j)<=exponent(lambda,j))
        {
          g.myPushBack(power(indet(P,j),exponent(PP(components[i]),j)));
        }
      }
      f.push_back(g);
    }
  }


//INTERFACE FUNCTIONS FOR THE MVT PACKAGE

//1- OUTPUT OF MAYER-VIETORIS TREES
//1a- THE OUTPUT IS BY MULTIDEGREE
  void MVTPrint(const PPVector& I)
  {
    MultidegreeMap output_list;
    MayerVietorisTree(output_list, I);
    PrintMultidegreeMap(output_list);
  }


//1b- THE OUTPUT IS BY DIMENSIONS

  void MVTPrintDim(const PPVector& I)
  {
    MultidegreeMap output_list;
    MayerVietorisTree(output_list, I);
    PrintMultidegreeMapDim(output_list);
  }


//2- OUTPUT PSEUDO BETTI DIAGRAM
  void MVTPseudoBettiDiagram(const PPVector& I)
  {
    MultidegreeMap decided, undecided;
    BettiPseudoDiagram dec,undec;
    int reg1,pd1,reg2,pd2;
    MayerVietorisTree(decided, I);
    ReduceMVTree(decided,undecided);
    GradedBettis(decided,dec,reg1,pd1);
    GradedBettis(undecided,undec,reg2,pd2);
    PseudoBettiDiagram(dec, undec);
  }

//3- REGULARITY AND BOUNDS
  int MVTRegularityLowerBound(const PPVector& I)
  {
    MultidegreeMap decided, undecided;
    BettiPseudoDiagram dec;
    int reg1,pd1;
    MayerVietorisTree(decided, I);
    ReduceMVTree(decided,undecided);
    GradedBettis(decided,dec,reg1,pd1);
    return reg1;
  }
  int MVTRegularityUpperBound(const PPVector& I)
  {
    MultidegreeMap decided, undecided;
    BettiPseudoDiagram dec,undec;
    int reg1,pd1,reg2,pd2;
    MayerVietorisTree(decided, I);
    ReduceMVTree(decided,undecided);
    GradedBettis(decided,dec,reg1,pd1);
    GradedBettis(undecided,undec,reg2,pd2);
    if (reg2>reg1) return reg2;
    else return reg1;
  }

  int MVTRegularity(const PPVector& I)
  {
    MultidegreeMap decided, undecided;
    BettiPseudoDiagram dec,undec;
    int reg1,pd1,reg2,pd2;
    MayerVietorisTree(decided, I);
    ReduceMVTree(decided,undecided);
    GradedBettis(decided,dec,reg1,pd1);
    GradedBettis(undecided,undec,reg2,pd2);
    if (reg2>reg1) {return -1;}
    else return reg1;
  }


// int MVTRegularity(const PPVector& I) //Innefficient, computes twice!
// {
// int r1,r2;
// r1=MVTRegularityUpperBound(I);
// r2=MVTRegularityLowerBound(I);
// if (r1>r2) return -1;
//   else return r1;
// }



//4- PROJECTIVE DIMENSION AND BOUNDS
  int MVTProjDimLowerBound(const PPVector& I)
  {
    MultidegreeMap decided, undecided;
    BettiPseudoDiagram dec;
    int reg1,pd1;
    MayerVietorisTree(decided, I);
    ReduceMVTree(decided,undecided);
    GradedBettis(decided,dec,reg1,pd1);
    return pd1;
  }
  int MVTProjDimUpperBound(const PPVector& I)
  {
    MultidegreeMap decided, undecided;
    BettiPseudoDiagram dec,undec;
    int reg1,pd1,reg2,pd2;
    MayerVietorisTree(decided, I);
    ReduceMVTree(decided,undecided);
    GradedBettis(decided,dec,reg1,pd1);
    GradedBettis(undecided,dec,reg2,pd2);
    if (pd2>pd1) return pd2;
    else return pd1;
  }

  int MVTProjDim(const PPVector& I)
  {
    MultidegreeMap decided, undecided;
    BettiPseudoDiagram dec,undec;
    int reg1,pd1,reg2,pd2;
    MayerVietorisTree(decided, I);
    ReduceMVTree(decided,undecided);
    GradedBettis(decided,dec,reg1,pd1);
    GradedBettis(undecided,dec,reg2,pd2);
    if (pd2>pd1) return -1;
    else return pd1;
  }

//---- handy interface -----------------------------------------------
  std::vector<RingElem> MayerVietorisTreeN1(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))  CoCoA_THROW_ERROR1(ERR::ReqSparsePolyRing);
    if (!AreGensMonomial(I))  CoCoA_THROW_ERROR1(ERR::ReqMonomialGens);
    ///  if (!gens(I).empty())
    {
      // convert input into PPVector, operate, convert back
      const SparsePolyRing P = RingOf(I);
      PPVector gensI(PPM(P), NewDivMaskEvenPowers());
      const std::vector<RingElem>& l1 = gens(I);
      for (vector<RingElem>::const_iterator it=l1.begin(); it!=l1.end() ; ++it)
        PushBack(gensI, LPP(*it));
      PPVector betti(PPM(P), DMR(gensI));
      MayerVietorisTreeN1(betti, gensI);
      std::vector<RingElem> res;
      convert(res, P, betti);
      return res;
    }
  }


/*
  -- by Eduardo Saenz de Cabezon -------------------------------------
  Define MVT5(X)
  Computation := "Mayer_Vietoris_Tree";
  Arguments := [Tagged(Gens(X),"PP_LIST")];
  Return $cocoa5.OperationCommunication(Computation, Arguments);
  EndDefine; -- MVT5

  Define MayerVietorisTreeN1(X)
  Computation := "MVT_N_minus_one";
  Arguments := [Tagged(Gens(X),"PP_LIST")];
  Return $cocoa5.OperationCommunication(Computation, Arguments);
  EndDefine; -- MayerVietorisTreeN1

  Define MVTReg5(X)
  Computation := "MVT_Regularity";
  Arguments := [Tagged(Gens(X),"PP_LIST")];
  Return $cocoa5.OperationCommunication(Computation, Arguments);
  EndDefine; -- MVTReg5

  Define MVTRegUpperBound5(X)
  Computation := "MVT_Regularity_Upper_Bound";
  Arguments := [Tagged(Gens(X),"PP_LIST")];
  Return $cocoa5.OperationCommunication(Computation, Arguments);
  EndDefine; -- MVTRegUpperBound5

  Define MVTRegLowerBound5(X)
  Computation := "MVT_Regularity_Lower_Bound";
  Arguments := [Tagged(Gens(X),"PP_LIST")];
  Return $cocoa5.OperationCommunication(Computation, Arguments);
  EndDefine; -- MVTRegLowerBound5

  Define MVTProjDim5(X)
  Computation := "MVT_ProjDim";
  Arguments := [Tagged(Gens(X),"PP_LIST")];
  Return $cocoa5.OperationCommunication(Computation, Arguments);
  EndDefine; -- MVTProjDim5

  Define MVTProjDimUpperBound5(X)
  Computation := "MVT_ProjDim_Upper_Bound";
  Arguments := [Tagged(Gens(X),"PP_LIST")];
  Return $cocoa5.OperationCommunication(Computation, Arguments);
  EndDefine; -- MVTProjDimUpperBound5

  Define MVTProjDimLowerBound5(X)
  Computation := "MVT_ProjDim_Lower_Bound";
  Arguments := [Tagged(Gens(X),"PP_LIST")];
  Return $cocoa5.OperationCommunication(Computation, Arguments);
  EndDefine; -- MVTProjDimLowerBound5

  Define IsMonomialIrreducible5(X)
  Computation := "Monomial_Is_Irreducible";
  Arguments := [Tagged(Gens(X),"PP_LIST")];
  Return $cocoa5.OperationCommunication(Computation, Arguments);
  EndDefine; -- MonomialIsIrreducible5

  Define IsMonomialPrime5(X)
  Computation := "Monomial_Is_Prime";
  Arguments := [Tagged(Gens(X),"PP_LIST")];
  Return $cocoa5.OperationCommunication(Computation, Arguments);
  EndDefine; -- MonomialIsPrime5

  Define IsMonomialPrimary5(X)
  Computation := "Monomial_Is_Primary";
  Arguments := [Tagged(Gens(X),"PP_LIST")];
  Return $cocoa5.OperationCommunication(Computation, Arguments);
  EndDefine; -- MonomialIsPrimary5


  Define MappingConeRes5(X)
  Computation := "Effective_Mapping_Cone";
  Arguments := [Tagged(Gens(X),"PP_LIST")];
  Return $cocoa5.OperationCommunication(Computation, Arguments);
  EndDefine; -- MappingConeRes5
*/

}  // end of namespace CoCoA

// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/TmpMayerVietorisTree.C,v 1.35 2024/08/01 16:55:51 bigatti Exp $
// $Log: TmpMayerVietorisTree.C,v $
// Revision 1.35  2024/08/01 16:55:51  bigatti
// Summary: error codes: Req(Elem)SparsePolyRing;
// now using CoCoA_THROW_ERROR1/2
//
// Revision 1.34  2024/07/23 16:43:34  bigatti
// Summary: error codes: ERR::ReqMonomialGens
//
// Revision 1.33  2022/02/18 14:12:01  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.32  2021/07/19 11:25:37  abbott
// Summary: Layout
//
// Revision 1.31  2020/06/17 15:49:29  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.30  2020/02/28 17:33:54  bigatti
// -- added some verbosity
//
// Revision 1.29  2020/02/11 16:56:42  abbott
// Summary: Corrected last update (see redmine 969)
//
// Revision 1.28  2020/02/11 16:12:20  abbott
// Summary: Added some checks for bad ostream (even to mem fns myOutput); see redmine 969
//
// Revision 1.27  2018/05/18 16:38:51  bigatti
// -- added include SparsePolyOps-RingElem.H
//
// Revision 1.26  2018/05/18 12:24:47  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.25  2018/05/17 15:55:26  bigatti
// -- renamed VectorOperations --> VectorOps
//
// Revision 1.24  2018/04/18 14:29:42  abbott
// Summary: Commented out unused fns; changed indentation
//
// Revision 1.23  2018/03/15 14:18:06  bigatti
// -- added files SparsePolyOps-ideal.H and SparsePolyOps-involutive.H
//
// Revision 1.22  2017/03/15 16:42:27  bigatti
// -- changed copyright
//
// Revision 1.21  2016/06/06 11:46:11  bigatti
// -- minor change in printout
//
// Revision 1.20  2016/05/24 07:59:09  bigatti
// -- minor changes in special_insert
//
// Revision 1.19  2016/05/23 12:47:55  bigatti
// -- resurrected ReducedMayerVietorisTree
//
// Revision 1.18  2014/09/04 12:17:05  bigatti
// -- handy interface with ideal (should it be here or in MonomialIdeal?)
//
// Revision 1.17  2014/08/16 13:24:01  abbott
// Summary: Commented out an apparently unused fn
// Author: JAA
//
// Revision 1.16  2014/07/31 14:45:19  abbott
// Summary: Merged io.H and UtilsTemplate.H into new header VectorOperations.H
// Author: JAA
//
// Revision 1.15  2014/07/14 15:08:31  abbott
// Summary: Changed include of tmp.H into UtilsTemplate.H
// Author: JAA
//
// Revision 1.14  2014/04/30 16:24:15  abbott
// Summary: Replaced X.size() by len(X); Replaced size_t by long
// Author: JAA
//
// Revision 1.13  2012/05/28 09:18:20  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.12  2011/12/23 14:57:11  bigatti
// -- changed and --> && , or --> ||
//
// Revision 1.11  2011/03/11 11:06:14  bigatti
// -- changed size_t --> long
// -- changed size --> len
//
// Revision 1.10  2009/10/28 15:42:02  bigatti
// -- minor fixes for compilation warnings (int --> size_t)
//
// Revision 1.9  2009/10/28 10:12:22  bigatti
// -- updates by Eduardo Saenz-de-Cabezon
//
// Revision 1.8  2009/07/02 16:35:14  abbott
// Minor change: improved layout.  Maybe also something to keep the compiler quiet.
//
// Revision 1.7  2009/06/10 08:19:16  bigatti
// -- additions at the CoCoASchool
//
// Revision 1.6  2009/01/08 13:11:08  bigatti
// -- fixed MVTN1 bug (myOutPPsPtr is now initialized in myReadArgs)
//
// Revision 1.5  2008/07/17 16:35:55  bigatti
// -- use I instead on J in HasHomology
//
// Revision 1.4  2008/07/16 10:00:52  bigatti
// -- added: Betti diagram (by Eduardo Saenz-de-Cabezon)
//
// Revision 1.3  2008/07/04 12:12:56  bigatti
// -- added general MayerVietorisTree
//
// Revision 1.2  2008/07/04 09:11:04  bigatti
// -- new PPVector class
//
// Revision 1.1  2008/05/27 16:21:26  bigatti
// -- first import
//
