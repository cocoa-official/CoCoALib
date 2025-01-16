#ifndef CoCoA_TmpMVT_H
#define CoCoA_TmpMVT_H

//   Copyright (c)  2008-2009  Anna Bigatti and Eduardo Sáenz-de-Cabezón

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


// Source code for computations using Mayer-Vietoris trees
// Includes computation of Irreducible Components using N-1st Koszul homology

#include "CoCoA/BigInt.H"  // for position_t

#include <iosfwd>
// using std::ostream;
#include <list>
/// using list
#include <map>
// using map;
#include <vector>
// using vector;

/**********************************************************/
namespace CoCoA
{ 
  class PPVector;  // fwd decl -- defined in TmpPPVector.H
  class PPMonoidElem; // fwd decl -- defined in PPMonoid.H
  class RingElem; // fwd decl -- defined in ring.H
  class ideal; // fwd decl -- defined in ideal.H

  typedef BigInt position_t;

//**************************************************************//
//              LISTS OF DIMS                                   //
//**************************************************************//

//This is the structure to encode a list of the dimensions and positions in which a multidegree appears in MVT(I);
//each element of the list has the form (d,L_d) where d is a dimension, L_d is the list of positions in which the corresponding multidegree appears in relevant nodes in the tree.

  typedef std::map<int, std::list<position_t> > ListOfDims;

//***********************************************************//
//  CONTAINER FOR THE REPEATED AND NOT REPATED MULTIDEGREES  //
//***********************************************************//

//To store more efficiently lists of multidegrees, and the dimensions in which they apepar in Mayer-Vietoris trees, we use maps. Here all the necessary types and functions are defined

  typedef std::map<PPMonoidElem, ListOfDims> MultidegreeMap;

  std::ostream& operator<<(std::ostream& out, const MultidegreeMap& L);



//**************************************************************//
//              BETTI PSEUDO-DIAGRAM                            //
//**************************************************************//

//This is the structure to encode a Betti diagram that gives an "approximation" of the correct one

  typedef std::map<int,std::vector<int> > BettiPseudoDiagram;



  std::vector<RingElem> MayerVietorisTreeN1(const ideal& I);

  void MayerVietorisTreeN1(PPVector& betti, const PPVector& I);
  void MayerVietorisTreeN1_gen(PPVector& betti, const PPVector& I, const PPMonoidElem&, int);
  void MayerVietorisTree(MultidegreeMap& betti, const PPVector& I);

  void MayerVietorisTreeExtremal(MultidegreeMap& output_list, std::list<PPMonoidElem>& extremal, const PPVector& I); //for the CoCoA school 2009


  void ReducedMayerVietorisTree(MultidegreeMap& betti, const PPVector& I);
  void Bettis(std::vector<int>& bettis, const MultidegreeMap& the_list);
  void ReduceMVTree(MultidegreeMap& the_map, MultidegreeMap& undecided_map);
  int ResSize(const std::vector<int>& R);
  void GradedBettis(const MultidegreeMap& decided,BettiPseudoDiagram& graded_bettis, int& reg, int& pd);
  void PrintBettiDiagram(const BettiPseudoDiagram& graded_bettis);
  void PseudoBettiDiagram(const BettiPseudoDiagram&, const BettiPseudoDiagram&);
  void PrintMultidegreeMap(const MultidegreeMap& myMap);
  void PrintMultidegreeMapDim(const MultidegreeMap& myMap);
  int regularity(const MultidegreeMap& decided);
  void irreducible(std::vector<PPVector>&, const PPVector&, const PPMonoidElem&);

//INTERFACE FUNCTIONS

void MVTPrint(const PPVector& );
void MVTPrintDim(const PPVector& );
void MVTPseudoBettiDiagram(const PPVector& );
int MVTRegularityLowerBound(const PPVector& );
int MVTRegularityUpperBound(const PPVector& );
int MVTRegularity(const PPVector& );
//void MVTRegularity(int& ,const PPVector& );
int MVTProjDimLowerBound(const PPVector& );
int MVTProjDimUpperBound(const PPVector& );
int MVTProjDim(const PPVector& );

}

#endif /* TmpMVT_H */

