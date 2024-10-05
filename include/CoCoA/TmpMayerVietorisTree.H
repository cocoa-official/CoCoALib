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

// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/TmpMayerVietorisTree.H,v 1.12 2024/03/28 22:09:41 abbott Exp $
// $Log: TmpMayerVietorisTree.H,v $
// Revision 1.12  2024/03/28 22:09:41  abbott
// Summary: Corrected include guards (redmine 1811)
//
// Revision 1.11  2022/02/18 14:11:50  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.10  2014/09/04 12:17:05  bigatti
// -- handy interface with ideal (should it be here or in MonomialIdeal?)
//
// Revision 1.9  2014/08/16 13:23:24  abbott
// Summary: Cleaned up include directives; added fwd decl; exposed a print operator
// Author: JAA
//
// Revision 1.8  2011/08/14 15:52:17  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.7  2009/10/28 15:44:05  bigatti
// -- changed MVTProjectiveDimension into MVTProjDim
//
// Revision 1.6  2009/10/28 10:12:22  bigatti
// -- updates by Eduardo Saenz-de-Cabezon
//
// Revision 1.5  2009/06/10 08:19:16  bigatti
// -- additions at the CoCoASchool
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

#endif /* TmpMVT_H */

