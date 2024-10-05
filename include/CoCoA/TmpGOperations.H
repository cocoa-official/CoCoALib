#ifndef CoCoA_GOperations_H
#define CoCoA_GOperations_H

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

#include "CoCoA/TmpGTypes.H"
#include "CoCoA/TmpGReductor.H"

//#include <vector>
// using vector;

// Operations depending on Groebner Bases


namespace CoCoA
{

  class ConstRefPPMonoidElem; // forward declaration -- defined in PPMonoid.H
  class RingElemAlias;        // forward declaration -- defined in ring.H
  typedef const RingElemAlias& ConstRefRingElem;     // forward declaration -- defined in ring.H
  class ModuleElem;           // forward declaration -- defined in PPMonoid.H
  class RingElem;             // forward declaration -- defined in ring.H

  void ComputeGBasis(PolyList& outGB, PolyList& outMinGens, const PolyList& , const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());
  void ComputeGBasisTrunc(PolyList& outGB, PolyList& outMinGens, long& TruncDeg, const PolyList& , const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());
  void ComputeGBasisSelfSatCore(PolyList& outGB, const PolyList& inGens, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());
  void ComputeGBasisRealSolve(PolyList& outGB, const PolyList& inGens, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());
  //  void ComputeSATMixGBasis(PolyList& outGB, const PolyList& inGens);
  
  void ComputeGBasis(VectorList& outGB, VectorList& outMinGens, const VectorList& inGens, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());
  void ComputeGBasisFrameWork(PolyList&, const PolyList&);

// The second argument is the list of indexes of the
// indeterminates that will be eliminated
void ComputeElim(VectorList&,const VectorList&,ConstRefPPMonoidElem, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());///< Not yet implemented
void ComputeElim(PolyList&,const PolyList&,ConstRefPPMonoidElem, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());

  void ComputeSyz(VectorList&, const FreeModule&, const VectorList&, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());
  void ComputeSyz(VectorList&, const FreeModule&, const PolyList&, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());

void ComputeIntersection(VectorList&,const VectorList&,const VectorList&, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());
void ComputeIntersection(PolyList&,const PolyList&,const PolyList&, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());

void ComputeColonByPrincipal(PolyList&, const VectorList&,const VectorList&, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());
void ComputeColonByPrincipal(PolyList&, const PolyList&, ConstRefRingElem f, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());
void ComputeColonByPrincipal(VectorList&, const VectorList&,const PolyList&, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());

void ComputeColon(PolyList&,const VectorList&,const VectorList&, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());
void ComputeColon(PolyList&,const PolyList&,const PolyList&, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit() );
void ComputeColon(VectorList&,const VectorList&,const PolyList&, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());

void ComputeSaturationByPrincipal(VectorList&,const VectorList&,const PolyList&, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());
void ComputeSaturationByPrincipal(PolyList&, const PolyList&, ConstRefRingElem, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());

void ComputeSaturation(PolyList&,const PolyList&,const PolyList&, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());
void ComputeSaturation(VectorList&,const VectorList&,const PolyList&, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());

void ComputeHomogenization(VectorList&,const VectorList&,const PolyList&, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());///< Not yet implemented
void ComputeHomogenization(PolyList&,const PolyList&,const PolyList&, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());///< Not yet implemented

bool RadicalMembership(const PolyList&,ConstRefRingElem, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());

void ComputeLT(VectorList&,const VectorList&, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());///< Not yet implemented
void ComputeLT(PolyList&,const PolyList&, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());

  RingElem ComputeElimFirst(const PolyList& inPL, ConstRefPPMonoidElem inds, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());

}// end namespace cocoa


// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/TmpGOperations.H,v 1.27 2024/05/28 12:53:40 bigatti Exp $
// $Log: TmpGOperations.H,v $
// Revision 1.27  2024/05/28 12:53:40  bigatti
// Summary: added GBasisTrunc, myGBasisTrunc, ComputeGBasisTrunc
//
// Revision 1.26  2023/07/04 09:18:01  abbott
// Summary: Changed TimeOut to Timeout
//
// Revision 1.25  2022/02/18 14:11:50  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.24  2022/02/14 14:09:36  bigatti
// Summary: added date for original author
//
// Revision 1.23  2022/02/11 09:49:34  abbott
// Summary: Updated copyright notices (redmine 855)
//
// Revision 1.22  2018/06/27 12:15:18  abbott
// Summary: Renamed RealSolveCore to RealSolve
//
// Revision 1.21  2018/06/27 08:50:39  abbott
// Summary: Revised to work with new CpuTimeLimit
//
// Revision 1.20  2018/05/25 09:24:46  abbott
// Summary: Major redesign of CpuTimeLimit (many consequences)
//
// Revision 1.19  2017/11/29 17:41:48  bigatti
// -- added GBasisRealSolveCore
//
// Revision 1.18  2017/11/24 17:46:40  bigatti
// -- renamed GBasisSelfSat --> GBasisSelfSatCore
// -- added GBasisSelfSat in cpkg5
//
// Revision 1.17  2017/11/23 12:37:20  bigatti
// -- added GBasisSelfSat
//
// Revision 1.16  2017/04/18 09:26:37  bigatti
// -- removed StatLevel argument (now using VerbosityLevel)
//
// Revision 1.15  2015/12/04 15:22:24  bigatti
// -- renamed ComputeSSaturation into ComputeSaturation
//
// Revision 1.14  2015/05/20 13:06:41  bigatti
// -- renamed ComputeCColon --> ComputeColon
//
// Revision 1.13  2015/03/04 10:17:34  bigatti
// -- added (..)ElimFirst
//
// Revision 1.12  2014/03/26 15:23:44  bigatti
// -- added MinGens for submodules
//
// Revision 1.11  2013/06/12 08:54:07  bigatti
// -- added computation of MinGens (in ComputeGBasis)
// -- changed some "the" in "out"/"in" in argument names
//
// Revision 1.10  2013/02/21 17:32:28  bigatti
// -- changed syntax for ComputeSyz
//
// Revision 1.9  2013/01/31 11:43:22  bigatti
// -- added Stats field to ComputeXXXGBasis for returning statistics
//
// Revision 1.8  2012/10/24 12:07:30  abbott
// Revised fwd decl of ConstRefRingElem.
//
// Revision 1.7  2012/10/03 12:22:01  bigatti
// -- added homog
//
// Revision 1.6  2011/12/07 15:54:34  bigatti
// -- renamed ambiguous "operator<" and hidden into anonymous namespace
// -- renamed ambiguous "operator==" into AreEqualLPPs (used by Saturation)
//
// Revision 1.5  2011/12/05 16:32:11  bigatti
// -- fixed bug about saturation (by non-principal ideal)
//
// Revision 1.4  2008/09/19 11:34:16  bigatti
// -- new mechanism for passing verbosity level (or StatLevel)
//    [only partially tested]
//
// Revision 1.3  2007/11/09 10:45:52  bigatti
// -- [caboara] preparation for self-saturating algorithm
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
// Revision 1.10  2007/03/08 18:42:06  cocoa
// Cleaned up whitespace.
//
// Revision 1.9  2007/03/07 22:24:03  bigatti
// -- reintroduced TmpGTypes.H (waiting for a better solution)
//
// Revision 1.8  2007/03/07 17:04:31  cocoa
// -- several changes by M.Caboara: more operations on ideals,
//    exception cleaner, coding conventions, WSugar, dynamic
//
// Revision 1.7  2006/12/21 13:48:33  cocoa
// Made all increment/decrement calls prefix (except where the must be postfix).
//
// Revision 1.6  2006/12/06 17:11:16  cocoa
// -- removed #include "config.H"
//
// Revision 1.5  2006/11/24 17:17:58  cocoa
// -- reorganized includes of header files
//
// Revision 1.4  2006/11/09 17:37:50  cocoa
// -- just some doxygen comments
//
// Revision 1.3  2006/10/06 16:46:17  cocoa
// -- syzygies for non-homogenous polynomials (Max)
// -- wip: evolution of Groebner Framework (Max)
//
// Revision 1.2  2006/10/06 14:04:16  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.1.1.1  2006/05/30 11:39:36  cocoa
// Imported files
//
// Revision 1.11  2006/05/16 08:59:16  cocoa
// -- added function for interactive Groebner
//
// Revision 1.10  2006/05/11 16:00:22  cocoa
// -- fixed spelling of "homogenize"
//
// Revision 1.9  2006/04/27 14:01:11  cocoa
// -- tidied up include files (using GTypes.H)
//
// Revision 1.8  2006/04/21 16:47:06  cocoa
// -- new syntax for ComputeGBasis by Max
//
// Revision 1.7  2006/04/11 16:22:40  cocoa
// -- added: Elim, LT
//

#endif
