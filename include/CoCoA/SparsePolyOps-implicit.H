#ifndef CoCoA_TmpImplicit_H
#define CoCoA_TmpImplicit_H

//   Copyright (c)  2014  John Abbott, and Anna M. Bigatti

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

#include "CoCoA/ring.H"

#include <vector>

namespace CoCoA
{

  RingElem ImplicitDirect(const std::vector<RingElem>& ParamDescr);
  RingElem ImplicitDirectLPP(const std::vector<RingElem>& ParamDescr);
  RingElem ImplicitDirectLPP2(const std::vector<RingElem>& ParamDescr);
  RingElem ImplicitDirectOrd2(const std::vector<RingElem>& ParamDescr);
  RingElem ImplicitDirectWithCond(const std::vector<RingElem>& ParamDescr, const std::vector<RingElem>& relations);
  RingElem ImplicitDirectWithCondLPP(const std::vector<RingElem>& ParamDescr, const std::vector<RingElem>& relations);

  RingElem ImplicitByPoints(const std::vector<RingElem>& ParamDescr);
  RingElem ImplicitByPoints2(const std::vector<RingElem>& ParamDescr);
  RingElem ImplicitByPoints3(const std::vector<RingElem>& ParamDescr);

  RingElem SliceCore(const std::vector<RingElem>& ParameDescr,
                     long RecDepth,
                     const std::string& FinalCalls);
  RingElem SliceCoreQQ(const std::vector<RingElem>& ParamDescr,
                       long RecDepth,
                       const std::string& FinalCall);
  

} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/SparsePolyOps-implicit.H,v 1.3 2022/02/18 14:11:49 abbott Exp $
// $Log: SparsePolyOps-implicit.H,v $
// Revision 1.3  2022/02/18 14:11:49  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.2  2021/01/07 14:44:46  abbott
// Summary: Corrected copyright
//
// Revision 1.1  2018/04/06 15:02:31  bigatti
// -- renamed TmpImplicit.H
//
// Revision 1.9  2016/01/26 13:49:14  bigatti
// -- added ImplicitDirectOrd2
//
// Revision 1.8  2015/03/09 12:11:30  bigatti
// -- added SliceCoreQQ
//
// Revision 1.7  2014/12/10 12:07:56  bigatti
// -- simplified SliceCore syntax
//
// Revision 1.6  2014/11/28 15:16:11  bigatti
// -- added arg to SliceCore for algorithm of final calls
//
// Revision 1.5  2014/11/27 11:29:07  abbott
// Summary: Added ImplicitDirectWithCondLPP
// Author: JAA
//
// Revision 1.4  2014/11/17 10:34:53  abbott
// Summary: Added ImplicitDirectLPP2 and ImplicitDirectWithCond
// Author: JAA
//
// Revision 1.3  2014/11/13 15:59:13  abbott
// Summary: Added ImplicitDirectLPP & ImplicitByPoints2
// Author: JAA
//
// Revision 1.2  2014/11/11 09:17:18  bigatti
// -- added SliceCore
//
// Revision 1.1  2014/11/04 18:10:12  abbott
// Summary: new tmp code for implicitization
// Author: JAA
//
//

#endif
