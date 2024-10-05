#ifndef CoCoA_MatrixForOrdering_H
#define CoCoA_MatrixForOrdering_H

//   Copyright (c)  2008,2015  John Abbott,  Anna M. Bigatti

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

#include <vector>
// using std::vector;

namespace CoCoA
{

  class ConstMatrixView; // fwd decl -- defined in matrix.H
  class ConstMatrix;     // fwd decl -- defined in matrix.H
  class MachineInt;      // fwd decl -- defined in MachineInt.H
  class matrix;          // fwd decl -- defined in matrix.H

  bool IsTermOrdering(const ConstMatrixView& M);
  bool IsPositiveGrading(const ConstMatrixView& M);  ///< is M a positive grading?
  //  bool IsNonNegGrading(const ConstMatrixView& M);    // removed 2019-09

  ConstMatrix LexMat(const MachineInt& n);
  ConstMatrix XelMat(const MachineInt& n);
  ConstMatrix RevLexMat(const MachineInt& n);
  ConstMatrix StdDegRevLexMat(const MachineInt& n);
  ConstMatrix StdDegLexMat(const MachineInt& n);

  matrix MakeTermOrdMat(ConstMatrixView M);
  matrix MakeTermOrdMat(ConstMatrixView M, const MachineInt& GrDim);

  matrix ElimMat(const std::vector<long>& IndetsToElim, const MachineInt& NumIndets);
  matrix ElimMat(const std::vector<long>& IndetsToElim, const ConstMatrixView& GradingM);
  matrix ElimHomogMat(const std::vector<long>& IndetsToElim, const ConstMatrixView& GradingM);

} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/MatrixForOrdering.H,v 1.26 2022/02/18 14:11:46 abbott Exp $
// $Log: MatrixForOrdering.H,v $
// Revision 1.26  2022/02/18 14:11:46  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.25  2022/02/04 21:41:06  abbott
// Summary: Changed name MakeTermOrd to MakeTermOrdMat (redine 854)
//
// Revision 1.24  2021/01/07 14:15:38  abbott
// Summary: Corrected copyright
//
// Revision 1.23  2019/09/25 14:25:05  bigatti
// -- just some resorting of the code for IsPositiveGrading
// -- IsNonNegGrading is now hidden in anonymous namespace
//
// Revision 1.22  2018/06/15 08:46:48  abbott
// Summary: Added IsNonNegGrading
//
// Revision 1.21  2016/09/22 15:33:37  bigatti
// -- renamed HomogElimMat into ElimHomogMat
// -- improved readability for ElimHomogMat/ElimMat (removed auxiliary functions)
//
// Revision 1.20  2015/12/08 14:05:11  abbott
// Summary: Renamed NewMatCompleteOrd to MakeTermOrd
//
// Revision 1.19  2015/12/01 16:53:15  abbott
// Summary: Added new MatCompleteOrd with 2nd arg
//
// Revision 1.18  2015/12/01 15:57:04  abbott
// Summary: Commented out NewPositiveMat, NewMatMinimize, NewIntegerOrdMat
//
// Revision 1.17  2015/12/01 13:34:44  abbott
// Summary: Changed arg order in ElimMat and HomogElimMat; doc is out-of-date!!
//
// Revision 1.16  2015/11/30 21:53:55  abbott
// Summary: Major update to matrices for orderings (not yet complete, some tests fail)
//
// Revision 1.15  2015/04/13 15:36:55  abbott
// Summary: Added decls of "MatXel", "MatStdDegRevLex", "MatStdDegLex", "MatRevLex"
// Author: JAA
//
// Revision 1.14  2014/07/15 16:53:19  abbott
// Summary: Corrected name of include guard
// Author: JAA
//
// Revision 1.13  2013/05/27 13:10:41  abbott
// Removed unnecessary include of config.H
//
// Revision 1.12  2013/02/14 17:34:35  bigatti
// -- cleaned up code for elimination matrices
//
// Revision 1.11  2012/03/30 17:28:09  bigatti
// -- added NewIntegerOrdMat
// -- accepting and returning matrices over QQ
//
// Revision 1.10  2012/02/08 17:21:07  bigatti
// -- changed: Z,Q -> ZZ,QQ
// -- code reorganization (more work to do...)
//
// Revision 1.9  2011/12/05 16:54:40  bigatti
// -- changed: MachineInteger --> MachineInt (just in comment)
//
// Revision 1.8  2011/05/26 11:57:28  bigatti
// -- added IsPositiveGrading with one arg
//
// Revision 1.7  2011/04/26 10:10:31  bigatti
// -- added NewMatCompleteOrd
//
// Revision 1.6  2011/03/23 17:29:54  bigatti
// -- added NewDenseMatStdDegLex
//
// Revision 1.5  2011/03/21 07:50:51  bigatti
// -- added NewDenseMatXel, NewDenseMatStdDegRevLex
//
// Revision 1.4  2011/03/08 16:10:16  abbott
// Changed size_t into long.
//
// Revision 1.3  2009/09/22 13:35:55  bigatti
// -- following coding conventions in function names Matrix --> Mat
// -- forced all matrices to be over RingZ
//
// Revision 1.2  2008/12/17 12:11:52  abbott
// Changed type from long to MachineInt in operations which use a machine integer
// in place of a RingElem.  The change is "superficial" but affects many files.
//
// Revision 1.1  2008/04/21 11:23:11  abbott
// Separated functions dealing with matrices and PPOrderings into a new file.
// Added matrix norms, and completed adjoint.
//
//

#endif
