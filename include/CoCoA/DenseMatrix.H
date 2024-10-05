#ifndef CoCoA_DenseMatrix_H
#define CoCoA_DenseMatrix_H

//   Copyright (c)  2005-2009,2011  John Abbott,  Anna M. Bigatti

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
//using std::vector;


namespace CoCoA
{

  class BigInt;           // fwd decl -- defined in BigInt.H
  class BigRat;           // fwd decl -- defined in BigRat.H
  class ConstMatrixView;  // fwd decl -- defined in matrix.H
  class RingElem;         // fwd decl -- defined in ring.H
  class matrix;           // fwd decl -- defined in matrix.H
  class ring;             // fwd decl -- defined in ring.H

  matrix NewDenseMat(const ring& R, long NumRows, long NumCols);
  matrix NewDenseMat(const ring& R, const std::vector< std::vector<long> >& M);
  matrix NewDenseMat(const ring& R, const std::vector< std::vector<BigInt> >& M);
  matrix NewDenseMat(const ring& R, const std::vector< std::vector<BigRat> >& M);
  matrix NewDenseMat(const ring& R, const std::vector< std::vector<RingElem> >& M);
  matrix NewDenseMatTranspose(const ring& R, const std::vector< std::vector<long> >& M);
  matrix NewDenseMatTranspose(const ring& R, const std::vector< std::vector<BigInt> >& M);
  matrix NewDenseMatTranspose(const ring& R, const std::vector< std::vector<BigRat> >& M);
  matrix NewDenseMatTranspose(const ring& R, const std::vector< std::vector<RingElem> >& M);
  matrix NewDenseMat(const ConstMatrixView& M);

} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/DenseMatrix.H,v 1.17 2022/02/18 14:11:45 abbott Exp $
// $Log: DenseMatrix.H,v $
// Revision 1.17  2022/02/18 14:11:45  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.16  2021/01/07 14:11:26  abbott
// Summary: Corrected copyright
//
// Revision 1.15  2020/02/28 08:55:14  abbott
// Summary: Added NewDenseMatTranspose
//
// Revision 1.14  2015/11/30 21:53:55  abbott
// Summary: Major update to matrices for orderings (not yet complete, some tests fail)
//
// Revision 1.13  2011/11/09 13:38:31  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.12  2011/11/07 10:28:25  bigatti
// -- forward --> fwd
//
// Revision 1.11  2011/10/07 12:27:21  bigatti
// -- added fwd decl of BigRat
//
// Revision 1.10  2011/10/05 09:23:46  abbott
// Added pseudo ctor from vector of vector of BigRat.
//
// Revision 1.9  2011/10/04 15:36:47  abbott
// Added new pseudo ctors for DenseMatrix (from vector of vector)
//
// Revision 1.8  2011/03/21 13:19:35  bigatti
// -- removed useless fwd declaration
//
// Revision 1.7  2011/03/08 17:54:54  bigatti
// -- changed: args for rows and cols are now  long  instead of  MachineInt
//
// Revision 1.6  2011/03/04 16:16:53  bigatti
// -- changed: args of type MachineInt instead of size_t
//
// Revision 1.5  2009/09/24 13:39:12  abbott
// Removed unnecessary forward declaration.
// Added missing logging keyword.
//
//

#endif
