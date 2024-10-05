#ifndef CoCoA_MatrixFp_H
#define CoCoA_MatrixFp_H

//   Copyright (c)  2016  John Abbott,  Anna M. Bigatti

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

#include "CoCoA/config.H"
#include "CoCoA/SmallFpImpl.H"

#include <iosfwd>
//using std::ostream;
#include <vector>
//using std::vector;

namespace CoCoA
{

  class matrix; // fwd decl
  class ConstMatrixView; // fwd decl


class MatrixFp
{
public:
  MatrixFp(const SmallFpImpl& ModP, int nrows, int ncols);
  explicit MatrixFp(const ConstMatrixView& M);
//  MatrixFp(const SmallFpImpl& Modp, const MatrixFpNonRed& M);
  friend long NumRows(const MatrixFp& M);
  friend long NumCols(const MatrixFp& M);
  friend const SmallFpImpl& ModPArith(const MatrixFp& M);
  friend std::ostream& operator<<(std::ostream& out, const MatrixFp& M);
  SmallFpImpl::value& operator()(long r, long c);
  const SmallFpImpl::value& operator()(long r, long c) const;
  std::vector<SmallFpImpl::value>& operator[](long r);
  const std::vector<SmallFpImpl::value>& operator[](long r) const;
//  friend void SetEntry(MatrixFp& M, long r, long c, SmallFpImpl::value x);
private: // data members
  long myRows;
  long myCols;
  const SmallFpImpl& myArith;
  std::vector< std::vector< SmallFpImpl::value > > myM;
};

// accessor fns
  inline long NumRows(const MatrixFp& M) { return M.myRows; }
  inline long NumCols(const MatrixFp& M) { return M.myCols; }
  inline const SmallFpImpl& ModPArith(const MatrixFp& M) { return M.myArith; }
  inline SmallFpImpl::value& MatrixFp::operator()(long r, long c) { return myM[r][c]; }
  inline const SmallFpImpl::value& MatrixFp::operator()(long r, long c) const { return myM[r][c]; }

  inline std::vector<SmallFpImpl::value>& MatrixFp::operator[](long r) { return myM[r]; }
  inline const std::vector<SmallFpImpl::value>& MatrixFp::operator[](long r) const  { return myM[r]; }

  std::ostream& operator<<(std::ostream& out, const MatrixFp& M);

  matrix LinSolveFp(ConstMatrixView M, ConstMatrixView rhs);
  matrix InverseFp(ConstMatrixView M);
  matrix LinKerFp(ConstMatrixView M_orig);
  
  class LinDepFp
  {
  public:
    LinDepFp(const SmallFpImpl& ModP, long dim);
    bool myAppendVec(std::vector<SmallFpImpl::value>& v);
//    bool IsLinDep() const { return };
    const std::vector<SmallFpImpl::value>& myLinReln() { /*CoCoA_ASSERT(???);*/ return myLinRelnVec;}

  private:
    const SmallFpImpl& myArith;
    long myDim;
    int myNumVecs;
    std::vector< std::vector<SmallFpImpl::value> > myM;
    std::vector<SmallFpImpl::NonRedValue> myRow;
    std::vector<int> myColIndices;
    std::vector<SmallFpImpl::value> myLinRelnVec;
  };

  long/*SmallFpImpl::value*/ det(MatrixFp& M);


  class RowRednFp
  {
  public:
    RowRednFp(const SmallFpImpl& ModP, long dim);
    void myAppendVec(std::vector<SmallFpImpl::value>& v);
//    bool IsLinDep() const { return };
//    const std::vector<SmallFpImpl::value>& myLinReln() { /*CoCoA_ASSERT(???);*/ return myLinRelnVec;}
    long myRank() const;
    SmallFpImpl::value myDet() const;
  private:
    const SmallFpImpl& myArith;
    long myDim;
    int myNumVecs;
    SmallFpImpl::value myDiagProd;
    std::vector< std::vector<SmallFpImpl::value> > myM;
    std::vector<SmallFpImpl::NonRedValue> myRow;
//    std::vector<std::vector<RingElem> > myRowRepr;
    std::vector<int> myColIndices;
//    std::vector<SmallFpImpl::value> myLinRelnVec;

  };

} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/MatrixFp.H,v 1.6 2023/11/11 20:27:14 abbott Exp $
// $Log: MatrixFp.H,v $
// Revision 1.6  2023/11/11 20:27:14  abbott
// Summary: Added printing for MatrixFp
//
// Revision 1.5  2023/10/25 19:21:36  abbott
// Summary: Added InverseFp
//
// Revision 1.4  2022/02/18 14:11:46  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.3  2021/01/07 14:16:10  abbott
// Summary: Corrected copyright
//
// Revision 1.2  2018/02/15 14:40:01  abbott
// Summary: update
//
// Revision 1.1  2018/01/17 12:41:26  abbott
// Summary: new headers
//
//

#endif
