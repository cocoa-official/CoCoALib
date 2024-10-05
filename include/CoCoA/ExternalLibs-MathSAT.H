#ifndef CoCoA_mathsat_H
#define CoCoA_mathsat_H

//   Copyright (c)  2017 John Abbott, Anna M. Bigatti
//   Authors: 2017 Anna M. Bigatti

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

#include "CoCoA/PREPROCESSOR_DEFNS.H"

#ifdef CoCoA_WITH_MATHSAT

#include "mathsat.h"
#include "CoCoA/matrix.H"
#include "CoCoA/ring.H"

namespace CoCoA
{
  class matrix; // forward declaration -- defined in matrix.H 
  class ConstMatrixView; // forward declaration -- defined in MatrixView.H 


namespace MathSAT
{

  enum RelOp {eq0, neq0, leq0, lt0};

  RelOp ToRelOp(const std::string& s);

  class env
  {
  public:
    env();
    ~env();

    void myInitializeOrCheck(long NumIndets); // might throw

    void myAddEq0(ConstRefRingElem Lin);
    void myAddNeq0(ConstRefRingElem Lin);
    void myAddLeq0(ConstRefRingElem Lin);
    void myAddLt0(ConstRefRingElem Lin);
    
    void myAddEq0(ConstMatrixView M);
    void myAddNeq0(ConstMatrixView M);
    void myAddLeq0(ConstMatrixView M);
    void myAddLt0(ConstMatrixView M);

    matrix myLinSolve() const;
    bool IamSatisfiable() const;

  private: // data members
    long myNumIndetsValue;
    //    std::vector<msat_term> myIndets;
    std::vector<std::string> myIndetNames;
    mutable msat_env myEnvValue; // need "mutable"?
  };

  void AddConstraint(env& E, RelOp, ConstMatrixView M);
  void AddConstraint(env& E, RelOp, ConstRefRingElem Lin);
  void AddEq0(env& E, ConstRefRingElem Lin);
  void AddNeq0(env& E, ConstRefRingElem Lin);
  void AddLeq0(env& E, ConstRefRingElem Lin);
  void AddLt0(env& E, ConstRefRingElem Lin);
  
  void AddEq0(env& E, ConstMatrixView M);
  void AddNeq0(env& E, ConstMatrixView M);
  void AddLeq0(env& E, ConstMatrixView M);
  void AddLt0(env& E, ConstMatrixView M);
  
  matrix LinSolve(const env& E);




//   // raw interface with msat_env
//   void AddLinEq0(msat_env& MSE, ConstRefRingElem Lin);
//   void AddLinNeq0(msat_env& MSE, ConstRefRingElem Lin);
//   void AddLinLeq0(msat_env& MSE, ConstRefRingElem Lin);
//   void AddLinLt0(msat_env& MSE, ConstRefRingElem Lin);

//   void AddLinEq0(msat_env& MSE, ConstMatrixView M);
//   void AddLinNeq0(msat_env& MSE, ConstMatrixView M);
//   void AddLinLeq0(msat_env& MSE, ConstMatrixView M);
//   void AddLinLt0(msat_env& MSE, ConstMatrixView M);

//   matrix LinSolve(const msat_env& MSE, long NumIndets);
//   bool IsSatisfiable(const msat_env& MSE);
 
} // end of namespace MathSAT
} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/ExternalLibs-MathSAT.H,v 1.10 2022/02/18 14:11:45 abbott Exp $
// $Log: ExternalLibs-MathSAT.H,v $
// Revision 1.10  2022/02/18 14:11:45  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.9  2018/06/25 10:09:24  bigatti
// -- changed: myInitialize into myInitializeOrCheck
//
// Revision 1.8  2018/05/28 08:47:25  bigatti
// -- added myIndetNames
// -- added myInitialize
//
// Revision 1.7  2017/08/08 14:25:49  bigatti
// -- updated for SC2
//
// Revision 1.6  2017/07/24 14:52:55  bigatti
// -- cleaned up design for MathSAT wrapper class
//
// Revision 1.5  2017/07/14 09:31:15  bigatti
// -- new class MathSAT::env.  Consequent changes
//
// Revision 1.4  2017/07/12 16:44:53  bigatti
// -- developed experimental code for MathSat and moved it from example
//    to ExternalLib-MathSAT.[CH]
//
// Revision 1.3  2017/04/27 16:18:50  bigatti
// -- changed copyright
//
// Revision 1.2  2017/03/13 12:17:08  abbott
// Summary: Include PREPROCESSOR_DEFNS.H (in header files which need the info)
//
// Revision 1.1  2017/02/24 08:20:13  bigatti
// -- first import
//

#endif // CoCoA_WITH_MATHSAT
#endif
