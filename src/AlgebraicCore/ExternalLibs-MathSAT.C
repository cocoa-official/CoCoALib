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


// Source code for MathSAT integration

#include "CoCoA/ExternalLibs-MathSAT.H"

#ifdef CoCoA_WITH_MATHSAT

#include "CoCoA/BigRat.H"
#include "CoCoA/error.H"
#include "CoCoA/DenseMatrix.H"
#include "CoCoA/matrix.H"
#include "CoCoA/time.H"
#include "CoCoA/SparsePolyIter.H"
#include "CoCoA/SparsePolyOps-RingElem.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/verbose.H"

// #include "CoCoA/VectorOps.H"  // just for debugging

#include <cstdlib>
using std::atoi;
//#include <map>
//using std::map;
// #include <vector>
// using std::vector;
#include <ostream>
// using std::ostream;
using std::endl;
#include <sstream>  // for 
#include <string>

namespace CoCoA
{

  namespace MathSAT  // "mathsat" is its own namespace
  {


    RelOp ToRelOp(const std::string& s)
    {
      if (s == "eq0") return eq0;
      if (s == "neq0") return neq0;
      if (s == "leq0") return leq0;
      if (s == "lt0") return lt0;
      CoCoA_THROW_ERROR2(ERR::ShouldNeverGetHere, "unknown relation "+s);
      return eq0; // just to keel the compiler quiet
    }

    //---------------------------------------------------------------
    namespace  // conversion functions, implementation at the end of the file
    {      

      /////////// improve these conversion......
//       std::string BRToString(ConstRefRingElem r)
//       {
//         std::ostringstream buffer;
//         buffer << r;
//         return buffer.str();
//       }


      /////////// improve these conversion......
//       msat_term IndetToMathSAT(msat_env MSE, int var)
//       {
//         std::ostringstream buf;
//         buf << "X" << var;
//         std::string name = buf.str();
//         msat_type rat = msat_get_rational_type(MSE);
//         msat_decl d = msat_declare_function(MSE, name.c_str(), rat);
//         msat_term x = msat_make_constant(MSE, d);
//         return x;
//       }


      msat_term IndetToMathSAT(msat_env MSE, const std::string& name)
      {
        msat_type rat = msat_get_rational_type(MSE);
        msat_decl d = msat_declare_function(MSE, name.c_str(), rat);
        msat_term x = msat_make_constant(MSE, d);
        return x;
      }
      
      
      msat_term RingElemToMathSAT(msat_env MSE, ConstRefRingElem L,
                                  const std::vector<std::string>& MSIndetNames)
      {
        msat_term res = msat_make_number(MSE, "0");
        mpq_t q;
        mpq_init(q);
        long v;
        for (SparsePolyIter iter=BeginIter(L); !IsEnded(iter); ++iter)
        {
          mpq_set(q, mpqref(ConvertTo<BigRat>(coeff(iter))));
          msat_term c = msat_make_mpq_number(MSE, q);
          if (IsOne(PP(iter)))  res = msat_make_plus(MSE, res, c);
          else
          {
            if (IsIndet(v, PP(iter)))
            {
              msat_term x = IndetToMathSAT(MSE, MSIndetNames[v]);
              res = msat_make_plus(MSE, res, msat_make_times(MSE, c, x));
            }
            else
              CoCoA_THROW_ERROR1(ERR::ReqIndet);
          }
        }
        return res;
      }
      

      msat_term RowToMathSAT(msat_env MSE, ConstMatrixView M, long i,
                             const std::vector<std::string>& MSIndetNames)
      {
        CoCoA_ASSERT_ALWAYS(i>=0 && i<NumRows(M));
        mpq_t q;
        mpq_init(q);
        long v = NumCols(M);
        --v;
        // change sign to constant term: 1 2 3 ==> 1*x+2*y 3 <=> 0
        mpq_set(q, mpqref(ConvertTo<BigRat>(M(i,v))));
        msat_term res = msat_make_mpq_number(MSE, q);
        for (--v; v>=0; --v)
        {
          mpq_set(q, mpqref(ConvertTo<BigRat>(M(i,v))));
          msat_term c = msat_make_mpq_number(MSE, q);
          msat_term x = IndetToMathSAT(MSE, MSIndetNames[v]);
          res = msat_make_plus(MSE, res, msat_make_times(MSE, c, x));
        }
        return res;
      }


//       msat_term RowToMathSAT(msat_env MSE, ConstMatrixView M, long i,
//                              const std::vector<msat_term>& MSIndet)
//       {
//         CoCoA_ASSERT_ALWAYS(i>=0 && i<NumRows(M));
//         mpq_t q;
//         mpq_init(q);
//         BigRat Q;
//         long v = NumCols(M);
//         --v;
//         // 1 2 3 ==> 1*x+2*y +3 <=> 0
//         IsRational(Q, M(i,v));
//         msat_term res = msat_make_mpq_number(MSE, mpqref(Q));
//         for (--v; v>=0; --v)
//         {
//           IsRational(Q, M(i,v));
//           if (!IsZero(Q))
//             res = msat_make_plus(MSE,
//                                  res,
//                                  msat_make_times(MSE,
//                                                  msat_make_mpq_number(MSE, mpqref(Q)),
//                                                  MSIndet[v]));
//         }
//         mpq_clear(q);
//         return res;
//       }


//       struct hash_eq
//       {
//         size_t operator()(msat_term t) const { return msat_term_id(t); }
//         bool operator()(msat_term a, msat_term b) const
//         { return msat_term_id(a) == msat_term_id(b); }
//       };

//       // std::unordered_map<msat_term, long, hash_eq, hash_eq> term_to_index;


//       map<string, long> MSIndetToCoCoAIndex(long NIndets)
//       {
//         map<string, long> res;
//         for (int i=0; i < NIndets; ++i)
//         {
//           std::ostringstream buf;
//           buf << "X" << var;
//           res[buf.str()] = var;
//         }
//         return res;
//       }


      long GetCoeffIndetIndex(msat_term x)
      {
        msat_decl d = msat_term_get_decl(x);
        char *name = msat_decl_get_name(d);
        long var = atoi(name+1);
        msat_free(name);
        return var;
      }
    } // end of anonymous namespace
    //---------------------------------------------------------------




    //---------------------------------------------------------------
    // --- env -------------------

    env::env()
    {
      myNumIndetsValue = 0;
      // create a MathSAT environment
      msat_config cfg = msat_create_config();
      msat_set_option(cfg, "model_generation", "true");
      msat_set_option(cfg, "preprocessor.toplevel_propagation", "false");
      // MathSAT debug ------------------------------------------
//       msat_set_option(cfg, "debug.api_call_trace", "1");
//       msat_set_option(cfg, "debug.api_call_trace_filename",
//                       "/Users/bigatti/Desktop/trace.smt2");
      // MathSAT debug ------------------------------------------
      myEnvValue = msat_create_env(cfg);
      msat_destroy_config(cfg);
    }


    env::~env()
    { msat_destroy_env(myEnvValue); }


    void env::myInitializeOrCheck(long NumIndets)
    {
      if (myNumIndetsValue != 0)
        CoCoA_ASSERT_ALWAYS(myNumIndetsValue == NumIndets);
      else
        myNumIndetsValue = NumIndets;
      for (long i=0; i<NumIndets; ++i)
      {
        std::ostringstream buf;
        buf << "X" << i;
        myIndetNames.push_back(buf.str());
      }
    }


//     void env::myInitializeOrCheck(long NumIndets)
//     {
//       myNumIndetsValue = NumIndets;
//       if (myNumIndetsValue != 0)
//         CoCoA_ASSERT_ALWAYS(myNumIndetsValue == NumIndets);
//       else
//         for (long i=0; i<NumIndets; ++i)
//         {
//           std::ostringstream buf;
//           buf << "X" << i;
//           msat_type rat = msat_get_rational_type(myEnvValue);
//           msat_decl d = msat_declare_function(myEnvValue, buf.str().c_str(), rat);
//           msat_term x = msat_make_constant(myEnvValue, d);
//           myIndets.push_back(x);
//         }
//     }


    //--- access to msat_env: polynomial ----------------------------
    void env::myAddEq0(ConstRefRingElem Lin)
    {
      VerboseLog VERBOSE("MathSAT::env::myAddEq0");
      if (VerbosityLevel()>=99)  VERBOSE(99) << Lin << " == 0" << endl;
      myInitializeOrCheck(NumIndets(owner(Lin)));
      msat_term L = RingElemToMathSAT(myEnvValue, Lin, myIndetNames);
      msat_term Z = msat_make_number(myEnvValue,"0");
      msat_assert_formula(myEnvValue, msat_make_equal(myEnvValue, L,Z));
    }


    void env::myAddNeq0(ConstRefRingElem Lin)
    {
      VerboseLog VERBOSE("MathSAT::env::myAddNeq0");
      if (VerbosityLevel()>=99)  VERBOSE(99) << Lin << " != 0" << endl;
      myInitializeOrCheck(NumIndets(owner(Lin)));
      msat_term L = RingElemToMathSAT(myEnvValue, Lin, myIndetNames);
      msat_term Z = msat_make_number(myEnvValue,"0");
      msat_assert_formula(myEnvValue, msat_make_not(myEnvValue, msat_make_equal(myEnvValue, L,Z)));
    }


    void env::myAddLeq0(ConstRefRingElem Lin)
    {
      VerboseLog VERBOSE("MathSAT::env::myAddLeq0");
      if (VerbosityLevel()>=99)  VERBOSE(99) << Lin << " <= 0" << endl;
      myInitializeOrCheck(NumIndets(owner(Lin)));
      msat_term L = RingElemToMathSAT(myEnvValue, Lin, myIndetNames);
      msat_term Z = msat_make_number(myEnvValue,"0");
      msat_assert_formula(myEnvValue, msat_make_leq(myEnvValue, L, Z));
    }


    void env::myAddLt0(ConstRefRingElem Lin)
    {
      VerboseLog VERBOSE("MathSAT::env::myAddLt0");
      if (VerbosityLevel()>=99)  VERBOSE(99) << Lin << " < 0" << endl;
      myInitializeOrCheck(NumIndets(owner(Lin)));
      msat_term L = RingElemToMathSAT(myEnvValue, Lin, myIndetNames);
      msat_term Z = msat_make_number(myEnvValue,"0");
      msat_assert_formula(myEnvValue, msat_make_not(myEnvValue, msat_make_leq(myEnvValue, Z,L)));
    }

    //--- access to msat_env: matrix ------------------------------
    void env::myAddEq0(ConstMatrixView M)
    {
      VerboseLog VERBOSE("MathSAT::env::myAddEq0");
      if (VerbosityLevel()>=99)  VERBOSE(99) << M << endl;
      myInitializeOrCheck(NumCols(M)-1);
      msat_term Z = msat_make_number(myEnvValue,"0");
      double t0=CpuTime();
      for (long i=0; i<NumRows(M); ++i)
      {
        msat_term L = RowToMathSAT(myEnvValue, M, i, myIndetNames);
        msat_assert_formula(myEnvValue, msat_make_equal(myEnvValue, L,Z));
      }
      VERBOSE(91) << "total for-loop time = " << CpuTime()-t0 << endl;
    }


    void env::myAddNeq0(ConstMatrixView M)
    {
      VerboseLog VERBOSE("MathSAT::env::myAddNeq0");
      myInitializeOrCheck(NumCols(M)-1);
      if (VerbosityLevel()>=99)  VERBOSE(99) << M << endl;
      msat_term Z = msat_make_number(myEnvValue,"0");
      double t0=CpuTime();
      for (long i=0; i<NumRows(M); ++i)
      {
        msat_term L = RowToMathSAT(myEnvValue, M, i, myIndetNames);
        msat_assert_formula(myEnvValue, msat_make_not(myEnvValue, msat_make_equal(myEnvValue, L,Z)));
      }
      VERBOSE(91) << "total for-loop time = " << CpuTime()-t0 << endl;
    }


    void env::myAddLeq0(ConstMatrixView M)
    {
      VerboseLog VERBOSE("MathSAT::env::myAddLeq0");
      myInitializeOrCheck(NumCols(M)-1);
      if (VerbosityLevel()>=99)  VERBOSE(99) << M << endl;
      msat_term Z = msat_make_number(myEnvValue,"0");
      double t0=CpuTime();
      for (long i=0; i<NumRows(M); ++i)
      {
        msat_term L = RowToMathSAT(myEnvValue, M, i, myIndetNames);
        msat_assert_formula(myEnvValue, msat_make_leq(myEnvValue, L, Z));
      }
      VERBOSE(91) << "total for-loop time = " << CpuTime()-t0 << endl;
    }


    void env::myAddLt0(ConstMatrixView M)
    {
      VerboseLog VERBOSE("MathSAT::env::myAddLt0");
      myInitializeOrCheck(NumCols(M)-1);
      if (VerbosityLevel()>=99)  VERBOSE(99) << M << endl;
      msat_term Z = msat_make_number(myEnvValue,"0");
      double t0=CpuTime();
      for (long i=0; i<NumRows(M); ++i)
      {
        msat_term L = RowToMathSAT(myEnvValue, M, i, myIndetNames);
        msat_assert_formula(myEnvValue, msat_make_not(myEnvValue, msat_make_leq(myEnvValue, Z,L)));
      }
      VERBOSE(91) << "total for-loop time = " << CpuTime()-t0 << endl;
    }


    //--- access to msat_env: LinSolve ------------------------------
    matrix env::myLinSolve() const
    {
      VerboseLog VERBOSE("MathSAT::env::myLinSolve");
      double t0=CpuTime();
      msat_result res = msat_solve(myEnvValue);
      VERBOSE(91) << "msat_solve time = " << CpuTime()-t0 << endl;
      if (res == MSAT_UNSAT) return NewDenseMat(RingQQ(), 0,0);
      matrix M = NewDenseMat(RingQQ(), myNumIndetsValue, 1);
      BigRat ans;
      t0=CpuTime();
      msat_model_iterator it = msat_create_model_iterator(myEnvValue);
      VERBOSE(91) << "msat_create_model_iterator time = "<< CpuTime()-t0 <<endl;
      //t0=CpuTime(); // negligible
      while (msat_model_iterator_has_next(it))
      {
        msat_term t, v;
        msat_model_iterator_next(it, &t, &v);
        msat_term_to_number(myEnvValue, v, mpqref(ans));
        long idx = MathSAT::GetCoeffIndetIndex(t);
        SetEntry(M,  idx,0,   ans);
      }
      //VERBOSE(91) << "while total time = " << CpuTime()-t0 << endl;
      return M;
    }


    bool env::IamSatisfiable() const
    { return msat_solve(myEnvValue) == MSAT_SAT; }


    //---------------------------------------------------------------

    void AddConstraint(env& E, RelOp rel, ConstMatrixView M)
    {
      switch (rel)
      {
      case eq0:  E.myAddEq0(M); break;
      case neq0: E.myAddNeq0(M); break;
      case leq0: E.myAddLeq0(M); break;
      case lt0:  E.myAddLt0(M); break;
      }
    }


    void AddConstraint(env& E, RelOp rel, ConstRefRingElem Lin)
    {
      switch (rel)
      {
      case eq0:  E.myAddEq0(Lin); break;
      case neq0: E.myAddNeq0(Lin); break;
      case leq0: E.myAddLeq0(Lin); break;
      case lt0:  E.myAddLt0(Lin); break;
      }
    }


    void AddEq0(env& E, ConstRefRingElem Lin)  { E.myAddEq0(Lin); }
    void AddNeq0(env& E, ConstRefRingElem Lin) { E.myAddNeq0(Lin); }
    void AddLeq0(env& E, ConstRefRingElem Lin) { E.myAddLeq0(Lin); }
    void AddLt0(env& E, ConstRefRingElem Lin)  { E.myAddLt0(Lin); }

    void AddEq0(env& E, ConstMatrixView M)  { E.myAddEq0(M); }
    void AddNeq0(env& E, ConstMatrixView M) { E.myAddNeq0(M); }
    void AddLeq0(env& E, ConstMatrixView M) { E.myAddLeq0(M); }
    void AddLt0(env& E, ConstMatrixView M)  { E.myAddLt0(M); }

    //---------------------------------------------------------------

    matrix LinSolve(const env& E)
    { return E.myLinSolve(); }

    //---------------------------------------------------------------


  } // namespace MathSAT
} // namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/ExternalLibs-MathSAT.C,v 1.19 2024/08/01 10:02:50 bigatti Exp $
// $Log: ExternalLibs-MathSAT.C,v $
// Revision 1.19  2024/08/01 10:02:50  bigatti
// Summary: now using CoCoA_THROW_ERROR1/2  (removing deprecated errors by string)
//
// Revision 1.18  2022/02/18 14:11:54  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.17  2021/01/07 15:07:02  abbott
// Summary: Corrected copyright
//
// Revision 1.16  2020/06/17 15:49:22  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.15  2019/09/16 17:39:25  abbott
// Summary: Corrected indentation
//
// Revision 1.14  2018/06/25 10:09:24  bigatti
// -- changed: myInitialize into myInitializeOrCheck
//
// Revision 1.13  2018/05/28 08:43:09  bigatti
// -- minor improvement to input to msat (vector of indetnames)
//
// Revision 1.12  2018/05/24 15:40:29  bigatti
// -- minor improvement using msat_make_mpq_number
// -- minor improvement with VERBOSE(99) for printing big objects (using if)
//
// Revision 1.11  2018/05/18 16:40:27  bigatti
// -- added include SparsePolyOps-RingElem.H
//
// Revision 1.10  2018/05/18 12:13:36  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.9  2018/05/17 16:04:04  bigatti
// -- added #include SparsePolyIter
//
// Revision 1.8  2018/01/15 12:41:37  abbott
// Summary: added missing include
//
// Revision 1.7  2017/08/08 14:25:48  bigatti
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
// Revision 1.2  2017/03/13 12:30:37  abbott
// Summary: Move #ifdef guard to after inclusion of header file (so that defns from PREPROCESSOR_DEFNS.H are visible)
//
// Revision 1.1  2017/02/24 08:20:31  bigatti
// -- first import
//

#endif
