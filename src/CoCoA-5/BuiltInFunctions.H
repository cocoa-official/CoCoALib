#ifndef CoCoA5_BuiltInFunctions_H
#define CoCoA5_BuiltInFunctions_H

//   Copyright (c) 2010-2012  John Abbott,  Anna M. Bigatti
//   Original author: 2010  Giovanni Lagorio (lagorio@disi.unige.it)
//
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

#if 0
#include <zip.h>
#endif

#include "CoCoA/library.H"
#include "Interpreter.H"


namespace CoCoA
{
  namespace InterpreterNS
  {

    // theValue was here up to 2021-05
    
    //---------------------------------------
    template <typename T> const T& RefTo(boost::intrusive_ptr<const RightValue>);  

    template <> const BigInt& RefTo<BigInt>(boost::intrusive_ptr<const RightValue>);
    template <> const BigRat& RefTo<BigRat>(boost::intrusive_ptr<const RightValue>);
    template <> const ModuleElem& RefTo<ModuleElem>(boost::intrusive_ptr<const RightValue>);
    template <> const RingElem& RefTo<RingElem>(boost::intrusive_ptr<const RightValue>);
    template <> const RingHom& RefTo<RingHom>(boost::intrusive_ptr<const RightValue>);
    template <> const ideal& RefTo<ideal>(boost::intrusive_ptr<const RightValue>);
    template <> const matrix& RefTo<matrix>(boost::intrusive_ptr<const RightValue>);
    template <> const module& RefTo<module>(boost::intrusive_ptr<const RightValue>);
    template <> const ring& RefTo<ring>(boost::intrusive_ptr<const RightValue>);
    template <> const std::string& RefTo<std::string>(boost::intrusive_ptr<const RightValue>);
    //---- impl ----
    template <> inline const BigInt& RefTo<BigInt>(boost::intrusive_ptr<const RightValue> x){ return intrusive_ptr_cast<const INT>(x)->theBigInt; }
    template <> inline const BigRat& RefTo<BigRat>(boost::intrusive_ptr<const RightValue> x){ return intrusive_ptr_cast<const RAT>(x)->theBigRat; }
    template <> inline const ModuleElem& RefTo<ModuleElem>(boost::intrusive_ptr<const RightValue> x){ return intrusive_ptr_cast<const MODULEELEM>(x)->theModuleElem; }
    template <> inline const RingHom& RefTo<RingHom>(boost::intrusive_ptr<const RightValue> x){ return intrusive_ptr_cast<const RINGHOM>(x)->theRingHom; }
    template <> inline const RingElem& RefTo<RingElem>(boost::intrusive_ptr<const RightValue> x){ return intrusive_ptr_cast<const RINGELEM>(x)->theRingElem; }
    template <> inline const ideal& RefTo<ideal>(boost::intrusive_ptr<const RightValue> x){ return intrusive_ptr_cast<const IDEAL>(x)->theIdeal; }
    template <> inline const matrix& RefTo<matrix>(boost::intrusive_ptr<const RightValue> x){ return intrusive_ptr_cast<const MAT>(x)->theMatrix; }
    template <> inline const module& RefTo<module>(boost::intrusive_ptr<const RightValue> x){ return intrusive_ptr_cast<const MODULE>(x)->theModule; }
    template <> inline const ring& RefTo<ring>(boost::intrusive_ptr<const RightValue> x){ return intrusive_ptr_cast<const RING>(x)->theRing; }
    template <> inline const std::string& RefTo<std::string>(boost::intrusive_ptr<const RightValue> x){ return intrusive_ptr_cast<const STRING>(x)->theString; }
    //-----------------------
    inline std::string ERRORMissingCode(const boost::intrusive_ptr<const RightValue>& x)
    { return "!!!CODE MISSING!!! for "+x->getType()->name; }



    //------------------------------------------------------------
    // The fn  builtIns  behaves like a global variable;  it is filled
    // by the ctor for AddBuiltIn, and the filled vector is used in the
    // the fn initBuiltInFunctions -- see BuiltInFunctions.C

    typedef std::pair<std::string, boost::intrusive_ptr<BuiltInFunction> > NameFunPair;

    std::vector<NameFunPair>& builtIns(); // declared in BuiltInFunctions.C

    struct AddBuiltIn
    {
      inline AddBuiltIn(const std::string &fname, boost::intrusive_ptr<BuiltInFunction> fun)
        {
          builtIns().push_back(std::make_pair(fname, fun));
        }
    };

    //------------------------------------------------------------

#define ARG(N) invocationExpression->args[N]

#define DECLARE_ARITYCHECK_FUNCTION(FNAME) bool builtin_check_arity_##FNAME(int nArg)

#define DECLARE_BUILTIN_FUNCTION(FNAME)                                 \
    intrusive_ptr<RightValue> builtin_##FNAME(const intrusive_ptr<const InvocationExpression> invocationExpression, RuntimeEnvironment * const runtimeEnv); \
    DECLARE_ARITYCHECK_FUNCTION(FNAME) ;                                \
    AddBuiltIn add_##FNAME(#FNAME, new BuiltInFunction(builtin_##FNAME, builtin_check_arity_##FNAME)); \
    intrusive_ptr<RightValue> builtin_##FNAME(const intrusive_ptr<const InvocationExpression> invocationExpression, RuntimeEnvironment * const runtimeEnv)

#define DECLARE_STD_BUILTIN_FUNCTION(FNAME, NARG)                       \
      DECLARE_ARITYCHECK_FUNCTION(FNAME) { return nArg==(NARG); }       \
      DECLARE_BUILTIN_FUNCTION(FNAME) {                                 \
      invocationExpression->checkNumberOfArgs(NARG);

#define END_STD_BUILTIN_FUNCTION }

#define DECLARE_MISSING_EXTLIB(FnName, ExtLibName)                      \
      DECLARE_ARITYCHECK_FUNCTION(FnName) { (void)nArg; return true; }  \
      DECLARE_BUILTIN_FUNCTION(FnName)                                  \
      {                                                                 \
        (void)invocationExpression;                                     \
        (void)runtimeEnv;                                               \
        CoCoA_THROW_ERROR(#FnName " undefined: CoCoA was built without ExternalLib " #ExtLibName, "BuiltInFunction"); \
        return VoidValue::theInstance;                                  \
      }

  } // namespace AST
} // namespace CoCoA

// Close include guard:
#endif
