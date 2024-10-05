#ifndef CoCoA5_BuiltInOneLiners_H
#define CoCoA5_BuiltInOneLiners_H

//   Copyright (c) 2012  John Abbott and Anna M. Bigatti
//   Original author: 2011  Giovanni Lagorio (lagorio@disi.unige.it)
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

#include "BuiltInFunctions.H"

namespace CoCoA
{
  namespace InterpreterNS
  {


#define DECLARE_COCOALIB_FUNCTION0(FUN)                 \
    DECLARE_STD_BUILTIN_FUNCTION(FUN, 0)                \
    (void)runtimeEnv; return Value::from(FUN()); }

#define DECLARE_COCOALIB_FUNCTION1(FUN, ARG1T)                          \
    DECLARE_STD_BUILTIN_FUNCTION(FUN, 1)                                \
    return Value::from(FUN(theValue(runtimeEnv->evalArgAs<ARG1T>(ARG(0))))); }

#define DECLARE_COCOALIB_FUNCTION2(FUN, ARG1T, ARG2T)                   \
    DECLARE_STD_BUILTIN_FUNCTION(FUN, 2)                                \
    return Value::from(FUN(                                             \
                         theValue(runtimeEnv->evalArgAs<ARG1T>(ARG(0))), \
                         theValue(runtimeEnv->evalArgAs<ARG2T>(ARG(1))) \
                            )); }

#define DECLARE_COCOALIB_FUNCTION3(FUN, ARG1T, ARG2T, ARG3T)            \
    DECLARE_STD_BUILTIN_FUNCTION(FUN, 3)                                \
    return Value::from(FUN(                                             \
                         theValue(runtimeEnv->evalArgAs<ARG1T>(ARG(0))), \
                         theValue(runtimeEnv->evalArgAs<ARG2T>(ARG(1))), \
                         theValue(runtimeEnv->evalArgAs<ARG3T>(ARG(2))) \
                            )); }

#define DECLARE_COCOALIB_FUNCTION4(FUN, ARG1T, ARG2T, ARG3T, ARG4T)     \
    DECLARE_STD_BUILTIN_FUNCTION(FUN, 4)                                \
    return Value::from(FUN(                                             \
                         theValue(runtimeEnv->evalArgAs<ARG1T>(ARG(0))), \
                         theValue(runtimeEnv->evalArgAs<ARG2T>(ARG(1))), \
                         theValue(runtimeEnv->evalArgAs<ARG3T>(ARG(2))), \
                         theValue(runtimeEnv->evalArgAs<ARG4T>(ARG(3))) \
                            )); }

#define DECLARE_COCOALIB_FUNCTION5(FUN, ARG1T, ARG2T, ARG3T, ARG4T, ARG5T) \
    DECLARE_STD_BUILTIN_FUNCTION(FUN, 5)                                \
    return Value::from(FUN(                                             \
                         theValue(runtimeEnv->evalArgAs<ARG1T>(ARG(0))), \
                         theValue(runtimeEnv->evalArgAs<ARG2T>(ARG(1))), \
                         theValue(runtimeEnv->evalArgAs<ARG3T>(ARG(2))), \
                         theValue(runtimeEnv->evalArgAs<ARG4T>(ARG(3))), \
                         theValue(runtimeEnv->evalArgAs<ARG5T>(ARG(4))) \
                            )); }

#define DECLARE_COCOALIBFORC5_FUNCTION0(FUN)                    \
    DECLARE_STD_BUILTIN_FUNCTION(FUN, 0)                        \
    (void)runtimeEnv; return Value::from(FUN##_forC5()); }

#define DECLARE_COCOALIBFORC5_FUNCTION1(FUN, ARG1T)                     \
    DECLARE_STD_BUILTIN_FUNCTION(FUN, 1)                                \
    return Value::from(FUN##_forC5(                                     \
                         theValue(runtimeEnv->evalArgAs<ARG1T>(ARG(0))) \
                                    )); }

#define DECLARE_COCOALIBFORC5_FUNCTION2(FUN, ARG1T, ARG2T)              \
    DECLARE_STD_BUILTIN_FUNCTION(FUN, 2)                                \
    return Value::from(FUN##_forC5(                                     \
                         theValue(runtimeEnv->evalArgAs<ARG1T>(ARG(0))), \
                         theValue(runtimeEnv->evalArgAs<ARG2T>(ARG(1))) \
                                    )); }

#define DECLARE_COCOALIBFORC5_FUNCTION3(FUN, ARG1T, ARG2T, ARG3T)       \
    DECLARE_STD_BUILTIN_FUNCTION(FUN, 3)                                \
    return Value::from(FUN##_forC5(                                     \
                         theValue(runtimeEnv->evalArgAs<ARG1T>(ARG(0))), \
                         theValue(runtimeEnv->evalArgAs<ARG2T>(ARG(1))), \
                         theValue(runtimeEnv->evalArgAs<ARG3T>(ARG(2))) \
                                    )); }

#define DECLARE_COCOALIBFORC5_FUNCTION4(FUN, ARG1T, ARG2T, ARG3T, ARG4T) \
    DECLARE_STD_BUILTIN_FUNCTION(FUN, 4)                                \
    return Value::from(FUN##_forC5(                                     \
                         theValue(runtimeEnv->evalArgAs<ARG1T>(ARG(0))), \
                         theValue(runtimeEnv->evalArgAs<ARG2T>(ARG(1))), \
                         theValue(runtimeEnv->evalArgAs<ARG3T>(ARG(2))), \
                         theValue(runtimeEnv->evalArgAs<ARG4T>(ARG(3))) \
                                    )); }

#define DECLARE_COCOALIBFORC5_FUNCTION5(FUN, ARG1T, ARG2T, ARG3T, ARG4T, ARG5T) \
    DECLARE_STD_BUILTIN_FUNCTION(FUN, 5)                                \
    return Value::from(FUN##_forC5(                                     \
                         theValue(runtimeEnv->evalArgAs<ARG1T>(ARG(0))), \
                         theValue(runtimeEnv->evalArgAs<ARG2T>(ARG(1))), \
                         theValue(runtimeEnv->evalArgAs<ARG3T>(ARG(2))), \
                         theValue(runtimeEnv->evalArgAs<ARG4T>(ARG(3))), \
                         theValue(runtimeEnv->evalArgAs<ARG5T>(ARG(4))) \
                                    )); }



  } // namespace AST
} // namespace CoCoA


// Close include guard:
#endif
