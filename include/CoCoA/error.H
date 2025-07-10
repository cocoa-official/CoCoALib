#ifndef CoCoA_error_H
#define CoCoA_error_H

//   Copyright (c)  2005-2020  John Abbott and Anna M. Bigatti
//   Authors:  2005-2010  John Abbott, 2010-2017 John Abbott, Anna M. Bigatti 

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



// All errors should be thrown by calling the function "error"
// as this makes debugging easier and doesn't adversely affect
// run-time performance.  The function "error" will print a brief
// message if CoCoA_DEBUG is set, and it will throw the error object.


#include "CoCoA/exception.H"

#include <string>
//using std::string;
#include <iosfwd>
//using std::ostream;

// UGLY WORKAROUND FOR SOLARIS BUG!!!
// (some header file defines the preprocessor symbol ERR to be 13)
#if defined(ERR)
#  undef ERR
#endif


namespace CoCoA
{

  namespace ERR
  {
    struct ID
    {
      ID(const char* const ErrName, const char* const DefaultMesg): myName(ErrName), myDefaultMesg(DefaultMesg) {};
      bool operator<(const ID& rhs) const noexcept;
      bool operator==(const ID& rhs) const noexcept;
      bool operator!=(const ID& rhs) const noexcept;
      const char* const myName;
      const char* const myDefaultMesg;
    };

    extern const ID
    LANGUAGE,          // Not really an error message, just the name of the language
      nonstandard,       // To create a nonstandard error message
      UNKNOWN,           // To signal an error of unknown type -- can this ever be used???
      ArgTooBig,         // Argument too large in a numerical function
      AssertFail,        // An assertion failed
      BadArg,            // The arg(s) given are unsuitable
      BadArraySize,      // Unsuitable array size (incompatible with other arguments)
      BadCodomain,       // Supplied RingHom has wrong codomain
    //      BadColIndex, // Column index too large in a matrix
      BadCompose,        // Attempt to compose maps with differing domain and codomain
    //      BadComptIndex,     // Index too large accessing component of ModuleElem
      BadConvert,        // Unable to convert value
    //      BadDegIndex, // Index too large accessing degree component
      BadDomain,         // Supplied RingHom has wrong domain
      BadGlobalSettings, // ambiguous, duplicate or incompatible global settings
    //      BadIndetIndex,     // Indeterminate index too large in a polynomial ring
      BadIndex,          // Index out of range
      BadIndetNames,     // Duplicate indet names or varied number of indices for a single name
      BadInducingHom,    // Inducing hom has wrong domain
      BadInducingHom2,   // Inducing hom has wrong codomain
      BadInducingHomKer, // Inducing hom has unsuitable kernel
      BadMatrixSetEntry, // Attempt to set a matrix entry where this is forbidden
      BadMatrixSize,     // Matrix of the wrong size for the requested operation
      BadModulus,        // modulus must be at least 2 and fit into a long
      BadNumBase,        // Base for number conversion must be between 2 and 36 (incl)
      BadNumIndets,      // Wrong number of indet names when creating a PPMonoid or PolyRing
      BadOpenMath,       // OpenMath input did not contain the expected tag
      BadPPMonoid,       // PPMonoid supplied is not of the right type
      BadPPMonoidHomArg, // Argument given to PPMonoidHom is not in the domain
      BadProbability,    // Probability P must satisfy  0 <= P <= 1
      BadPwrZero,        // Attempt to compute 0^N where N is not strictly positive
      BadPolyRingHomImages,// Unsuitable indet images for creating PolyRing hom (wrong number or wrong ring)
      BadQuot,           // Attempt to compute a quotient which does not exist in a given ring/monoid
      BadQuotRing,       // Attempt to quotient by ideal(1)
      BadRing,           // Used ring does not have required properties
      BadRingHomArg,     // Argument given to a RingHom is in the wrong ring
      BadPartialRingHomArg,// Partial RingHom is undefined for this argument
    //      BadRowIndex, // Row index too large in a matrix
      BadSmallFpChar,    // Characteristic (for small finite field) too large or not prime
      BadSymbolHead,     // Invalid symbolic name head (with illegal characters)
      BadSymbolSubscript,// Symbolic name: subscript too large or name does not have that many subscripts
      BadSymbolRange,    // Incompatible names given to range
    //      BLASFailed,        // A BLAS function failed
      CannotReconstruct,
      CanonicalHomFail,  // Unable to construct canonical homomorphism
      ConstMatEntry,     // Cannot assign to constant matrix entry
      DeadMemPool,       // Attempt to use a MemPool after it has been destroyed
      DivByZero,         // Division by zero or a zero-divisor
      ReqNonEmpty,       // List or vector must be non-empty
      ReqCoeffsInField,  // Coefficient ring must be a field
      EmbedBigRatFailed, // Cannot embed rational number into ring
      ExpTooBig,         // Exponent is too large
      ExternalLib,       // Error in external library
      GlobalManager1,    // No GlobalManager object has been created
      GlobalManager2,    // May not create more than one GlobalManager object
      IdealNotInRing,    // Incompatible ring and ideal
      IncompatArgs,      // Args given are incompatible
      IncompatDims,      // Incompatible dimensions
      InputFail,         // An input operation failed
      InsuffPrec,        // Insufficient precision when using RingTwinFloat values
    ////???            IntDivByNeg,       // Integer division/remainder by negative divisor
      InvertibleRingElem,// Non-invertible RingElem required
      IterEnded,         // Attempt to advance an iter which has already reached the end
      LapackFailed,      // A Lapack driver failed
    //      LogZero,           // Cannot compute log of zero
      MemPoolZero,       // Cannot use MemPool to manage blocks of zero size
      MissNumLibs,       // Numerical libraries not configured in
      MixedCoeffRings,   // Arithmetic operation between polys with different coeff rings
      MixedDegrees,      // Arithmetic operation between incompatible degrees
      MixedModules,      // Arithmetic operation between elements of different modules
      MixedPolyIters,    // Comparison between iterators over different polys
      MixedPPMs,         // Arithmetic operation between elements of different PPMonoids
      MixedRings,        // Arithmetic operation between elements of different rings
    //      MixedSizes,        // Operation between objects with different size ==> IncompatDims
    //      ModulusLT2,        // Modulus must be >= 2
      NegExp,            // Negative exponent
      NoGlobalMgr,       // GlobalManager must be created before using CoCoALib
      NotCommutative,    // Ring must be commutative
      NotDenseUPolyRing, // Ring must be a dense univariate polynomial ring, but is not
      NotElemFrF,        // Operation makes sense only for elements of a fraction field
      NotElemGCDDomain,  // Operation makes sense only for elements of a GCD domain
      ReqElemPolyRing,   // Operation makes sense only for elements of a polynomial ring
      NotElemQuotientRing,// Operation makes sense only for elements of a quotient ring
      ReqElemSparsePolyRing,// Operation makes sense only for elements of a sparse polynomial ring
      ReqChar0,          // Characteristic must be 0
      ReqFGModule,       // Module must be Finitely Generated, but is not
      ReqField,          // Ring must be a field, but is not
      NotFracField,      // Ring must be a FractionField, but is not
      ReqFreeModule,     // Module must be free, but is not
      ReqFullRank,       // Matrix must be full rank, but is not
      ReqModuleSpPR,     // Module must be on a SparsePolyRing, but is not
      ReqHomog,          // Input must be homogeneous
      ReqIndet,          // Expected an indeterminate
      NotIntegralDomain, // Ring must be an integral domain, but is not
      NotInvMatrix,      // Matrix is not invertible
      ReqMonomialGens,   // Ideal must be generated by monomials
      ReqNonNegative,    // Value must be non-negative
      ReqNonNegativeGrading,// Grading is not non-negative
      ReqOrdDom,         // Ring must be ordered, but is not
      ReqPolyRing,       // Ring must be a polynomial ring, but is not
      ReqPositive,       // Value is not positive
      ReqPositiveGrading,// Grading is not positive
      NotQuotientRing,   // Ring must be a quotient ring, but is not
      NotRingTwinFloat,  // Operation valid only over RingTwinFloat
      ReqSparsePolyRing, // Ring must be a sparse polynomial ring, but is not
      ReqSquareMatrix,   // Matrix must be square, but is not
      ReqTermOrdering,   // Ordering is not a term-ordering
      NotTrueGCDDomain,  // Ring must be a true GCD domain (and not a field), but is not
      NotUnit,           // Operation is valid only for invertible elements
      ReqUnivariate,     // Polynomial is not univariate
      ReqNonZero,        // Value must be non-zero
      ReqNonZeroGradingDim, // GradingDim (grading dimension) is 0
      ReqGradingDim1,    // GradingDim (grading dimension) must be 1
      ReqNonZeroModulus, // Numerical operation specified a zero modulus
      ReqNonZeroRingElem,// Non-zero RingElem required
      ReqZeroDim,        // Ideal must be 0-dimensional
      NullPtr,           // Null pointer passed where forbidden
      NYI,               // NOT YET IMPLEMENTED
      OBSOLESCENT,       // obsolescent fn called (to avoid this error give option AllowObsolescentFns to GlobalManager)
      OutOfRange,        // Argument is out of range (too big or too small)
      PolyIterEnded,     // Attempt to use an off-the-end iterator
      PPOrder,           // PP is not in the right order
      ShouldNeverGetHere,// Execution should never have reached here
      TimedOut;          // Computation exceeded given time limit
  }

  /*-----------------------------------------------------------------*/
  /** \include error.txt  */
  /*-----------------------------------------------------------------*/
  // Struct to allow error context to be passed to a fn
  // (e.g. if several fns share the same arg checking code)
  struct ErrorContext
  {
  public:
    ErrorContext(const char* fn, const char* file, long line): FnName(fn), FileName(file), LineNum(line) {}
  public: // data members
    const char* FnName;
    const char* FileName;
    long LineNum;
  };

  std::ostream& operator<<(std::ostream& out, const ErrorContext& ErrCtx);

// Macro to create temporary object to convey error context
#define CoCoA_ERROR_CONTEXT  ::CoCoA::ErrorContext(__func__, __FILE__, __LINE__)



  class ErrorInfo: public exception
  {
  public:
    ErrorInfo(const std::string& mesg, const std::string& func);
    ErrorInfo(const std::string& mesg, const ErrorContext& ErrCtx);
    ErrorInfo(const std::string& mesg, const std::string& func, const char* file, unsigned long line);
    ErrorInfo(ERR::ID id, const std::string& func);
    ErrorInfo(ERR::ID id, const ErrorContext& ErrCtx);
    ErrorInfo(ERR::ID id, const std::string& MoreDetails, const ErrorContext& ErrCtx);
    ErrorInfo(ERR::ID id, const std::string& func, const char* file, unsigned long line);
    ~ErrorInfo() {}
///???    const ERR::ID& id() const { return myID; }
    bool operator==(const ERR::ID& id) const noexcept { return myID == id; }
    bool operator!=(const ERR::ID& id) const noexcept { return myID != id; }
    void myOutputSelf(std::ostream& out) const override /*NOT final*/;
  private: // Data members (incl. myMessage & myContext from CoCoA::exception)
    friend void ANNOUNCE(std::ostream& out, const ErrorInfo& err);
    const ERR::ID myID;
    const std::string myFile;
    const unsigned long myLine;
  private: // impl detail
    static std::string ourJoinTogether(const std::string& MainDescr, const std::string& MoreDetails);
  };

  std::ostream& operator<<(std::ostream& out, const ErrorInfo& err);
  void ANNOUNCE(std::ostream& out, const ErrorInfo& err);


// Macro to avoid having to type explicitly __FILE__ and __LINE__
#define CoCoA_STATIC_ERROR_MESG(ErrName, ErrID, FnName)  static const ::CoCoA::ErrorInfo ErrName(ErrID, FnName, __FILE__, __LINE__)

// Helpful Macro: create error obj and throw it (puts in file name and line num)
// !!! NOT YET -- NEEDS MORE THOUGHT !!!   #define CoCoA_THROW_ERROR(ErrID, FnName)  ::CoCoA::ThrowException(::CoCoA::ErrorInfo(ErrID, CoCoA_ERROR_CONTEXT))
#define CoCoA_THROW_ERROR(ErrID, FnName)  ::CoCoA::ThrowException(::CoCoA::ErrorInfo(ErrID, FnName, __FILE__, __LINE__))
#define CoCoA_THROW_ERROR1(ErrID)  ::CoCoA::ThrowException(::CoCoA::ErrorInfo(ErrID, CoCoA_ERROR_CONTEXT))
#define CoCoA_THROW_ERROR2(ErrID, MoreDetails)  ::CoCoA::ThrowException(::CoCoA::ErrorInfo(ErrID, MoreDetails, CoCoA_ERROR_CONTEXT))
#define CoCoA_THROW_ERROR_WITH_CONTEXT2(ErrID, context)  ::CoCoA::ThrowException(::CoCoA::ErrorInfo(ErrID, context))
#define CoCoA_THROW_ERROR_WITH_CONTEXT3(ErrID, MoreDetails, context)  ::CoCoA::ThrowException(::CoCoA::ErrorInfo(ErrID, MoreDetails, context))

// Macro below is obsolete: remove it by Dec 2025
#define CoCoA_ERROR  OBSOLETE___USE_CoCoA_THROW_ERROR

  namespace ErrorLanguage
  {
    void english(); // Set english as the language for error messages.
    void italian(); // Set italian as the language for error messages.
  }



} // end of namespace CoCoA

#endif
