//   Copyright (c)  2005-2017  John Abbott and Anna M. Bigatti
//   Authors:  2005-2017  John Abbott, Anna M. Bigatti

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


#include "CoCoA/error.H"
#include "CoCoA/PREPROCESSOR_DEFNS.H"
#include "CoCoA/verbose.H"

#include <iostream>
using std::ostream;
using std::cerr;
using std::endl;
#include <map>
using std::map;
#include <utility>
using std::make_pair;

namespace CoCoA
{

  // This is simply a forward declaration; the real code is later on
  namespace ErrorLanguage
  {
    static const char* id2message(const ERR::ID& id);
  } // end of namespace ErrorLanguage


  /*static*/ std::string ErrorInfo::ourJoinTogether(const std::string& MainDescr, const std::string& MoreDetails)
  {
    if (MoreDetails.empty())  return MainDescr;
    return MainDescr + ": " + MoreDetails;
  }


  ErrorInfo::ErrorInfo(const std::string& mesg, const std::string& func):
      exception(mesg, func),
      myID(ERR::nonstandard),
      myFile(""),
      myLine(0)
  {}


  ErrorInfo::ErrorInfo(const std::string& mesg, const ErrorContext& ErrCtx):
      exception(mesg, ErrCtx.FnName),
      myID(ERR::nonstandard),
      myFile(ErrCtx.FileName),
      myLine(ErrCtx.LineNum)
  {}


  ErrorInfo::ErrorInfo(const std::string& mesg, const std::string& func, const char* file, unsigned long line):
      exception(mesg, func),
      myID(ERR::nonstandard),
      myFile(file),
      myLine(line)
  {}


  ErrorInfo::ErrorInfo(ERR::ID id, const std::string& func):
      exception(ErrorLanguage::id2message(id), func),
      myID(id),
      myFile(""),
      myLine(0)
  {}


  ErrorInfo::ErrorInfo(ERR::ID id, const ErrorContext& ErrCtx):
      exception(ErrorLanguage::id2message(id), ErrCtx.FnName),
      myID(id),
      myFile(ErrCtx.FileName),
      myLine(ErrCtx.LineNum)
  {}


  ErrorInfo::ErrorInfo(ERR::ID id, const std::string& MoreDetails, const ErrorContext& ErrCtx):
      exception(ourJoinTogether(ErrorLanguage::id2message(id), MoreDetails), ErrCtx.FnName),
      myID(id),
      myFile(ErrCtx.FileName),
      myLine(ErrCtx.LineNum)
  {}


  ErrorInfo::ErrorInfo(ERR::ID id, const std::string& func, const char* file, unsigned long line):
      exception(ErrorLanguage::id2message(id), func),
      myID(id),
      myFile(file),
      myLine(line)
  {}


  void ErrorInfo::myOutputSelf(std::ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams
    out << "CoCoA::ErrorInfo(" << myID.myName << ", \"" << message(*this) << "\")";
  }

  std::ostream& operator<<(std::ostream& out, const ErrorInfo& err)
  {
    err.myOutputSelf(out);
    return out;
  }


  void ANNOUNCE(std::ostream& out, const ErrorInfo& err)
  {
    if (!out) return;  // short-cut for bad ostreams (but in this case...???)
    // Save stream flags to restore them upon exit.
    const std::ios::fmtflags OrigFlags = out.flags();
    out.setf(std::ios::dec, std::ios::basefield); // Want line no. in decimal!
    out << endl
        << endl
        << "***************************************************************************" << endl
        << "****CoCoA ERROR****  ErrCode: " << err.myID.myName << endl
        << "****CoCoA ERROR****  Message: " << message(err) << endl;
    if (!context(err).empty())
      out << "****CoCoA ERROR****  Context: " << context(err) << endl;
    if (!err.myFile.empty())
    {
      out << "****CoCoA ERROR****  File:    " << err.myFile << endl;
      out << "****CoCoA ERROR****  Line:    " << err.myLine << endl;
    }
    out << "***************************************************************************" << endl;
    out << endl;
    out.setf(OrigFlags);
  }


  ///////////////////////////////////////////////////////////////////////////
  // The error IDs and their default messages...

  namespace ERR
  {
// Nasty macro because someone wants the error code available as a string!! (sigh)
#define DEFINE_ERROR(ErrCode, ErrMesg) const ID ErrCode("CoCoA::ERR::" #ErrCode, ErrMesg)

    DEFINE_ERROR(LANGUAGE, "English");
    DEFINE_ERROR(nonstandard, "Error code corresponding to an error different from the standard CoCoA errors");
    DEFINE_ERROR(UNKNOWN, "UNKNOWN ERROR CODE PRODUCED -- please inform the CoCoA Team of this");
    DEFINE_ERROR(ArgTooBig, "Argument to a numerical function too large (value would be too big)");
    DEFINE_ERROR(AssertFail, "Assertion failed");
    DEFINE_ERROR(BadArg, "The arg(s) given are unsuitable");
    DEFINE_ERROR(BadArraySize, "Array size is incompatible with other arguments");
    DEFINE_ERROR(BadCodomain, "RingHom has wrong codomain");
    //    DEFINE_ERROR(BadColIndex, "Column index out of range");
    DEFINE_ERROR(BadCompose, "Attempt to compose maps with differing domain and codomain");
    //    DEFINE_ERROR(BadComptIndex, "Index too large accessing component of ModuleElem");
    DEFINE_ERROR(BadConvert, "Unable to convert value");
    //    DEFINE_ERROR(BadDegIndex, "Index too large accessing degree component");
    DEFINE_ERROR(BadDomain, "RingHom has wrong domain");
    DEFINE_ERROR(BadGlobalSettings, "Ambiguous, duplicate or incompatible global settings");
    //    DEFINE_ERROR(BadIndetIndex, "Indeterminate index out of range");
    DEFINE_ERROR(BadIndetNames, "Duplicate indet names or varied number of indices for a single name");
    DEFINE_ERROR(BadIndex, "Index out of range");
    DEFINE_ERROR(BadInducingHom, "Inducing hom has wrong domain");
    DEFINE_ERROR(BadInducingHom2, "Inducing hom has wrong codomain");
    DEFINE_ERROR(BadInducingHomKer, "Inducing hom has unsuitable kernel");
    DEFINE_ERROR(BadMatrixSetEntry, "Attempt to set a matrix entry where this is forbidden");
    DEFINE_ERROR(BadMatrixSize, "Matrix size is wrong for this operation");
    DEFINE_ERROR(BadModulus, "Modulus must be at least 2 and fit into a long");
    DEFINE_ERROR(BadNumBase, "Base for number conversion must be between 2 and 36 (incl)");
    DEFINE_ERROR(BadNumIndets, "Wrong number of indet names when creating a PPMonoid");
    DEFINE_ERROR(BadOpenMath, "OpenMath input did not contain the expected tag");
    DEFINE_ERROR(BadPPMonoid, "PPMonoid is not of the right type");
    DEFINE_ERROR(BadPPMonoidHomArg, "Argument to PPMonoidHom must be in domain");
    DEFINE_ERROR(BadPolyRingHomImages, "Indet images must be in codomain");  // rename ReqElemCodomain "must be in codomain"?
    DEFINE_ERROR(BadProbability, "Probability P must satisfy  0 <= P <= 1");
    DEFINE_ERROR(BadPwrZero, "Power of zero must be positive");
    DEFINE_ERROR(BadQuot, "Inexact division (i.e. quotient does not exist in ring or monoid)");
    DEFINE_ERROR(BadQuotRing, "Attempt to quotient by ideal(1)");
    DEFINE_ERROR(BadRing, "Unsuitable ring");
    DEFINE_ERROR(BadRingHomArg, "Argument to RingHom must be in domain");
    DEFINE_ERROR(BadPartialRingHomArg, "Partial RingHom is undefined for this argument");
    //    DEFINE_ERROR(BadRowIndex, "Row index out of range");
    DEFINE_ERROR(BadSmallFpChar, "Characteristic (for small finite field) too large or not prime");
    DEFINE_ERROR(BadSymbolHead, "Illegal character in symbol head");
    DEFINE_ERROR(BadSymbolSubscript, "Symbol name: subscript too large or name does not have that many indices");
    DEFINE_ERROR(BadSymbolRange, "Incompatible symbols given to range");
    //    DEFINE_ERROR(BLASFailed, "A BLAS function failed");
    DEFINE_ERROR(CannotReconstruct, "Unable to reliably reconstruct");
    DEFINE_ERROR(CanonicalHomFail, "Unable to construct canonical homomorphism");
    DEFINE_ERROR(ConstMatEntry, "Cannot assign to constant matrix entry");
    DEFINE_ERROR(DeadMemPool, "Attempt to use a MemPool after it has been destroyed");
    DEFINE_ERROR(DivByZero, "Division by zero or by a zero-divisor");
    DEFINE_ERROR(ReqCoeffsInField, "Coefficient ring must be a field");
    DEFINE_ERROR(EmbedBigRatFailed, "Cannot embed rational number into ring");
    DEFINE_ERROR(ReqNonEmpty, "List, or vector, must be non-empty");
    DEFINE_ERROR(ExpTooBig, "Exponent is too large");
    DEFINE_ERROR(ExternalLib, "Error in external library");
    DEFINE_ERROR(GlobalManager1, "You must create a GlobalManager object before using CoCoALib features");
    DEFINE_ERROR(GlobalManager2, "Attempt to create more than one GlobalManager object");
    DEFINE_ERROR(IdealNotInRing, "Incompatible ring and ideal: ideal is not in the given ring");
    DEFINE_ERROR(IncompatArgs, "Incompatible arguments");
    DEFINE_ERROR(IncompatDims, "Incompatible dimensions");
    DEFINE_ERROR(InputFail, "An input operation failed (e.g. premature EOF or wrong type of input)");
    DEFINE_ERROR(InsuffPrec, "RingTwinFloat: insufficient precision (e.g. error growth impedes further computation)");
///???    DEFINE_ERROR(IntDivByNeg, "Integer division by negative divisor with non-zero remainder");
    DEFINE_ERROR(InvertibleRingElem, "Non-invertible RingElem required");
    DEFINE_ERROR(IterEnded, "Attempt to advance an iter which has already reached the end");
    DEFINE_ERROR(LapackFailed, "A Lapack driver failed");
    //    DEFINE_ERROR(LogZero, "Cannot compute log of zero");
    DEFINE_ERROR(MemPoolZero, "Cannot use MemPool to manage blocks of zero size");
    DEFINE_ERROR(MissNumLibs, "Numerical libraries not configured in");
    DEFINE_ERROR(MixedCoeffRings, "Arithmetic operation between polys with different coeff rings");
    DEFINE_ERROR(MixedDegrees, "Arithmetic operation between incompatible degrees");
    DEFINE_ERROR(MixedModules, "Elements must be in the same module; mixed module operations are forbidden");
    DEFINE_ERROR(MixedPolyIters, "Comparison between iterators over different polys");
    DEFINE_ERROR(MixedPPMs, "Elements must be in the same PPMonoid; mixed PPMonoid operations are forbidden");
    DEFINE_ERROR(MixedRings, "Elements must be in the same ring; mixed ring operations are forbidden");
    //    DEFINE_ERROR(MixedSizes, "Operation between objects with different size"); ==> IncompatDims
    //    DEFINE_ERROR(ModulusLT2, "Modulus must be >= 2");
    DEFINE_ERROR(NegExp, "Exponent must be non-negative");
    DEFINE_ERROR(NoGlobalMgr, "GlobalManager must be created before using CoCoALib");
    DEFINE_ERROR(NotCommutative, "Ring is not commutative");
    DEFINE_ERROR(NotDenseUPolyRing, "Ring must be a dense univariate polynomial ring");
    DEFINE_ERROR(NotElemFrF, "RingElem is not in a fraction field");
    DEFINE_ERROR(NotElemGCDDomain, "RingElem is not in a GCD domain");
    DEFINE_ERROR(ReqElemPolyRing, "RingElem is not in a polynomial ring");
    DEFINE_ERROR(NotElemQuotientRing, "RingElem is not in quotient ring");
    DEFINE_ERROR(ReqElemSparsePolyRing, "RingElem is not in a sparse polynomial ring");
    DEFINE_ERROR(ReqFGModule, "Module must be FGModule, but is not");
    DEFINE_ERROR(ReqField, "Ring must be a field");
    DEFINE_ERROR(NotFracField, "Ring must be a FractionField");
    DEFINE_ERROR(ReqFreeModule, "Module must be a free module");
    DEFINE_ERROR(ReqFullRank, "Matrix must be full rank");
    DEFINE_ERROR(ReqChar0, "Characteristic must be 0");
    DEFINE_ERROR(ReqModuleSpPR, "Module must be a over a SparsePolyRing");
    DEFINE_ERROR(ReqHomog, "Input must be homogeneous");
    DEFINE_ERROR(ReqIndet, "Expected an indeterminate");
    DEFINE_ERROR(NotIntegralDomain, "Ring is not an integral domain");
    DEFINE_ERROR(NotInvMatrix, "Matrix must be invertible");
    DEFINE_ERROR(ReqMonomialGens, "Ideal generators must be monomial");
    DEFINE_ERROR(ReqNonNegative, "Value must be non-negative");
    DEFINE_ERROR(ReqNonNegativeGrading, "Grading must be non-negative");
    DEFINE_ERROR(ReqNonZero, "Value must be non-zero");
    DEFINE_ERROR(ReqNonZeroGradingDim, "GradingDim (grading dimension) must be non 0");
    DEFINE_ERROR(ReqGradingDim1, "GradingDim (grading dimension) must be 1");
    DEFINE_ERROR(ReqNonZeroModulus, "A zero modulus was specified for a numerical operation");
    DEFINE_ERROR(ReqNonZeroRingElem, "RingElem must be non-zero");
    DEFINE_ERROR(ReqOrdDom, "Ring is not an ordered domain");
    DEFINE_ERROR(ReqPolyRing, "Ring must be a polynomial ring");
    DEFINE_ERROR(ReqPositive, "Value must be positive");
    DEFINE_ERROR(ReqPositiveGrading, "Grading must be positive");
    DEFINE_ERROR(NotQuotientRing, "Ring must be a quotient ring");
    DEFINE_ERROR(NotRingTwinFloat, "Operation valid only over RingTwinFloat");
    DEFINE_ERROR(ReqSparsePolyRing, "Ring must be a sparse polynomial ring");
    DEFINE_ERROR(ReqSquareMatrix, "Matrix must be square");
    DEFINE_ERROR(NotTrueGCDDomain, "Ring is not a true GCD domain (e.g. it is a field)");
    DEFINE_ERROR(ReqTermOrdering, "Ordering must be a term-ordering (i.e. all indets>1)");
    DEFINE_ERROR(NotUnit, "Cannot invert non-unit");
    DEFINE_ERROR(ReqUnivariate, "Polynomial must be univariate");
    DEFINE_ERROR(ReqZeroDim, "Ideal must be 0-dimensional");
    DEFINE_ERROR(NullPtr, "Null pointer passed where forbidden");
    DEFINE_ERROR(NYI, "NOT YET IMPLEMENTED -- please be patient, we're working on it");
    DEFINE_ERROR(OBSOLESCENT, "obsolescent fn called (to avoid this error give option AllowObsolescentFns to GlobalManager)");
    DEFINE_ERROR(OutOfRange, "Argument is out of range (too big or too small)");
    DEFINE_ERROR(PPOrder, "PP is not in the right order");
    DEFINE_ERROR(PolyIterEnded, "Attempt to use an off-the-end iterator");
    DEFINE_ERROR(ShouldNeverGetHere, "Execution should never get here -- please inform the CoCoA Team");
    DEFINE_ERROR(TimedOut, "Computation exceeded given time limit");
  } // end of namespace ERR


  // These three could be inline, but I don't believe speed is important here.
  // I have defined them outside namespace ERR for improved readability
  // (well, I find ERR::ID usefully explicit).
  bool ERR::ID::operator<(const ERR::ID& rhs) const noexcept
  {
    return myName < rhs.myName; // simply compare pointer values
  }

  bool ERR::ID::operator==(const ERR::ID& rhs) const noexcept
  {
    return myName == rhs.myName; // simply compare pointer values
  }

  bool ERR::ID::operator!=(const ERR::ID& rhs) const noexcept
  {
    return myName != rhs.myName; // simply compare pointer values
  }


//-------------------------------------------------------

  namespace ErrorLanguage
  {

    typedef map<ERR::ID, const char*> MsgTable_t;

    /*const*/ MsgTable_t* GlobalErrorMsgTablePtr = nullptr;

    static const char* id2message(const ERR::ID& id)
    {
      if (GlobalErrorMsgTablePtr == nullptr) return id.myDefaultMesg;  // default message (in english)
      if (GlobalErrorMsgTablePtr->count(id) > 0)
        return (*GlobalErrorMsgTablePtr)[id];                // return translated message

      // There is no translation of the message associated to id, so use the default message.
      return id.myDefaultMesg; // ??? add apology for absence of translated message???
    }


    static void SetLanguage(MsgTable_t& tbl)
    {
      if (tbl.count(ERR::LANGUAGE) == 0) return; // silently refuse to set language if there is no language name
      GlobalErrorMsgTablePtr = &tbl;
    }


    static void set(MsgTable_t& tbl, ERR::ID id, const char* const msg)
    {
      if (tbl.insert(std::make_pair(id, msg)).second) return;
      cerr << "WARNING: Slot for error `" << id.myName << "' already occupied\n"
           << "WARNING: Failed to register message `" << msg << "'" << endl;
      throw ERR::ShouldNeverGetHere;
    }


    ///////////////////////////////////////////////////////////////////////////
    // Selecting english is a special cases as it is the default language.

    void english()
    {
      // Setting the pointer should not leak memory as it points to static
      // data when it is not nullptr.
      GlobalErrorMsgTablePtr = nullptr;
    }


    ///////////////////////////////////////////////////////////////////////////
    // Some Italian error messages are listed here.
    void italian()
    {
      static MsgTable_t MsgTable;
      if (MsgTable.empty())
      {
        set(MsgTable, ERR::LANGUAGE,     "Italian/Italiano");
        set(MsgTable, ERR::nonstandard,  "Codice d'errore che non appartiene all'insieme degli errori standard");
        set(MsgTable, ERR::UNKNOWN,      "CODICE D'ERRORE SCONOSCIUTO -- prego di segnalarlo al CoCoA Team");
        set(MsgTable, ERR::ShouldNeverGetHere,      "Accaduto un errore grave nella libreria CoCoA -- prego di segnalarlo al CoCoA Team");
      }
      SetLanguage(MsgTable);
    }

  } // end of namespace ErrorLanguage



  std::ostream& operator<<(std::ostream& out, const ErrorContext& ErrCtx)
  {
    if (!out) return out;
    out << "ErrorContext(FileName=" << ErrCtx.FileName << ",  LineNum=" << ErrCtx.LineNum << ",  FnName=" << ErrCtx.FnName << ")";
      return out;
  }

} // end of namespace CoCoA
