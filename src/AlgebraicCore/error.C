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
    DEFINE_ERROR(IncompatArgs, "Args given are incompatible");
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


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/error.C,v 1.93 2024/09/04 15:06:36 bigatti Exp $
// $Log: error.C,v $
// Revision 1.93  2024/09/04 15:06:36  bigatti
// Summary: Error codes: removed MixedSizes (into IncompatDims),
// removed BLASFailed, added ExternalLibs and ReqGradingDim1,
// changed ExpectedCoeffsInField into ReqCoeffsInField
//
// Revision 1.92  2024/08/28 15:21:56  bigatti
// Summary: Error renamed: Not ==> ReqSquareMatrix, ReqFullRank, ReqTermOrdering
//
// Revision 1.91  2024/08/28 13:23:34  bigatti
// Summary: added ReqZeroDim
//
// Revision 1.90  2024/08/02 16:06:39  bigatti
// Summary: added comment about elem codomain
//
// Revision 1.89  2024/08/01 16:47:35  bigatti
// Summary: error codes: Req(Elem)SparsePolyRing;   string message for BadRing
//
// Revision 1.88  2024/07/31 14:15:40  abbott
// Summary: Added defn of ReqHomog
//
// Revision 1.87  2024/07/23 16:38:30  bigatti
// Summary: error codes: ERR::ReqMonomialGens, and using new CoCoA_THROW_ERROR1/2
//
// Revision 1.86  2024/07/15 17:12:42  bigatti
// Summary: Error codes:  added CoCoA_THROW_ERROR1(ErrID);
// Removed BadRowIndex, BadColIndex;  Renamed  ERR::Empty --> ERR::ReqNonEmpty
//
// Revision 1.85  2024/07/14 10:27:51  abbott
// Summary: Added new ctor with extra strig arg (MoreDetails); also new macro
//
// Revision 1.84  2024/07/11 15:44:09  bigatti
// Summary: error codes: NotField to ReqField
//
// Revision 1.83  2024/07/03 15:39:46  bigatti
// Summary: removed unused error code ModulusLT2
//
// Revision 1.82  2024/07/03 14:46:49  bigatti
// Summary: error codes: Not... changed into ReqFGModule, ReqFreeModule, ReqModuleSpPR
//
// Revision 1.81  2024/07/03 09:44:11  bigatti
// Summary: error codes (redmine 92):
// removed BadIntedIndex, BadDegIndex, BadComptIndex
//
// Revision 1.80  2024/07/02 20:42:07  abbott
// Summary: Renaming of errors (redmine 308)
//
// Revision 1.79  2024/07/02 15:37:50  bigatti
// Summary: Changed error codes: LogZero into ReqNonZero
// and Not... into ReqNonNegative, ReqNonNegativeGrading, ReqPositive, ReqPositiveGrading
//
// Revision 1.78  2024/07/02 10:00:09  bigatti
// Summary: error codes, first batch:
// ReqUnivariate, ReqNonZero, ReqNonZeroGradingDim, ReqNonZeroModulus,  ReqNonZeroRingElem
//
// Revision 1.77  2024/03/28 09:11:22  bigatti
// Summary: just a comment: suggestion to rename NotNonZero into ReqNonZero
//
// Revision 1.76  2023/02/23 20:52:00  abbott
// Summary: Added new ErrorContext struct & macro
//
// Revision 1.75  2022/02/18 14:12:02  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.74  2021/02/10 19:40:01  abbott
// Summary: Added noexcept (sometimes instead of throw()) -- redmine 1572
//
// Revision 1.73  2020/06/19 19:42:59  abbott
// Summary: Cleaned; removed cruft
//
// Revision 1.72  2020/06/19 19:39:21  abbott
// Summary: Now all throws go through new template fn ThrowException; seems to work
//
// Revision 1.71  2020/06/19 14:51:45  abbott
// Summary: Eliminated macro CoCoA_THROW (caused trouble)
//
// Revision 1.70  2020/02/11 16:56:42  abbott
// Summary: Corrected last update (see redmine 969)
//
// Revision 1.69  2020/02/11 16:12:20  abbott
// Summary: Added some checks for bad ostream (even to mem fns myOutput); see redmine 969
//
// Revision 1.68  2019/10/15 11:54:41  abbott
// Summary: Changed 0 into nullptr (where appropriate)
//
// Revision 1.67  2019/03/19 11:07:08  abbott
// Summary: Replaced 0 by nullptr where appropriate
//
// Revision 1.66  2018/04/18 14:14:08  abbott
// Summary: Added new errors: IncompatArgs and NullPtr
//
// Revision 1.65  2018/03/20 11:38:36  bigatti
// -- added error "ExpectedCoeffsInField"
//
// Revision 1.64  2018/03/14 14:30:22  abbott
// Summary: Added new error OutOfRange
//
// Revision 1.63  2017/12/12 14:19:34  abbott
// Summary: Corrected use of CoCoA_DEBUG CPP symbol
//
// Revision 1.62  2017/11/29 20:34:17  abbott
// Summary: Added SturmSeq and NumRealRoots
//
// Revision 1.61  2017/09/06 14:10:05  abbott
// Summary: Added ERR::TimedOut
//
// Revision 1.60  2017/09/06 11:56:30  abbott
// Summary: Changed ERR::SERIOUS into ERR::ShouldNeverGetHere
//
// Revision 1.59  2017/07/22 12:57:42  abbott
// Summary: Updated rtn type of myOutputSelf
//
// Revision 1.58  2017/05/11 08:46:46  bigatti
// -- added error ERR::CannotReconstruct
// -- cleaned up code accordingly
//
// Revision 1.57  2017/05/09 13:51:30  bigatti
// -- changed BadRingHomArg2 --> BadPartialRingHomArg
//
// Revision 1.56  2017/04/26 09:15:42  bigatti
// -- simplified messages for symbols
//
// Revision 1.55  2017/04/21 15:45:15  bigatti
// -- removed word "supplied"
// -- simplified BadArraySize
//
// Revision 1.54  2017/04/21 13:46:48  bigatti
// -- changed some "is not" into "must be"
//
// Revision 1.53  2017/03/13 12:22:12  abbott
// Summary: Added include of PREPROCESSOR_DEFNS.H, because file uses CoCoA_DEBUG
//
// Revision 1.52  2016/11/11 14:15:34  abbott
// Summary: Added short-cut to operator<< when ostream is in bad state
//
// Revision 1.51  2016/11/04 20:40:03  abbott
// Summary: Renamed ERR::OBSOLETE to ERR::OBSOLESCENT
//
// Revision 1.50  2016/11/03 12:29:58  abbott
// Summary: Added file for obsolescent fns; also there is a global flag saying whether to give error if calling one.
//
// Revision 1.49  2016/09/08 07:53:18  bigatti
// -- changed mesg in BadIndetIndex
//
// Revision 1.48  2015/06/26 14:57:00  abbott
// Summary: Now ErrorInfo derives from CoCoA::exception
// Author: JAA
//
// Revision 1.47  2015/03/04 11:25:42  bigatti
// -- changed message for BadPolyRingHomImages
//
// Revision 1.46  2014/07/31 14:45:19  abbott
// Summary: Merged io.H and UtilsTemplate.H into new header VectorOperations.H
// Author: JAA
//
// Revision 1.45  2014/07/10 14:33:43  abbott
// Summary: Removed unused error code NotRingZ (or NotRingZZ, there was a typo!)
// Author: JAA
//
// Revision 1.44  2014/05/14 07:36:26  bigatti
// -- fixed typo
//
// Revision 1.43  2013/05/27 17:05:38  bigatti
// -- changed name to error for free module over SparsePolyRing
//
// Revision 1.42  2013/05/20 15:58:17  abbott
// Added new error code NotPositive.
//
// Revision 1.41  2013/02/21 16:56:49  bigatti
// -- added ERR:Empty
//
// Revision 1.40  2013/02/13 09:08:38  bigatti
// -- added ZeroGradingDim
//
// Revision 1.39  2013/01/21 13:37:29  abbott
// Added new error IncompatDims.
//
// Revision 1.38  2012/12/17 14:32:48  abbott
// Improved message for GlobalManager1.
//
// Revision 1.37  2012/05/24 14:49:23  bigatti
// -- changed symbol "index" into "subscripts"
//
// Revision 1.36  2012/05/22 10:02:37  abbott
// Removed IsGCDDomain; substituted by IsTrueGCDDomain.
// Added IsFractionFieldOfGCDDomain.
//
// Revision 1.35  2012/03/26 11:51:35  abbott
// Rephrased error message associated with ERR::NotIndet.
//
// Revision 1.34  2012/03/16 15:51:35  abbott
// Added new error NotNonZero.
//
// Revision 1.33  2012/02/08 16:16:55  bigatti
// -- changed: Z,Q -> ZZ,QQ
//
// Revision 1.32  2011/08/24 10:32:04  bigatti
// -- renamed QQ --> BigRat
// -- sorted #include
//
// Revision 1.31  2011/06/23 16:02:26  abbott
// Minor changes to two error message strings.
//
// Revision 1.30  2011/03/16 13:22:43  abbott
// Added new error "BadModulus".
//
// Revision 1.29  2011/01/31 14:11:11  bigatti
// -- changed BadRowIndex, BadColIndex
//
// Revision 1.28  2011/01/19 16:12:18  bigatti
// -- added ERR::ModulusLT2
//
// Revision 1.27  2011/01/17 17:41:01  bigatti
// -- improved error message for BadQuot
//
// Revision 1.26  2010/12/17 16:09:15  abbott
// Changed interface to ANNOUNCE function: now it requires the ostream to passed in.
// A few other minor changes.
//
// Revision 1.25  2010/10/22 14:03:03  abbott
// Major change to GMPAllocator -- it is now set/activated by the GlobalManager.
// This is a Friday afternoon check-in... hope to check in cleaner code in the
// next few days.
//
// Revision 1.24  2010/07/09 17:05:01  abbott
// Added new error for PPMonoid homs.
//
// Revision 1.23  2010/04/27 16:07:50  bigatti
// -- added 2 error messages for DynamicBitset
//
// Revision 1.22  2010/03/18 13:54:42  abbott
// Added new error for OpenMath input problems.
//
// Revision 1.21  2010/02/03 15:24:39  bigatti
// -- added NotMonomialGens
//
// Revision 1.20  2009/12/23 18:53:51  abbott
// Major change to conversion functions:
//   convert(..) is now a procedure instead of a function
//   IsConvertible(..) replaces the former convert(..) function
//   Added new NumericCast conversion function (placeholder for BOOST feature)
//   Consequent changes in code which uses these features.
//
// Revision 1.19  2009/12/11 11:46:32  abbott
// Changed fn  convert  into  IsConvertible.
// Added template procedure  convert.
// New version because change is not backward compatible.
//
// Revision 1.18  2009/12/03 17:26:34  abbott
// Renamed EFGModule to FGModule.
// Renamed ModuleBase member fns  myInit -> myNew, myKill -> myDelete.
// Removed some cruft (old code that was not used by anyone).
//
// Revision 1.17  2009/10/29 18:47:58  abbott
// Added IterEnded error code.
//
// Revision 1.16  2009/07/02 16:32:10  abbott
// Consequential changes stemming from new class QQ, and modified interface to the member
// function RingBase::myIsRational.  Also some new conversion functions.
//
// Revision 1.15  2009/05/20 14:25:30  abbott
// Added new error for InputFailCheck.
//
// Revision 1.14  2008/12/16 21:11:42  abbott
// Improved error message for bad symbol index -- there are two possible causes.
//
// Revision 1.13  2008/11/20 10:50:13  abbott
// Minor change to printed message in ANOUNCE.
//
// Revision 1.12  2008/10/08 13:54:05  abbott
// New minor version to reflect the backward-incompatible changes to CoCoA errors.
//
// Revision 1.11  2008/10/08 09:48:19  abbott
// Final(?) minor refinement to CoCoA::ERR:ID internal layout.
//
// Revision 1.10  2008/10/08 08:20:12  abbott
// Minor refinement to new error implementation & corr. change to doc.
//
// Revision 1.9  2008/10/07 15:45:22  abbott
// Changed ErrorInfo objects so they include the name of their own error ID.
// Changed catch statements to catch const objects.
// Removed calls to the member fn which accessed the error ID member of an
// ErrorInfo; now you simply compare directly with the error ID (makes the
// code easier to read).
//
// Revision 1.8  2008/04/21 12:51:53  abbott
// Added some new error codes.
//
// Revision 1.7  2008/03/12 16:39:31  bigatti
// -- changed: ERR:ZeroPoly --> ERR::ZeroRingElem
// -- added: ERR:InvertibleRingElem
//
// Revision 1.6  2007/10/30 17:14:06  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.5  2007/10/05 14:35:16  bigatti
// -- added error: NotDenseUPolyRing
//
// Revision 1.4  2007/06/21 21:29:47  abbott
// Changed name of RingFloat into RingTwinFloat.
//
// Revision 1.3  2007/05/21 12:57:28  abbott
// New class for passing machine integers as args; includes some simple
// operations on machine integers (cmp, gcd, IsNegative,...).  Operations
// between ZZ and machine integers modified to use the new class.  Inexact
// integer division (of a ZZ) by a negative value now triggers an error;
// new error for reporting inexact integer division by a negative value.
//
// Revision 1.2  2007/03/28 10:06:13  abbott
// Now gives error when you use RingZ() or RingQ() without creating GlobalManager.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.12  2007/03/08 18:22:28  cocoa
// Just whitespace cleaning.
//
// Revision 1.11  2007/03/08 11:07:12  cocoa
// Made pseudo ctors for polynomial rings more uniform.  This allowed me to
// remove an include of CoCoA/symbol.H  from the RingDistrM*.H files, but then
// I had to put the include in several .C files.
//
// Revision 1.10  2007/03/05 21:25:57  cocoa
// Forgot to check these in a few minutes ago.
//
// Revision 1.9  2007/03/05 16:17:11  bigatti
// -- clenup for numerical code (and 3 new error codes)
//
// Revision 1.8  2007/03/03 14:07:23  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.7  2007/03/02 10:47:53  cocoa
// First stage of RingZ modifications -- tests do not compile currently, Anna will fix this.
//
// Revision 1.6  2007/02/12 17:40:19  bigatti
// -- added MissNumLibs
//
// Revision 1.5  2007/01/13 14:14:34  cocoa
// Overhaul of RingHom code: it nows uses SmartPtrIRC, and printing is more logical.
// Have not yet updated the documentation.
//
// Revision 1.4  2006/12/21 13:48:32  cocoa
// Made all increment/decrement calls prefix (except where the must be postfix).
//
// Revision 1.3  2006/08/07 21:23:25  cocoa
// Removed almost all publicly visible references to SmallExponent_t;
// changed to long in all PPMonoid functions and SparsePolyRing functions.
// DivMask remains to sorted out.
//
// Revision 1.2  2006/07/17 19:16:53  cocoa
// New, better errors :-)
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.9  2006/04/14 13:51:40  cocoa
// Added a class for generating random bits (RandomBitStream in the file random.H).
//
// Revision 1.8  2006/04/07 16:44:52  cocoa
// Considerably updated MemPool design -- it works, and I'm about to test
// its efficiency against the old one.
//
// Revision 1.7  2006/04/04 13:01:02  cocoa
// -- added: NotNonNegativeGrading, NotPositiveGrading
//
// Revision 1.6  2006/03/27 12:21:25  cocoa
// Minor silly changes to reduce number of complaints from some compiler or other.
//
// Revision 1.5  2006/03/12 21:28:33  cocoa
// Major check in after many changes
//
// Revision 1.4  2005/12/31 12:22:17  cocoa
// Several minor tweaks to silence the Microsoft compiler:
//  - added some missing #includes and using directives
//  - moved some function defns into the right namespace
//  - etc.
//
// Revision 1.3  2005/11/24 16:09:38  cocoa
// -- added operator[] for ModuleElem
//
// Revision 1.2  2005/10/26 14:34:51  cocoa
// Changed ANNOUNCE: it now sets sane ostream flags for itself, and resets
// the original flags before exiting.  No more line numbers in hexadecimal!
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.8  2005/09/28 11:50:34  cocoa
// -- new code for graded modules
//
// Revision 1.7  2005/08/08 16:36:32  cocoa
// Just checking in before going on holiday.
// Don't really recall what changes have been made.
// Added IsIndet function for RingElem, PPMonoidElem,
// and a member function of OrdvArith.
// Improved the way failed assertions are handled.
//
// Revision 1.6  2005/07/19 15:30:20  cocoa
// A first attempt at iterators over sparse polynomials.
// Main additions are to SparsePolyRing, DistrMPoly*.
// Some consequential changes to PPMonoid*.
//
// Revision 1.5  2005/07/15 16:34:33  cocoa
// Added iterators for sparse polynomials.
// The code compiles (and the old tests still run).
// It'd Friday evening -- I'm going home before
// getting any ideas about making the iterator code run.
//
// Revision 1.4  2005/07/08 15:09:28  cocoa
// Added new symbol class (to represent names of indets).
// Integrated the new class into concrete polynomial rings
// and PPMonoid -- many consequential changes.
// Change ctors for the "inline" sparse poly rings: they no
// longer expect a PPMonoid, but build their own instead
// (has to be a PPMonoidOv).
//
// Revision 1.3  2005/07/01 16:08:15  cocoa
// Friday check-in.  Major change to structure under PolyRing:
// now SparsePolyRing and DUPolyRing are separated (in preparation
// for implementing iterators).
//
// A number of other relatively minor changes had to be chased through
// (e.g. IndetPower).
//
// Revision 1.2  2005/06/27 14:55:24  cocoa
// Cleaned up some more PPMonoid code.
//
// Revision 1.1.1.1  2005/05/03 15:47:31  cocoa
// Imported files
//
// Revision 1.8  2005/04/29 15:42:02  cocoa
// Improved documentation for GMPAllocator.
// Added example program for GMPAllocator.
// Added example program for simple ops on polynomials.
// Added two new ctors for (principal) ideals (from long, and from ZZ).
// Added (crude) printing for PPMonoids.
// Updated library.H (#included GMPAllocator.H).
//
// Revision 1.7  2005/04/27 16:14:56  cocoa
// Cleaned up example programs -- added "free use" permit.
// Changed a couple of ErrorInfo object names, and added
// ERR::NotTrueGCDDomain.
//
// Revision 1.6  2005/04/20 15:40:47  cocoa
// Major change: modified the standard way errors are to be signalled
// (now via a macro which records filename and line number).  Updated
// documentation in error.txt accordingly.
//
// Improved the documentation in matrix.txt (still more work to be done).
//
// Revision 1.5  2005/04/19 14:06:03  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.4  2005/03/11 16:44:18  cocoa
// New abstract class structure for matrices.
// New types of special matrix.
//
// Revision 1.3  2005/03/02 18:46:41  cocoa
// Added new types ConstRefMatrix, and RefMatrix following along
// the lines of ConstRefRingElem and RefRingElem.  The semantics
// should be a bit clearer now.
//
// Revision 1.2  2005/02/22 12:50:40  cocoa
// -- added: NotFullRank, NotIntDom, NotTermOrdering
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.17  2004/11/25 16:14:21  cocoa
// (1) Fixed definition of specialization of std::swap template function
//     so that it compiles with gcc 3.4.3
// (2) Implemented monomial function for polynomial rings.
// (3) Added one(PPM) and PPM->myOne() functions.
//
// Revision 1.16  2004/11/19 17:43:50  cocoa
// -- added matrix related errors
//
// Revision 1.15  2004/11/19 16:15:51  cocoa
// (a) Removed unused error message about degree of zero;
//     replaced it by a more general message about needing a
//     non-zero polynomial (for various ops such as LC, LPP).
// (b) Added some missing arg checking in LC, LPP and deg
//     (for elems of a PolyRing).
// (c) Updated some commented out code in GPair and GPoly.
//
// Revision 1.14  2004/11/19 15:14:09  cocoa
// (a) Added new check to MemPool so that it can signal an
//     error if one tries to use a MemPool after it has been
//     destroyed.
// (b) Improved makefile in TEST/ so that it checks output,
//     and prints useful messages if the test fails.
// (c) Tidied ring.txt a bit (still more to do).
//
// Revision 1.13  2004/11/12 15:49:29  cocoa
// Tidying prior to 0.90 release.
// (a) documentation improved (or marked as poor)
// (b) sundry minor improvements to the code
//
// Revision 1.12  2004/11/11 11:56:09  cocoa
// (1) Tidied makefiles, and introduced common.mki
// (2) Improved several tests, and cleaned them so that they
//     handle sanely any otherwise unhandled exceptions.
//
// Revision 1.11  2004/11/09 15:53:57  cocoa
// -- added BadQuotRing("Attempt to quotient by ideal(1)")
//
// Revision 1.10  2004/11/05 15:37:56  cocoa
// Added a couple of new errors.
//
// Revision 1.9  2004/11/03 18:01:55  cocoa
// -- clarified behaviour of error()
// -- added ANNOUNCE function
// -- minor tidying
//
// Revision 1.7  2004/07/27 16:03:38  cocoa
// Added IsCommutative test and IamCommutative member function
// to all rings.  Tidied geobuckets a little.
//
// Revision 1.6  2004/07/14 16:40:42  cocoa
// Separated RingFpLog from its implementation which now resides in
// a new class: SmallFpLogImpl.  This is analogous to the change made
// to RingFp yesterday.
//
// Some tidying and other sundry minor changes.
//
// Revision 1.5  2004/07/13 16:32:26  cocoa
// First stage of major revamp of ring elements.
// Implementation of RingFp has been split into "ring interface"
// and "algorithms plus data structures".
//
// Revision 1.4  2004/06/29 17:10:22  cocoa
// Partially tidied use of "protected" and "private" in various
// base classes.  Checking in at the end of the day -- it works,
// and I wouldn't want it to be lost next time point's disk
// misbehaves.
//
// Revision 1.3  2004/05/27 16:14:02  cocoa
// Minor revision for new coding conventions.
//
// Revision 1.2  2004/05/24 15:52:13  cocoa
// Major update:
//   new error mechanism
//   many fixes
//   RingHoms almost work now
//   RingFloat much improved
//
// Revision 1.1.1.1  2003/09/24 12:55:43  cocoa
// Imported files
//
// Revision 1.4  2003/06/23 16:15:21  abbott
// Minor cleaning prior to public release.
//
// Revision 1.3  2002/11/13 15:41:57  abbott
// The definition of CoCoA::error has been improved.  It now preceeds
// each output line with "CoCoA ERROR", it outputs on std::cerr
// (instead of std::cout), and it will print out the location information
// whenever present.  The function finally throws the exception whereas
// previously it just called exit -- throwing should be cleaner since
// stack unwinding occurs.
//
// Revision 1.2  2001/12/07 15:56:43  abbott
// error function now takes an exception as an argument, instead of a string.
//
// Revision 1.1  2001/11/26 19:13:33  abbott
// Initial revision
//
