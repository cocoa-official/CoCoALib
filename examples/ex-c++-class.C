// Copyright (c) 2017  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This is a (harder) example showing the implemention of a simple C++ class.\n"
  "It shows class definition, and use of an object of that class.            \n"
  "The class contains a general factorization of an integer.                 \n";

const string LongDescription =
  "This is a (harder) example showing the implemention of a simple C++ class.\n"
  "It shows class definition, and use of an object of that class.            \n"
  "The class contains a general factorization of an integer.  The class data \n"
  "members are `private' so their values can be seen or changed only by      \n"
  "`member functions' or `friend functions'.                                 \n";
 
//----------------------------------------------------------------------

namespace CoCoA
{

  // A simple class for representing (general) factorizations of BigInts.
  // The class ensures that:
  // each factor is not -1, 0 or +1
  // each factor has (exactly) one multiplicity
  // each multiplicity is positive

  // EXERCISE:
  // can you modify it so that it also ensures that factors are prime?


  //---- CLASS DECLARATION ---
  // constructor/destructor,  member functions, friend functions, data members 

  class BigIntFactorization
  {
  public: // PUBLIC data/functions may be accessed by anyone
    //-- CONSTRUCTORs --
    BigIntFactorization(): myFactorVec(), myMultiplicityVec() { }

    //-- MEMBER FUNCTIONs --
    // read-only (const) access to data members:  [inline definition]
    const std::vector<BigInt>& myFactors() const { return myFactorVec; }
    const std::vector<long>& myMultiplicities() const { return myMultiplicityVec; }
    // safe modifier for data members:
    void myAppend(const BigInt& fac, long mult);

    //-- FRIEND FUNCTIONs (can access private data/functions) --
    // same as myProduct, but with functional syntax
    friend BigInt product(const BigIntFactorization& FacInfo);

  private: // PRIVATE data/functions accessed only by friend/member functions
    //-- DATA MEMBERs --
    std::vector<BigInt> myFactorVec;     // k-th entry contains k-th factor
    std::vector<long> myMultiplicityVec; // k-th entry contains multiplicity of k-th factor

    // PRIVATE functions can be called only by friend/member functions
    BigInt myProduct() const;
  };

  //---- FUNCTION DECLARATIONs (non-member functions: friends and other) ----
  void AppendFacPow(BigIntFactorization& FacInfo, const BigInt& fac, long mult);
  BigInt product(const BigIntFactorization& FacInfo);

  //----------------------------------------------------------------
  //---- FUNCTION DEFINITIONs (member, friend, other functions) ----

  void AppendFacPow(BigIntFactorization& FacInfo, const BigInt& fac, long mult)
  { FacInfo.myAppend(fac, mult); }


  BigInt product(const BigIntFactorization& FacInfo)
  { return FacInfo.myProduct(); }
  

  void BigIntFactorization::myAppend(const BigInt& fac, long mult)
  {
    if (abs(fac) < 2) CoCoA_THROW_ERROR("bad factor", "factorization::myAppend");
    if (mult <= 0) CoCoA_THROW_ERROR(ERR::NegExp, "factorization::myAppend");
    myFactorVec.push_back(fac);
    myMultiplicityVec.push_back(mult);
  }


  BigInt BigIntFactorization::myProduct() const // const: members not changed
  {
    BigInt ans(1);
    for (long i=0; i<len(myFactorVec); ++i)
      ans *= power(myFactorVec[i], myMultiplicityVec[i]);
    return ans;
  }

  // PRINTING: Function for printing out a BigIntFactorization.
  // [ CoCoALib has a function for printing vector<ANY-TYPE> ]
  std::ostream& operator<<(std::ostream& out, const BigIntFactorization& FacInfo)
  {
    return out << "BigIntFactorization(myFactors=" << FacInfo.myFactors()
               << ", myMultiplicities=" << FacInfo.myMultiplicities() << ")";
  }




  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    // The next lines create FacInfo, an object of type BigIntFactorization.
    // Then we "append" some factors to the factorization it represents.
    BigIntFactorization FacInfo;  // initially an "empty" factorization
    BigInt N(12345);

    AppendFacPow(FacInfo, N,3);         // "friend" function call syntax
    FacInfo.myAppend(BigInt(54321), 2); // "member" function call syntax

    cout << "FacInfo = " << FacInfo << endl;
    cout << "product(FacInfo) = " << product(FacInfo) << endl;
  }

} // end of namespace CoCoA

//----------------------------------------------------------------------
// Use main() to handle any uncaught exceptions and warn the user about them.
int main()
{
  try
  {
    CoCoA::program();
    return 0;
  }
  catch (const CoCoA::ErrorInfo& err)
  {
    cerr << "***ERROR***  UNCAUGHT CoCoA error";
    ANNOUNCE(cerr, err);
  }
  catch (const std::exception& exc)
  {
    cerr << "***ERROR***  UNCAUGHT std::exception: " << exc.what() << endl;
  }
  catch(...)
  {
    cerr << "***ERROR***  UNCAUGHT UNKNOWN EXCEPTION" << endl;
  }

  CoCoA::BuildInfo::PrintAll(cerr);
  return 1;
}

//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-c++-class.C,v 1.6 2022/02/13 09:57:00 abbott Exp $
// $Log: ex-c++-class.C,v $
// Revision 1.6  2022/02/13 09:57:00  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.5  2020/06/17 15:49:20  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.4  2017/07/08 19:07:02  abbott
// Summary: Removed comment out (dodgy) code for reporting unhandled interrupts
//
// Revision 1.3  2017/02/23 10:48:58  abbott
// Summary: Fixed a typo
//
// Revision 1.2  2017/02/23 08:32:49  bigatti
// -- many more comments
// -- added product/myProduct
// -- added EXERCISE
//
// Revision 1.1  2017/02/22 12:31:27  abbott
// Summary: New example: a simple C++ class
//
// Revision 1.11  2016/11/18 18:05:15  abbott
// Summary: Added commented out code to catch InterruptReceived
//
// Revision 1.10  2015/06/29 14:23:19  abbott
// Summary: added missing CoCoA:: prefix
// Author: JAA
//
// Revision 1.9  2015/06/29 13:25:54  bigatti
// -- code in namespace CoCoA
//
// Revision 1.8  2015/06/25 14:19:02  abbott
// Summary: Added call to CoCoA::BuildInfo::Printall
// Author: JAA
//
// Revision 1.7  2013/05/28 07:07:04  bigatti
// -- added "cout << boolalpha": useful for testing
//
// Revision 1.6  2012/11/30 14:04:55  abbott
// Increased visibility of comment saying "put your code here".
//
// Revision 1.5  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.4  2008/10/07 12:12:54  abbott
// Removed useless commented out #include.
//
// Revision 1.3  2007/05/31 16:06:16  bigatti
// -- removed previous unwanted checked-in version
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.9  2007/03/07 11:51:40  bigatti
// -- improved test alignment
//
// Revision 1.8  2007/03/03 14:15:45  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.7  2007/03/02 17:46:40  bigatti
// -- unique RingZ and RingQ
// -- requires foundations.H ;  foundations blah;  (thik of a better name)
//
// Revision 1.6  2007/03/02 10:47:53  cocoa
// First stage of RingZ modifications -- tests do not compile currently, Anna will fix this.
//
// Revision 1.5  2007/03/01 13:52:59  bigatti
// -- minor: fixed typo
//
// Revision 1.4  2007/02/28 15:15:56  bigatti
// -- minor: removed quotes in description
//
// Revision 1.3  2007/02/12 16:27:43  bigatti
// -- added strings ShortDescription and LongDescription for indexing
//
// Revision 1.2  2007/02/10 18:44:03  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.1.1.1  2006/05/30 11:39:36  cocoa
// Imported files
//
// Revision 1.1  2006/03/12 21:28:34  cocoa
// Major check in after many changes
//
