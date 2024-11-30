// Copyright (c) 2017  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This is a (harder) example showing the implementation of a simple C++ class.\n"
  "It shows class definition, and use of an object of that class.            \n"
  "The class contains a general factorization of an integer.                 \n";

const string LongDescription =
  "This is a (harder) example showing the implementation of a simple C++ class.\n"
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
