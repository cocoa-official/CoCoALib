// Copyright (c) 2017  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This is a more advanced example showing the implementation of a simple \n"
  "C++ class.  It illustrates class definition, and use of an object of   \n"
  "that class.  The example class represents a general factorization of   \n"
  "an integer (of unlimited size).                                        \n";

const string LongDescription =
  "This is a more advanced example showing the implementation of a simple \n"
  "C++ class.  It illustrates class definition, and use of an object of   \n"
  "that class.  The example class represents a general factorization of   \n"
  "an integer (of unlimited size).  Following standard recommendations,   \n"
  "the class data members are `private' so their values can be seen or    \n"
  "changed only by `member functions' or `friend functions'.              \n";
 
//----------------------------------------------------------------------

namespace CoCoA
{

  // MORE ADVANCED EXAMPLE: defining a C++ class
  // (this is just a refresher, not a complete tutorial)

  // We define a simple class for representing factorizations of BigInts.
  // The class ensures that:
  // (*) each factor is not -1, 0 or +1;
  // (*) each multiplicity is positive.

  // EXERCISE:
  //   Can you modify the class so that it also ensures that the factors
  //   are pairwise coprime or even distinct and (probably) prime?
  // HINT: you need to modify just 1 function.


  //---- CLASS DECLARATION ---
  // constructor/destructor,  member functions, friend functions, data members 

  class BigIntFactorization
  {
  public: // PUBLIC data/functions may be accessed by anyone
    //-- CONSTRUCTORs --
    BigIntFactorization(): myFactorVec(), myMultiplicityVec() { }
    BigIntFactorization(const BigIntFactorization&) = default;

    //-- ASSIGNMENT--
    BigIntFactorization& operator=(const BigIntFactorization&) = default;

    //-- MEMBER FUNCTIONs --
    // Read-only (const) accessors to data members:  [inline definition]
    const std::vector<BigInt>& myFactors() const  { return myFactorVec; }
    const std::vector<long>& myMultiplicities() const  { return myMultiplicityVec; }
    // Safe modifier for data members:
    void myAppend(const BigInt& fac, long mult);
    // Compute some "property" of the object:
    BigInt myProduct() const;

  private: // PRIVATE data/functions accessed only by friend/member functions
    //-- PRIVATE DATA MEMBERs --
    std::vector<BigInt> myFactorVec;     // k-th entry contains k-th factor
    std::vector<long> myMultiplicityVec; // k-th entry contains multiplicity of k-th factor

    // -- PRIVATE MEMBER FUNCTIONS --
    // No private member functions
  };

  //---- FUNCTION DECLARATIONs (non-member functions: friends and other) ----
  void AppendFacPow(BigIntFactorization& FacInfo, const BigInt& fac, long mult);
  BigInt product(const BigIntFactorization& FacInfo);


  //----------------------------------------------------------------
  //---- FUNCTION DEFINITIONs (member, friend, other functions) ----

  void BigIntFactorization::myAppend(const BigInt& fac, long mult)
  {
    if (abs(fac) < 2)
      CoCoA_THROW_ERROR2(ERR::BadArg, "factor must not be -1, 0 or +1");
    if (mult <= 0)
      CoCoA_THROW_ERROR1(ERR::NegExp);
    myFactorVec.push_back(fac);
    myMultiplicityVec.push_back(mult);
  }

  // Alternative syntax -- delegate argument checks to the member function
  void AppendFacPow(BigIntFactorization& FacInfo, const BigInt& fac, long mult)
  { FacInfo.myAppend(fac, mult); }



  BigInt BigIntFactorization::myProduct() const // const: members not changed
  {
    BigInt ans(1);
    for (long i=0; i<len(myFactorVec); ++i)
      ans *= power(myFactorVec[i], myMultiplicityVec[i]);
    return ans;
  }

  // Alternative syntax
  BigInt product(const BigIntFactorization& FacInfo)
  { return FacInfo.myProduct(); }
  

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

    AppendFacPow(FacInfo, N,3);         // standard function call syntax
    FacInfo.myAppend(BigInt(54321), 2); // member function call syntax

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
