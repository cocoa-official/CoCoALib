// Copyright (c) 2008  John Abbott,  Anna Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "Example of use of power products in different PPMonoids.   \n"
  "Program exhibiting timings of the different implementations.\n";
  const string LongDescription =
    "The implementations of PPMonoids are optimized for different uses:  \n"
    "PPMonoidEv:   stores the Exponent vector                            \n"
    "              good for accessing the exponents, slow for ordering   \n"
    "PPMonoidOv:   stores the Order vector                               \n"
    "              good for ordering, slow for accessing the exponents   \n"
    "PPMonoidEvOv: stores the Exponent vector and the Order vector       \n"
    "              good for accessing the exponents and for ordering     \n"
    "              uses more memory and takes more time to assign        \n"
    "PPMonoidBigEv: stores the Exponent vector as BigInt's               \n"
    "              necessary if you use big exponents (>2^10)            \n"
    "              obviously slow   \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void TryPPMonoid(const PPMonoid PPM)
  {
    //  cout << PPM << endl;
    const vector<PPMonoidElem>& x = indets(PPM);

    // at least 20 indets!
    PPMonoidElem t1 = product(x);
    PPMonoidElem t2 = x[0] * x[1] * x[20];

    double t0;  // for CpuTime

    constexpr long NumIters = 50000;

    // If your optimizer is too clever it might realise these operations
    // are "useless" (return value not used and no side effect)
    // so you will get constant timings
    t0 = CpuTime();
    for (long i=1; i < NumIters; ++i)    t2 < t1;
    cout << " <    -->          " << CpuTime()-t0 << " s" << endl;
    t0 = CpuTime();
    for (long i=1; i < NumIters; ++i)    t1 * t2;
    cout << " *    -->          " << CpuTime()-t0 << " s" << endl;
    t0 = CpuTime();
    for (long i=1; i < NumIters; ++i)    gcd(t1, t2);
    cout << " gcd  -->          " << CpuTime()-t0 << " s" << endl;  
    for (long i=1; i < NumIters; ++i)    IsDivisible(t1, t2);
    cout << " IsDivisible  -->  " << CpuTime()-t0 << " s" << endl;  
  }


  void TryOrd(const PPOrdering ord)
  {
    const vector<symbol> x = SymbolRange("x", 0, NumIndets(ord)-1);
    cout << "-- =========================================== --" << endl;
    cout << "   " << ord << endl;
    cout << "-- =========================================== --" << endl;
    cout << "---------   PPMonoidBigEv  ----------" << endl;
    TryPPMonoid(NewPPMonoidEv(x, ord, PPExpSize::big));
    cout << "----------   PPMonoidEv  -----------" << endl;
    TryPPMonoid(NewPPMonoidEv(x, ord));
    cout << "----------   PPMonoidOv  -----------" << endl;
    TryPPMonoid(NewPPMonoidOv(x, ord));
    cout << "---------   PPMonoidEvOv  ----------" << endl;
    TryPPMonoid(NewPPMonoidEvOv(x, ord));
  }


  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;
    constexpr int n = 40;
    vector<long> y(1,1);  // y = [1]
    PPOrdering MatOrd = NewMatrixOrdering(ElimMat(y, n), 1);

    TryOrd(lex(n));
    TryOrd(StdDegRevLex(n));
    TryOrd(MatOrd);
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
