// Copyright (c) 2005,2007  John Abbott, Anna Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "Example of use of power products and PPMonoids.     \n"
  "Program exhibiting most functions on power products.\n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    PPMonoid PPM = NewPPMonoidEv(SymbolRange("x",0,3), lex);
    // Indets will be printed are x[0], x[1], x[2] & x[3].

    // We can ask the PPM to create powers of its own indets...
    PPMonoidElem t1 = indet(PPM,0) * IndetPower(PPM,1,3) * IndetPower(PPM,2,7);

    // or, for handy access to the indeterminates in PPM...
    const vector<PPMonoidElem>& x = indets(PPM);

    PPMonoidElem t2 = x[0] * power(x[1],7) * power(x[2],3);

    cout << "initial power products are\nt1 = " << t1 << "\nt2 = " << t2 << endl;
    cout << endl;

    cout << "t1 == t2 gives " << (t1 == t2) << endl;
    cout << "t1 != t2 gives " << (t1 != t2) << endl;
    cout << "t1 < t2  gives " << (t1 < t2) << endl;
    cout << "t1 <= t2 gives " << (t1 <= t2) << endl;
    cout << "t1 > t2  gives " << (t1 > t2) << endl;
    cout << "t1 >= t2 gives " << (t1 >= t2) << endl;
    cout << endl;

    cout << "t1 * t2  gives " << t1*t2 << endl;
    cout << "IsDivisible(t1, t2) gives " << IsDivisible(t1, t2) << endl;
    cout << "We CANNOT compute t2 / t1" << endl;
    cout << endl;

    cout << "colon(t1, t2) gives " << colon(t1, t2) << endl;
    cout << "colon(t2, t1) gives " << colon(t2, t1) << endl;
    cout << "gcd(t1, t2)   gives " << gcd(t1, t2) << endl;
    cout << "lcm(t1, t2)   gives " << lcm(t1, t2) << endl;
    cout << "power(t1, 5)  gives " << power(t1, 5) << endl;
    cout << "IsCoprime(t1, t2) gives " << IsCoprime(t1, t2) << endl;
    cout << endl;

    cout << "StdDeg(t1) gives " << StdDeg(t1) << endl;
    cout << "wdeg(t1) gives " << wdeg(t1) << endl;
    cout << "[note: Lex is automatically ungraded (i.e. graded over Z^0)]" << endl;
    cout << "[see: ex-OrderingGrading1.C]" << endl;
    cout << endl;

    cout << "We can find the power to which each indeterminate appears:" << endl;
    cout << "exponent(t1, 0) gives " << exponent(t1, 0) << endl;
    cout << "exponent(t1, 1) gives " << exponent(t1, 1) << endl;
    cout << "exponent(t1, 2) gives " << exponent(t1, 2) << endl;
    cout << "exponent(t1, 3) gives " << exponent(t1, 3) << endl;
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
