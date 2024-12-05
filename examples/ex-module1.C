// Copyright (c) 2005  John Abbott,  Anna Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "Example program illustrating the creation of free modules, and \n"
  "some operations on them. \n";

const string LongDescription =
  "Please note that the module code is still rather young. \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void trial(ring R)
  {
    const int NumComps = 4;
    const FreeModule F = NewFreeModule(R, NumComps);
    const vector<ModuleElem>& e = gens(F);

    const ModuleElem u = e[0] + 2*e[1];
    const ModuleElem v = 4*e[0] + 3*e[3];
    const RingElem a = 2 * one(RingOf(F));  // RingOf(F) is R

    cout << "---- F = " << F << " ----" << endl;
    cout << "u[0] = " << u[0] << ";  u[1] = " << u[1] << endl;
    cout << "u + v = " << u + v << endl;
    cout << "u * a = " << u * a << endl;
    cout << "a * v = " << a * v << endl;

    cout << "u == v is " << (u == v) << endl;
  }


  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;
    cout << boolalpha;

    trial(RingZZ());
    trial(RingQQ());
    trial(NewRingTwinFloat(32));
    trial(NewZZmod(2));
    trial(NewZZmod(1048576)); // ring has zero divisors
  }

} // end of namespace CoCoA


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
