// Copyright (c) 2005  John Abbott,  Anna Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "Example program illustrating the creation of polynomial modules\n"
  "with ordering and shifts. \n";

const string LongDescription =
  "Please note that the module code is still rather young. \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void trial(FreeModule F)
  {
    const vector<ModuleElem>& e = gens(F);
    const SparsePolyRing P = RingOf(F);
    const vector<RingElem>& x = indets(P);

    const ModuleElem u = power(x[1],5)*e[0] + x[0]*e[1];
    const ModuleElem v = x[2]*x[2]*e[1] + x[1]*x[2]*e[2];

    cout << "---- F = " << F << " ----" << endl;
    cout << "-- " << ordering(F) << endl;
    //  cout << "-- " << ordering(PPM(P)) << endl;
    cout << "u[0] = " << u[0] << "  and  u[1] = " << u[1] << endl;
    cout << "u = " << u << endl;
    cout << "LPosn(u) = " << LPosn(u) << std::endl;
    cout << "LPP(u) = " << LPP(u) << endl;
    cout << "wdeg(u) = " << wdeg(u) << endl;
    if (GradingDim(P)==0)
      cout << "IsHomog(u) undefined because GradingDim is 0" << endl;
    else
      cout << "IsHomog(u) = " << IsHomog(u) << std::endl;
    cout << "v = " << v << endl;
    cout << "LPosn(v) = " << LPosn(v) << std::endl;

    cout << endl;

    // submodules
    module M = submodule(u-v, u+v);
    module N = submodule(u);
    cout << "M = "  << M << endl;
    cout << "IsElem(u, M) = "  << IsElem(u, M) << endl;
    cout << "IsContained(N, M) = "  << IsContained(N, M) << endl;
    cout << "TidyGens(M) = "  << TidyGens(M) << endl;
  }


  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;
    cout << boolalpha;

    ring QQ = RingQQ();
    SparsePolyRing PLex    = NewPolyRing(QQ, symbols("x,y,z"), lex);
    SparsePolyRing PDegLex = NewPolyRing(QQ, symbols("x,y,z"), StdDegLex);
    long n = 4;
    trial(NewFreeModule(PLex, n));
    trial(NewFreeModule(PDegLex, n));
    //  trial(NewFreeModule(PDegLex, n, PosnOrd)); // Not Yet Implemented
    trial(NewFreeModule(PDegLex,  n, OrdPosn));
    trial(NewFreeModule(PLex,     n, WDegPosnOrd));
    trial(NewFreeModule(PDegLex,  n, WDegPosnOrd));
    // with shifts
    std::vector<degree> sh(n, wdeg(one(PDegLex)));
    sh[1] = wdeg(power(indet(PDegLex,0),4));
    trial(NewFreeModule(PDegLex, sh, WDegPosnOrd));
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
