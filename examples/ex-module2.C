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

//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-module2.C,v 1.9 2022/02/13 09:57:01 abbott Exp $
// $Log: ex-module2.C,v $
// Revision 1.9  2022/02/13 09:57:01  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.8  2018/09/28 15:54:04  abbott
// Summary: Removed pseudo-ctors NewPolyRing which took just num of indets; now must specify their names
//
// Revision 1.7  2015/06/26 15:34:48  abbott
// Summary: Moved code into namespace CoCoA (see redmine 739)
// Author: JAA
//
// Revision 1.6  2014/07/30 13:59:40  abbott
// Summary: Changed BaseRing into RingOf
// Author: JAA
//
// Revision 1.5  2014/07/09 14:28:46  abbott
// Summary: Removed AsFreeModule and AsFGModule
// Author: JAA
//
// Revision 1.4  2014/07/08 12:47:26  abbott
// Summary: Removed AsPolyRing, AsSparsePolyRing, AsQuotientRing
// Author: JAA
//
// Revision 1.3  2013/08/02 16:40:02  bigatti
// -- LPos --> LPosn
//
// Revision 1.2  2013/06/06 05:43:08  bigatti
// -- added IsElem, IsContained, TidyGens
//
// Revision 1.1  2013/05/28 07:08:17  bigatti
// -- first import
//
// Revision 1.6  2013/01/23 14:43:36  bigatti
// -- added "v==w"
//
// Revision 1.5  2012/02/08 17:52:17  bigatti
// -- changed: Z,Q -> ZZ,QQ
//
// Revision 1.4  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.3  2008/04/21 12:33:45  abbott
// Made a very minor change.
//
// Revision 1.2  2007/06/21 21:29:47  abbott
// Changed name of RingFloat into RingTwinFloat.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.6  2007/03/03 14:15:45  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.5  2007/02/28 15:15:10  bigatti
// -- added sum and multiplication examples
//
// Revision 1.4  2007/02/26 17:40:34  bigatti
// -- getting ready for unique ring Z: using NewZmod(N), NewRingQ()
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
// Revision 1.1  2005/12/16 17:53:01  cocoa
// --- first import
//
