// Copyright (c) 2010  John Abbott,  Anna Bigatti
// Orig authors: 2010  Anna Bigatti, Christof Soeger
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This program shows what can be computed in CoCoALib using Normaliz:  \n"
  "a library for computations in affine monoids, vector configurations, \n"
  "lattice polytopes, and rational cones.";

const string LongDescription = 
  "This program shows what can be computed in CoCoALib using Normaliz:  \n"
  "a library for computations in affine monoids, vector configurations, \n"
  "lattice polytopes, and rational cones.";

//----------------------------------------------------------------------

namespace CoCoA
{

  // utility function while waiting for C++Ox vector initializer
  vector<BigInt> BigIntVec(long* CVector, long len)
  {
    vector<BigInt> v;
    for (long i=0; i<len; ++i)  v.push_back(BigInt(CVector[i]));
    return v;
  }


  void program()
  {
    GlobalManager CoCoAFoundations;

#ifndef CoCoA_WITH_NORMALIZ
    cout << "Normaliz library is not available to CoCoALib." << endl;
#else // NORMALIZ is available
    cout << ShortDescription << endl;
    cout << boolalpha; // so that bools print out as true/false

    long M[5][3] = {{2, 0, 0},
                    {1, 3, 1},
                    {3, 3, 1},
                    {0, 0, 2},
                    {1,-3, 1}};

    vector<vector<BigInt> > l;
    for (long i=0; i<5; ++i)  l.push_back(BigIntVec(M[i], 3));

    cout << "l -> " << l << endl;

    // waiting for Normaliz Team decisions on input
    // namespace Type {
    // enum InputType {
    // 	integral_closure,
    // 	normalization,
    // 	polytope,
    // 	rees_algebra,
    // 	lattice_ideal
    // 	hyperplanes,
    // 	equations,
    // 	congruences
    // };
    // } //end namespace Type

    //  libnormaliz::verbose = true;         //default: false
    libnormaliz::Type::InputType InputType = libnormaliz::Type::normalization;
    vector<vector<BigInt> > res;
    double t0;

    Normaliz::cone C3(InputType, l);
    cout << "Cone created: " << endl << C3 << endl << endl;
  
    t0 = CpuTime();
    res = Normaliz::HilbertBasis(C3);
    cout << "time HilbertBasis  -> " << CpuTime()-t0 << endl;
    cout << "HilbertBasis -> " << res << endl;

    cout << endl;
    cout << "More operations on Cone:" << endl;
    res = Normaliz::Deg1Elements(C3);
    cout << "Deg1Elements -> " << res << endl;
    res = Normaliz::SupportHyperplanes(C3);
    cout << "SupportHyperplanes -> " << res << endl;

    cout << "Hilbert series -> " << Normaliz::HilbertSeries(C3) << endl;
    //   Normaliz::Triangulation(res, c2);
    //   cout << "Triangulation -> " << res << endl;
    cout << "Cone is now: " << endl << C3 << endl << endl;

#endif // CoCoA_WITH_NORMALIZ
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-Normaliz1.C,v 1.28 2022/02/14 11:25:48 bigatti Exp $
// $Log: ex-Normaliz1.C,v $
// Revision 1.28  2022/02/14 11:25:48  bigatti
// Summary: added date for original authors
//
// Revision 1.27  2022/02/13 09:56:57  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.26  2015/09/03 10:17:11  bigatti
// -- changes by Christof Soeger (in Aarhus)
//
// Revision 1.25  2015/06/29 15:27:25  bigatti
// -- minor fix
//
// Revision 1.24  2015/06/29 12:20:43  bigatti
// -- moved in cocoa namespace
//
// Revision 1.23  2014/07/11 10:11:53  abbott
// Summary: Christof's cleaning of the example
// Author: JAA
//
// Revision 1.22  2014/01/30 09:55:55  bigatti
// -- removed DOS newlines
//
// Revision 1.21  2014/01/29 17:34:22  bigatti
// -- added HilbertSeries (by Christof Soeger)
//
// Revision 1.20  2012/10/08 13:53:45  bigatti
// -- more cleaning and updates by Christof Soeger
//
// Revision 1.19  2012/08/03 16:31:22  bigatti
// -- changed: procedural --> functional (by C.Soeger)
//
// Revision 1.18  2012/07/20 13:27:33  abbott
// Commented out the code which uses the (obsolescent) pseudo-ctors NewConeLong
// and NewConeBigInt.
//
// Revision 1.17  2012/07/19 17:11:13  abbott
// Added unified "NewCone" (user does not have to choose between long or BigInt).
// Cleaned up.
//
// Revision 1.16  2012/07/04 12:31:08  abbott
// Changed int into long; improved ShortDescription & LongDescription.
//
// Revision 1.15  2012/06/19 14:39:39  bigatti
// -- unix-ified
//
// Revision 1.14  2012/06/19 14:36:23  bigatti
// -- by C.Soeger
//
// Revision 1.13  2011/11/09 15:13:49  bigatti
// -- new syntax for HilbertBasis
//
// Revision 1.12  2011/10/12 15:50:17  abbott
// Simplified use of CPP macro CoCoA_WITH_NORMALIZ.
//
// Revision 1.11  2011/09/30 12:55:44  bigatti
// -- introduced namespace "Normaliz" and removed Normaliz from function names
// -- input of Normaliz functions in CoCoA-5 is now a matrix instead of
//    vector<vector<BigInt>>
//
// Revision 1.10  2011/08/23 12:04:04  bigatti
// -- updated after renaming ZZ --> BigInt
//
// Revision 1.9  2011/07/20 13:49:37  bigatti
// -- added "Normaliz" postfix to Normaliz function calls
//
// Revision 1.8  2011/07/20 12:45:12  bigatti
// -- new normaliz interface (not yet public)
//
// Revision 1.7  2011/02/17 16:50:04  bigatti
// -- getting ready for new official veson on Normaliz: added HVector, removed Triangulation
//
// Revision 1.6  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.5  2010/11/22 17:33:31  bigatti
// -- input with c-style array
//
// Revision 1.4  2010/10/12 11:22:54  bigatti
// -- TmpNormaliz.H simplified:
//    now Cone is a smart pointer to ConeBase
//    and the concrete classes are entirely in the .C file
// -- added Ht1Elements, SupportHyperplanes, Triangulation
// -- added some text in ex-Normaliz1 and 2
//
// Revision 1.3  2010/10/08 13:53:51  bigatti
// -- cleaned up
//
// Revision 1.2  2010/10/08 10:39:32  bigatti
// -- extended interface for normaliz
//
// Revision 1.1  2010/10/05 14:31:26  bigatti
// -- first import
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
