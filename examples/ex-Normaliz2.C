// Copyright (c) 2010  John Abbott, Anna Bigatti
// Orig authors: 2010  Anna Bigatti, Christof Soeger
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"
#include<fstream>
using std::ifstream;

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This example shows how to set some libnormaliz flags.";

const string LongDescription = 
  "This example reads a normaliz input file and computes the support hyperplanes.\n"
  "It also shows how to set some libnormaliz flags.";
//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
#ifndef CoCoA_WITH_NORMALIZ
    cout << "Normaliz library is not available to CoCoALib." << endl;
#else // NORMALIZ is available

    GlobalManager CoCoAFoundations;
    using namespace CoCoA::Normaliz;

    cout << ShortDescription << endl;
  
    //Set the output stream used to give libnormaliz warnings
    //(at the moment also error messages, but they will be transferred into the thrown Exception)
    libnormaliz::setErrorOutput(cerr);   //default: std::cerr
  
    //Set the output stream used to give libnormaliz verbose output;
    //verbose output will only be printed if the flag libnormaliz::verbose is set to true
    libnormaliz::setVerboseOutput(cout); //default: std::cout
    libnormaliz::verbose = true;         //default: false

    // Read normaliz input: num_rows, num_cols, then the matrix of integers row-by-row.
    long nr, nc;
    cin >> nr;
    cin >> nc;
    if (!cin || nr < 1 || nc < 1 || nr > 1000 || nc > 1000)
    { cerr << "Error: number of rows/cols must be between 1 and 1000.  Quitting!\n"; exit(1); }
    vector<vector<BigInt> > l(nr, vector<BigInt>(nc));

    for (long i=0; i<nr; i++)
      for (long j=0; j<nc; j++)
        cin >> l[i][j];

    cout << "The input is:\n" << l << endl;

    libnormaliz::Type::InputType InputType = libnormaliz::Type::integral_closure;

    cone c(InputType, l);
    long t0 = CpuTime();
    l = SupportHyperplanes(c);
    cout << "Time taken to compute support hyperplanes: " << CpuTime()-t0 << endl;
  
    cout << "The support hyperplanes are:\n" << l << endl;
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-Normaliz2.C,v 1.28 2022/02/14 11:25:49 bigatti Exp $
// $Log: ex-Normaliz2.C,v $
// Revision 1.28  2022/02/14 11:25:49  bigatti
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
// Revision 1.23  2014/07/11 11:07:52  abbott
// Summary: Corrected short/long description
// Author: JAA
//
// Revision 1.22  2014/03/21 11:17:09  abbott
// Summary: Replaced HilbertBasis by SupportHyperplanes (latter is much faster)
// Author: JAA
//
// Revision 1.21  2012/10/08 13:53:45  bigatti
// -- more cleaning and updates by Christof Soeger
//
// Revision 1.20  2012/07/19 17:11:13  abbott
// Added unified "NewCone" (user does not have to choose between long or BigInt).
// Cleaned up.
//
// Revision 1.19  2012/07/04 12:32:04  abbott
// Changed int into long; improved ShortDescription & LongDescription.
//
// Revision 1.18  2012/04/11 09:55:04  abbott
// Changed print stream to cout for mesg when Normliz is absent
// (previously was clog, but this caused the script to think the example had failed).
//
// Revision 1.17  2012/03/16 15:44:29  abbott
// Cleaned up main reading loop.
//
// Revision 1.16  2011/11/07 11:30:47  bigatti
// -- changed syntax for HilbertBasis
//
// Revision 1.15  2011/10/12 15:50:17  abbott
// Simplified use of CPP macro CoCoA_WITH_NORMALIZ.
//
// Revision 1.14  2011/09/30 16:08:45  abbott
// *** empty log message ***
//
// Revision 1.13  2011/09/30 15:59:24  bigatti
// -- moved "using namespace CoCoA::Normaliz" inside #ifdef
//
// Revision 1.12  2011/09/30 12:55:44  bigatti
// -- introduced namespace "Normaliz" and removed Normaliz from function names
// -- input of Normaliz functions in CoCoA-5 is now a matrix instead of
//    vector<vector<BigInt>>
//
// Revision 1.11  2011/08/23 12:04:04  bigatti
// -- updated after renaming ZZ --> BigInt
//
// Revision 1.10  2011/07/20 13:49:37  bigatti
// -- added "Normaliz" postfix to Normaliz function calls
//
// Revision 1.9  2011/07/20 12:45:12  bigatti
// -- new normaliz interface (not yet public)
//
// Revision 1.8  2011/02/17 17:12:57  bigatti
// -- reading input from cin
//
// Revision 1.7  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.6  2010/10/22 13:31:01  bigatti
// -- added comments about normaliz flags (by C.Soeger)
//
// Revision 1.5  2010/10/12 11:22:54  bigatti
// -- TmpNormaliz.H simplified:
//    now NormalizCone is a smart pointer to NormalizConeBase
//    and the concrete classes are entirely in the .C file
// -- added Ht1Elements, SupportHyperplanes, Triangulation
// -- added some text in ex-Normaliz1 and 2
//
// Revision 1.4  2010/10/08 16:33:27  bigatti
// -- added #include <fstream>
//
// Revision 1.3  2010/10/08 13:53:51  bigatti
// -- cleaned up
//
// Revision 1.2  2010/10/08 10:39:32  bigatti
// -- extended interface for normaliz
//
// Revision 1.1  2010/10/08 06:57:09  bigatti
// * ex-Normaliz2-1.in (Module): first import
//
// * ex-Normaliz2-2.in (Module): first import
//
// * ex-Normaliz2.C (Module): first import
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
