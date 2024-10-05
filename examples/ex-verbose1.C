// Copyright (c) 2016  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This program illustrates the use of VerboseLog for producing progress  \n"
  "reports on the various steps of an algorithm.                          \n";


const string LongDescription =
  "This program illustrates the use of VerboseLog for producing progress    \n"
  "reports on the various steps of an algorithm.  You must set the global   \n"
  "verbosity level using the fn SetVerbosityLevel (higher values mean more  \n"
  "messages are printed) immediately before the call you want to investigate\n"
  "and set it back to 0 afterwards!                                         \n";

//----------------------------------------------------------------------

namespace CoCoA
{

  // To illustrate the use of verbose logging we take the following
  // "heuristic" function to look for short polynomials in an ideal.
  //
  // At the start we create a local object of type VerboseLog, we;
  // call it VERBOSE.  Later in the function to print out logging
  // info when VerbosityLevel is at least 20 we write a line like:
  //   VERBOSE(20) << "Logging message" << endl;
  vector<RingElem> ShortIdealElems(const ideal& I)
  {
    VerboseLog VERBOSE("ShortIdealElems"); // arg is name of this fn
    const ring& P = RingOf(I);
    const vector<RingElem>& GB = GBasis(I);
    VERBOSE(25) << "GBasis = " << GB << endl;
    const int n = len(GB);
    vector<RingElem> ans;
    // Try all pairs of GB elements to see if LCM of LPPs give a short poly under NF
    for (int i=0; i < n; ++i)
      for (int j=i+1; j < n; ++j)
      {
        VERBOSE(20) << "Doing (" << i << ", " << j << ")" << endl;
        const RingElem t = monomial(P,lcm(LPP(GB[i]), LPP(GB[j])));
        const RingElem f = NF(t, I);
        VERBOSE(25) << "NumTerms(NF) = " << NumTerms(f) << endl;
        if (NumTerms(f) < 10)
        {
          ans.push_back(t-f);
          VERBOSE(22) << "Found short poly " << ans.back() << endl;
        }
      }
    VERBOSE(20) << "Finished: found " << len(ans) << " short polys." << endl;
    return ans;
  }

  
  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;
    cout << boolalpha; // so that bools print out as true/false

    //-------------------------------------------------------
    // First example
    ring P = NewPolyRing(RingQQ(), symbols("x,y,z"));
    ideal I = ideal(RingElem(P, "x^2*z-y*z^2"),
                    RingElem(P, "z^4-1"),
                    RingElem(P, "y*z^3 - (x+y+z+1)^2"));

    SetVerbosityLevel(20); // change arg to get more/less verbose logging.
    const vector<RingElem> ShortPolysI = ShortIdealElems(I);
    SetVerbosityLevel(0); // set verbosity to 0 (i.e. no more logging)

    cout << "The ideal I = " << I << endl
         << "contains the following short polynomials:" << endl
         << ShortPolysI << endl << endl;


    //-------------------------------------------------------
    // Second example
    ideal J = ideal(RingElem(P, "x^15 -z^15 -2*z^14 -2*z^13 -2*z^12 -2*z^11"
                                "-2*z^10 -2*z^9 -2*z^8 -2*z^7 -2*z^6"
                                "-2*z^5 -2*z^4 -2*z^3 -2*z^2 -2*z -1"),
                    RingElem(P, "y^15 -z^15 +2*z^14 -2*z^13 +2*z^12 -2*z^11"
                                "+2*z^10 -2*z^9 +2*z^8 -2*z^7 +2*z^6"
                                "-2*z^5 +2*z^4 -2*z^3 +2*z^2 -2*z +1"));

    SetVerbosityLevel(22); // change arg to get more/less verbose logging.
    const vector<RingElem> ShortPolysJ = ShortIdealElems(J);
    SetVerbosityLevel(0); // set verbosity to 0 (i.e. no more logging)
    cout << "The ideal J = " << J << endl
         << "contains the following short polynomials:" << endl
         << ShortPolysJ << endl << endl;
    
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-verbose1.C,v 1.6 2022/02/13 09:57:01 abbott Exp $
// $Log: ex-verbose1.C,v $
// Revision 1.6  2022/02/13 09:57:01  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.5  2019/03/15 16:33:15  bigatti
// minor fix
//
// Revision 1.4  2017/09/25 12:36:21  abbott
// Summary: Added helpful comments
//
// Revision 1.3  2017/04/27 15:24:51  bigatti
// -- changed ReadExpr --> RingElem
//
// Revision 1.2  2016/11/11 15:39:43  abbott
// Summary: Made clearer (and also longer)
//
// Revision 1.1  2016/11/11 13:20:37  abbott
// Summary: Added new example ex-verbose1
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
