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
