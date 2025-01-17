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
