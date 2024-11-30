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
