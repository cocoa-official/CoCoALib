// Copyright (c) 2005,2012  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "When computing over a finite field normally it is best to use the\n"
  "function NewZZmod to create the field.  However, for the curious \n"
  "it is possible to create small prime finite fields stating the   \n"
  "particular implementation method (there are 3 possibilities).      \n";

const string LongDescription =
  "This program compares the speeds of computing sums and products in   \n"
  "the three different implementations of small prime finite fields.    \n"
  "It performs the timing tests for different sizes of prime, and       \n"
  "illustrates that different implementations have different upper      \n"
  "limits for the characteristic -- these limits depend on the platform.\n"
  "It also gives the performance of the finite field created by NewZZmod\n";

//----------------------------------------------------------------------

namespace CoCoA
{

  // Measure how long it takes to build an NxN addition table.
  double TimeAdd(ring R, long N)
  {
    vector<RingElem> value;
    value.reserve(N); // not necessary but good style: reserves memory (C++ STL)
    for (long i=1; i <= N; ++i)
      value.push_back(RingElem(R, i));

    RingElem tmp(R);
    const double StartTime = CpuTime();
    for (long i=0; i < N; ++i)
      for (long j=0; j < N; ++j)
        tmp = value[i] + value[j];
    return CpuTime()-StartTime;
  }

  // Measure how long it takes to build an NxN multiplication table.
  double TimeMult(ring R, long N)
  {
    vector<RingElem> value; value.reserve(N);
    for (long i=1; i <= N; ++i)
      value.push_back(RingElem(R, i));

    RingElem tmp(R);
    const double StartTime = CpuTime();
    for (long i=0; i < N; ++i)
      for (long j=0; j < N; ++j)
        tmp = value[i] * value[j];
    return CpuTime()-StartTime;
  }


  void test(long p)
  {
    cout << "Timing test for characteristic " << p << endl;

    constexpr long N=1000;

    cout << "Time to compute addition table of size " << N << endl;

    cout << "   ZZmod:        " << TimeAdd(NewZZmod(p), N) << endl;

    cout << "   RingFp:       ";
    try { cout << TimeAdd(NewRingFp(p), N) << endl; }
    catch (...) { cout << "CHAR TOO BIG\n"; }

    cout << "   RingFpLog:    ";
    try { cout << TimeAdd(NewRingFpLog(p), N) << endl; }
    catch (...) { cout << "CHAR TOO BIG\n"; }

    cout << "   RingFpDouble: ";
    try { cout << TimeAdd(NewRingFpDouble(p), N) << endl; }
    catch (...) { cout << "CHAR TOO BIG\n"; }

    cout << endl;

    cout << "Time to compute multiplication table of size " << N << endl;

    cout << "   ZZmod:        " << TimeMult(NewZZmod(p), N) << endl;

    cout << "   RingFp:       ";
    try { cout << TimeMult(NewRingFp(p), N) << endl; }
    catch (...) { cout << "CHAR TOO BIG\n"; }

    cout << "   RingFpLog:    ";
    try { cout << TimeMult(NewRingFpLog(p), N) << endl; }
    catch (...) { cout << "CHAR TOO BIG\n"; }

    cout << "   RingFpDouble: ";
    try { cout << TimeMult(NewRingFpDouble(p), N) << endl; }
    catch (...) { cout << "CHAR TOO BIG\n"; }

    cout << "--------------------------------------------------" << endl;

  }


  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    for (long P = 256; P < power(2,30); P *= 16)
    {
      test(PrevPrime(P));
      if (P > numeric_limits<long>::max()/16) break;
    }
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
