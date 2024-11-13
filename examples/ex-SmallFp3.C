// Copyright (c) 2013-2015  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "**ADVANCED**  This example program is for advanced CoCoALib users. \n"
  "It shows how to use SmallFpImpl for efficient arithmetic in a small\n"
  "prime finite field.  Not for the faint-hearted!                    \n";

const string LongDescription =
  "**ADVANCED**  This example program is for advanced CoCoALib users.    \n"
  "SmallFpImpl enables you to perform arithmetic efficiently in a small  \n"
  "prime finite field.  The catch is that efficient use is not as simple \n"
  "as using RingElems directly.  We take as a specific illustrative      \n"
  "example the computation of an inner product.  Be prepared to spend    \n"
  "some time reading and comprehending the code!                         \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  // The impl of inner product for RingElem is pretty simple and clear:
  RingElem InnerProd_RINGELEM(const vector<RingElem>& u, const vector<RingElem>& v)
  {
    ring Fp = owner(u[0]);
    const long n = len(u);
    RingElem ans(Fp);
    for (long i=0; i < n; ++i)
      ans += u[i]*v[i];
    return ans;
  }


  // Now some faster functions using SmallFpImpl -- faster but less clear.
  // Handy typedef to make reading/writing the code simpler!!!
  typedef SmallFpImpl::value FpElem;


  // Here is a SIMPLE BUT INEFFICIENT impl using SmallFpImpl
  // (because every intermediate value is normalized).
  FpElem InnerProd_SLOW(const SmallFpImpl& Fp, const vector<FpElem>& u, const vector<FpElem>& v)
  {
    const long n = len(u);
    FpElem ans; // initially zero
    for (long i=0; i < n; ++i)
      ans = Fp.myAdd(ans, Fp.myMul(u[i],v[i]));
    return ans;
  }


  // We present 2 fast impls (which avoid many intermediate reductions):
  //   (A) is slightly clearer, while
  //   (B) is slightly faster.
  // You decide whether the extra complication of impl (B) is worth the speed gain!

  // Impl (A) for "fast" inner product;
  // it is much fiddlier than the RingElem implementation above!
  FpElem InnerProd_A(const SmallFpImpl& Fp, const vector<FpElem>& u, const vector<FpElem>& v)
  {
    const long n = len(u);
    CoCoA_ASSERT(len(v) == n);
    const long MaxSafeIters = Fp.myMaxIters();

    SmallFpImpl::NonRedValue ans; // initially zero
    long NextHalfNormalize = MaxSafeIters;
    for (long i=0; i < n; ++i)
    {
      ans += u[i]*v[i];  // <--- the actual computation, the rest is overhead!
      if (i < NextHalfNormalize) continue;
      NextHalfNormalize += MaxSafeIters; // overflow?
      ans = Fp.myHalfNormalize(ans);
    }
    return Fp.myNormalize(ans);
  }

  // Impl (B) for "fast" inner product;
  // it is harder to understand than (A), but is a bit faster (on my computer).
  FpElem InnerProd_B(const SmallFpImpl& Fp, const vector<FpElem>& u, const vector<FpElem>& v)
  {
    const long n = len(u);
    const long MaxSafeIters = Fp.myMaxIters();

    long i = 0; // loop counter
    long NextNormalize = 0;
    SmallFpImpl::NonRedValue ans;
    while (NextNormalize < n)
    {
      NextNormalize += MaxSafeIters;
      if (NextNormalize > n) NextNormalize = n;
      for (; i < NextNormalize; ++i)
        ans += u[i]*v[i];  // <--- the actual computation, the rest is overhead!
      ans = Fp.myHalfNormalize(ans);
    }
    return Fp.myNormalize(ans);
  }


  ///////////////////////////////////////////////////////
  void TimeTrialRingElem(long p)
  {
    ring Fp = NewZZmod(p);
    vector<RingElem> u; u.reserve(p*p);
    vector<RingElem> v; v.reserve(p*p);
    for (long i=0; i < p; ++i)
      for (long j=0; j < p; ++j)
      {
        u.push_back(RingElem(Fp,i));
        v.push_back(RingElem(Fp,j));
      }

    // Timing the computation of the inner product:
    const double StartTime = CpuTime();
    const RingElem InProd = InnerProd_RINGELEM(u, v);
    const double EndTime = CpuTime();
    cout << "Ans is " << InProd << endl;
    cout << "Using ring ZZmod(" << p << ") time is " << EndTime-StartTime << endl;
  }


  void TimeTrialSmallFp(long p)
  {
    if (!SmallFpImpl::IsGoodCtorArg(p)) return;//????
    SmallFpImpl Fp(p);
    // Create two vectors to work on
    vector<FpElem> u; u.reserve(p*p);
    vector<FpElem> v; v.reserve(p*p);
    for (long i=0; i < p; ++i)
      for (long j=0; j < p; ++j)
      {
        u.push_back(Fp.myReduce(i));
        v.push_back(Fp.myReduce(j));
      }

    // Timing method SLOW (normalize every intermediate result)
    const double StartTime_SLOW = CpuTime();
    const FpElem InProd_SLOW = InnerProd_SLOW(Fp, u, v);
    const double EndTime_SLOW = CpuTime();
    cout << "Ans is " << Fp.myExport(InProd_SLOW) << endl;
    cout << "Using impl (SLOW) for p=" << p << "  time is " << (EndTime_SLOW - StartTime_SLOW) << endl;

    // Timing method (A)
    const double StartTime_A = CpuTime();
    const FpElem InProd_A = InnerProd_A(Fp, u, v);
    const double EndTime_A = CpuTime();
    cout << "Ans is " << Fp.myExport(InProd_A) << endl;
    cout << "Using impl (A) for p=" << p << "  time is " << (EndTime_A - StartTime_A) << endl;

    // Timing method (B)
    const double StartTime_B = CpuTime();
    const FpElem InProd_B = InnerProd_B(Fp, u, v);
    const double EndTime_B = CpuTime();
    cout << "Ans is " << Fp.myExport(InProd_B) << endl;
    cout << "Using impl (B) for p=" << p << "  time is " << (EndTime_B - StartTime_B) << endl;

  }


  void program()
  {
    GlobalManager CoCoAFoundations;
    cout << ShortDescription << endl;

    const int p = NextPrime(2048);
    TimeTrialRingElem(p);
    TimeTrialSmallFp(p);
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
