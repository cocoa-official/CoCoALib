// Copyright (c) 2005-2012  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "Program showing how to homogenize a sparse polynomial using iterators.\n";

const string LongDescription =
  "This is just an example!  If you want to homogenize polynomials  \n"
  "you should use the library function \"homog\".                   \n";
// ----------------------------------------------------------------------

namespace CoCoA
{

  void HomogCheck(const SparsePolyRing& P, const vector<long>& h)
  {
    // Verify that wdeg(indet(h[i])) is (0,..,0,1,0,..,0) with 1 in i-th position
    const long D = GradingDim(P);
    CoCoA_ASSERT(len(h) == D);
    for (long i=0; i < D; ++i)
    {
      CoCoA_ASSERT(0 <= h[i] && h[i] < NumIndets(P));
      const degree d = wdeg(indet(P,h[i]));
      for (long j=0; j < D; ++j)
        if (i == j)
          CoCoA_ASSERT(IsOne(d[j]));
        else
          CoCoA_ASSERT(IsZero(d[j]));
    }
  }


  // This procedure assumes that args have been sanity checked.
  void MultiHomog(RingElem& hf, ConstRefRingElem f, const vector<long>& h)
  {
    CoCoA_ASSERT(owner(hf) == owner(f)); // explain CoCoA_ASSERT
    const SparsePolyRing P = owner(f);
    const long D = GradingDim(P);

    if (IsZero(f)) { hf = f; return; }// trivial case

    // Compute in H the "top" of the degrees of the PPs in f
    //??? if (D == 1) ... SPECIAL CASE
    degree H(D);
    for (SparsePolyIter iter=BeginIter(f); !IsEnded(iter); ++iter)
    {
      H = top(H, wdeg(PP(iter)));
    }

    // Now homogenize f.  Accumulate result into a geobucket for speed.
    geobucket ans(P);
    vector<long> expv(NumIndets(P));
    for (SparsePolyIter iter=BeginIter(f); !IsEnded(iter); ++iter)
    {
      const degree diff = H - wdeg(PP(iter));
      for (long j=0; j < D; ++j)
        if (!IsConvertible(expv[h[j]], diff[j]))
          CoCoA_THROW_ERROR("Exponent too big", "MultiHomog(ans,f,h)");
      RingElem HomogTerm = monomial(P, coeff(iter), PP(iter)*PPMonoidElem(PPM(P),expv));
      ans.myAddClear(HomogTerm, 1);
    }
    hf = 0; // NB cannot do this earlier in case hf aliases f.
    AddClear(hf, ans);
  }


  RingElem MultiHomog(ConstRefRingElem f, const vector<long>& h)
  {
    CoCoA_ASSERT(IsSparsePolyRing(owner(f)));
    const SparsePolyRing P = owner(f);
    HomogCheck(P, h);
    RingElem ans(P);
    MultiHomog(ans, f, h);
    return ans;
  }


  RingElem MultiHomog(ConstRefRingElem f, const vector<RingElem>& h)
  {
    if (!IsSparsePolyRing(owner(f)))
      CoCoA_THROW_ERROR(ERR::ReqSparsePolyRing, "homog(f,h)");
    const SparsePolyRing P = owner(f);

    const long k = GradingDim(P);
    vector<long> indices(k);
    if (len(h) != k)
      CoCoA_THROW_ERROR(ERR::BadArraySize, "homog(f,h)");
    for (long i=0; i < k; ++i)
    {
      if (owner(h[i]) != P)
        CoCoA_THROW_ERROR(ERR::MixedRings, "homog(f,h)");
      long index;
      if (!IsIndet(index, h[i]))
        CoCoA_THROW_ERROR(ERR::ReqIndet, "homog(f,h)");

      indices[k] = index;
    }
    HomogCheck(P, indices);
    RingElem ans(P);
    MultiHomog(ans, f, indices);
    return ans;
  }


  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    ring Fp = NewZZmod(32003);                 // coefficient ring
    SparsePolyRing Fpx = NewPolyRing(Fp, SymbolRange("x",0,3));
    const vector<RingElem>& x = indets(Fpx);
    vector<long> h;
    h.push_back(0);
    RingElem f = x[0] + 3*power(x[1],2) + 5*power(x[2],4) + 7*power(x[3],8);
    cout << "Original f  = " << f << endl;
    cout << "MultiHomog(f, h) =  " << MultiHomog(f, h) << endl;
    cout << "homog(f, x[0]) = " << homog(f, x[0]) << endl;
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
