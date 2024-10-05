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

//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-PolyIterator2.C,v 1.16 2024/08/01 21:22:39 abbott Exp $
// $Log: ex-PolyIterator2.C,v $
// Revision 1.16  2024/08/01 21:22:39  abbott
// Summary: Fixed so that it compiles: changed name os some error IDs
//
// Revision 1.15  2024/07/11 15:39:58  bigatti
// Summary: error codes: NotIndet to ReqIndet
//
// Revision 1.14  2022/02/13 09:56:58  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.13  2020/06/17 15:49:19  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.12  2018/09/28 15:54:03  abbott
// Summary: Removed pseudo-ctors NewPolyRing which took just num of indets; now must specify their names
//
// Revision 1.11  2016/11/18 18:04:11  abbott
// Summary: Changed name of loop variable (from i to iter) for SparsePolyIter.
//
// Revision 1.10  2015/06/29 15:41:17  bigatti
// *** empty log message ***
//
// Revision 1.9  2015/06/29 12:45:12  bigatti
// -- code in namespace CoCoA
//
// Revision 1.8  2014/07/08 12:47:26  abbott
// Summary: Removed AsPolyRing, AsSparsePolyRing, AsQuotientRing
// Author: JAA
//
// Revision 1.7  2014/03/05 10:01:57  bigatti
// - -some improvements
//
// Revision 1.6  2012/10/11 16:29:29  abbott
// Replaced RefRingElem by RingElem&; and corrected call to myAddClear in MultiHomog.
//
// Revision 1.5  2012/02/08 17:41:29  bigatti
// -- changed: Z,Q -> ZZ,QQ
//
// Revision 1.4  2011/03/10 17:15:50  abbott
// Replaced size_t by long.
//
// Revision 1.3  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.2  2009/12/11 11:46:32  abbott
// Changed fn  convert  into  IsConvertible.
// Added template procedure  convert.
// New version because change is not backward compatible.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.8  2007/03/08 14:38:07  cocoa
// Added new range function in symbol.H, and tidied many calls to PolyRing
// pseudo ctors (as a consequence).
//
// Revision 1.7  2007/03/03 14:15:45  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.6  2007/02/26 17:41:16  bigatti
// -- changed "homog" --> "MultiHomog"
// -- getting ready for unique ring Z: using NewZmod(N), NewRingQ()
//
// Revision 1.5  2007/02/22 17:31:20  bigatti
// -- added call to CoCoA::homog
//
// Revision 1.4  2007/02/12 15:45:15  bigatti
// -- added strings ShortDescription and LongDescription for indexing
//
// Revision 1.3  2007/02/10 18:44:04  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.2  2006/08/17 10:05:27  cocoa
// -- changed: SmallExponent_t --> long
//
// Revision 1.1.1.1  2006/05/30 11:39:36  cocoa
// Imported files
//
// Revision 1.1.1.1  2005/10/17 10:46:53  cocoa
// Imported files
//
// Revision 1.5  2005/10/14 15:25:07  cocoa
// Major tidying and cleaning to small prime finite fields.
// Several consequential changes.  Improved their documentation.
//
// Added Makefile and script to include/CoCoA/ directory to
// keep library.H up to date.
//
// Revision 1.4  2005/10/12 15:52:09  cocoa
// Completed test-RingFp1 and corrected/cleaned the SmallFp*
// and RingFp* files.
//
// Some minor tidying elsewhere.
//
// Revision 1.3  2005/10/07 17:08:09  cocoa
// Fixed an embarrassing buglet -- I ought to try compiling before
// I check (and there are plenty of other things I ought to do too...)
//
// Revision 1.2  2005/10/06 16:36:42  cocoa
// Added the capability find out build information at run-time.
// The Makefiles should be a little tidier too.
//
// Revision 1.1  2005/08/08 16:36:33  cocoa
// Just checking in before going on holiday.
// Don't really recall what changes have been made.
// Added IsIndet function for RingElem, PPMonoidElem,
// and a member function of OrdvArith.
// Improved the way failed assertions are handled.
//
// Revision 1.2  2005/07/19 15:30:20  cocoa
// A first attempt at iterators over sparse polynomials.
// Main additions are to SparsePolyRing, DistrMPoly*.
// Some consequential changes to PPMonoid*.
//
// Revision 1.1  2005/05/04 16:35:06  cocoa
// -- new examples
//
