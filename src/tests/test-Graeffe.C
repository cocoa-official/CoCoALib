
#include "CoCoA/assert.H"
#include "CoCoA/BuildInfo.H"
#include "CoCoA/error.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/SparsePolyOps-RingElem.H"
#include "CoCoA/SparsePolyOps-Graeffe.H"

#include <iostream>
using std::cerr;
using std::endl;

namespace CoCoA
{

  void TestGraeffeN()
  {
    SparsePolyRing P = NewPolyRing(RingZZ(), symbols("x")); // ZZ[x];

    RingElem f = RingElem(P, "3*x^2+4*x+2");
    GraeffeSeq graeffeSeq(f, 6);

    RingElem g = f;
    ++graeffeSeq;
    CoCoA_ASSERT_ALWAYS(GraeffeN(f, 1) == f);
    CoCoA_ASSERT_ALWAYS(*graeffeSeq == f);

    g = graeffe(f);
    ++graeffeSeq;
    CoCoA_ASSERT_ALWAYS(GraeffeN(f, 2) == g);
    CoCoA_ASSERT_ALWAYS(*graeffeSeq == g);

    g = graeffe3(f);
    ++graeffeSeq;
    CoCoA_ASSERT_ALWAYS(GraeffeN(f, 3) == g);
    CoCoA_ASSERT_ALWAYS(*graeffeSeq == g);

    g = graeffe(graeffe3(f));
    ++graeffeSeq;
    ++graeffeSeq;
    ++graeffeSeq;
    CoCoA_ASSERT_ALWAYS(GraeffeN(f, 6) == g);
    CoCoA_ASSERT_ALWAYS(*graeffeSeq == g);

    SparsePolyRing Pxy = NewPolyRing(RingZZ(), symbols("x,y")); // ZZ[x,y];
    RingElem m = RingElem(Pxy, "y^2+x+1");
    try
    {
      GraeffeN(m, 5);
      CoCoA_ASSERT_ALWAYS(!"Non-univariate polynomials should throw in GraeffeN!");
    }
    catch (const CoCoA::ErrorInfo &err)
    {
      if (err != ERR::ReqUnivariate) throw;
    }
  }

  void program()
  {
    GlobalManager CoCoAFoundations;

    TestGraeffeN();
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
  catch (const CoCoA::ErrorInfo &err)
  {
    cerr << "***ERROR***  UNCAUGHT CoCoA Error";
    ANNOUNCE(cerr, err);
  }
  catch (const std::exception &exc)
  {
    cerr << "***ERROR***  UNCAUGHT std::exception: " << exc.what() << endl;
  }
  catch (...)
  {
    cerr << "***ERROR***  UNCAUGHT UNKNOWN EXCEPTION" << endl;
  }

  CoCoA::BuildInfo::PrintAll(cerr);
  return 1;
}
