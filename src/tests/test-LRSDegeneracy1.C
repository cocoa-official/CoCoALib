
#include "CoCoA/assert.H"
#include "CoCoA/BuildInfo.H"
#include "CoCoA/error.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/SparsePolyOps-LRSDegeneracy.H"

#include <iostream>
using std::cerr;
using std::endl;

namespace CoCoA
{

  void TestIsLRSDegenerate()
  {
    SparsePolyRing P = NewPolyRing(RingZZ(), symbols("x")); // ZZ[x];

    RingElem f = RingElem(P, "x^4+3*x^3+2*x^2-3*x+1");
    CoCoA_ASSERT_ALWAYS(LRSDegeneracyOrder(f) == 6);

    RingElem g = RingElem(P, "x^4+2*x^3+4*x^2+4*x+4");
    CoCoA_ASSERT_ALWAYS(LRSDegeneracyOrder(g) == 2);

    RingElem h = RingElem(P, "x^4+x^3-3*x^2+x+1");
    CoCoA_ASSERT_ALWAYS(LRSDegeneracyOrder(h) == 0);

    RingElem j = RingElem(P, "(x-1)*(x+2)");
    CoCoA_ASSERT_ALWAYS(LRSDegeneracyOrder(j) == 0);

    try
    {
      SparsePolyRing Pxy = NewPolyRing(RingZZ(), symbols("x,y")); // ZZ[x,y];
      const RingElem m = RingElem(Pxy, "y^2+x+1");
      LRSDegeneracyOrder(m);
      CoCoA_ASSERT_ALWAYS(!"Non-univariate polynomials should throw in IsLRSDegenerate!");
    }
    catch (const CoCoA::ErrorInfo &err)
    {
      if (err != ERR::ReqUnivariate) throw;
    }
  }

  void program()
  {
    GlobalManager CoCoAFoundations;

    TestIsLRSDegenerate();
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
