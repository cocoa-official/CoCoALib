// Copyright (c) 2005  John Abbott,  Anna Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "Predefined and user-defined orderings and gradings \n"
  "on PPMonoid and PolyRing.                          \n";

const string LongDescription =
  "Each ordering is degree-compatible with grading over Z^GradingDim \n"
  "i.e. the grading is given by the first GradingDim rows            \n"
  "of the ordering matrix.                                           \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  //-- auxiliary ---------------------------------------------------------
  // This is an ad hoc function for converting a basic C structure
  // into a CoCoA matrix.

  // convention: a function containing a "new" should be named "New.."
  matrix NewZMatrixFromC(int cmat[4][4])
  {
    matrix M(NewDenseMat(RingZZ(),4,4));
  
    for (int i=0; i < 4; ++i)
      for (int j=0; j < 4; ++j)
        SetEntry(M, i, j, cmat[i][j]);
    return M;
  }


  //-- TestOrdering ------------------------------------------------------
  // behaviour of different orderings and gradings on PPMonoid and PolyRing

  void TestOrdering(const PPOrdering& ord)
  {
    cout << "  GradingDim = " << GradingDim(ord) << endl;
    cout << "  NumIndets  = " << NumIndets(ord) << endl;

    if (NumIndets(ord)<3) 
      CoCoA_THROW_ERROR("ord has less than 3 indets", "TestOrdering");

    // Indet names are x[0], x[1], ...
    vector<symbol> IndetNames = SymbolRange("x", 0, NumIndets(ord)-1);
    PPMonoid PPM = NewPPMonoidEvOv(IndetNames, ord);

    // For handy access to the indeterminates in PPM
    vector<PPMonoidElem> x;
    for (long i=0; i < NumIndets(PPM); ++i)
      x.push_back(indet(PPM, i));

  
    PPMonoidElem t1 = power(x[0],4) * power(x[1],3) * x[2];
    PPMonoidElem t2 = power(x[0],2) * power(x[1],9);
    cout << "  t1 = " << t1 << endl;
    cout << "  t2 = " << t2 << endl;
 
    cout << "wdeg(t1)   gives  " << wdeg(t1) << endl;
    cout << "wdeg(t2)   gives  " << wdeg(t2) << endl;
    cout << "t1 < t2   gives  " << (t1 < t2) << endl;

    // Now create Zx a multigraded ring of polynomials with coefficients in Z;
    // NewPolyRing_DMPI works better than NewPolyRing for matrix ordering.
    PolyRing Zx = NewPolyRing_DMPI(RingZZ(), IndetNames, ord);

    RingElem f = indet(Zx,0) + indet(Zx,3);
    cout << "  f = " << f << endl;  

    cout << "deg(f) gives " << deg(f) << endl;
    cout << "------------------------------" << endl;
  }


  //-- program --------------------------------------------------------------
  // we run TestOrdering on predefined and user-defined orderings

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    // user-defined ordering and grading
    const int GradingDim = 2;

    // the first 2 rows represent the degree matrix
    int M[4][4] = {{2, 0, 0, 3},
                   {1, 2, 4, 0},
                   {1, 0, 0, 0},
                   {0, 1, 0, 0}};

    PPOrdering MatOrd4 = NewMatrixOrdering(NewZMatrixFromC(M), GradingDim);

    // TestOrdering calls
    cout << "  ordering: lex(4)" << endl;
    TestOrdering(lex(4));
  
    cout << "  ordering: StdDegLex(4)" << endl;
    TestOrdering(StdDegLex(4));
  
    cout << "  ordering: StdDegRevLex(4)" << endl;
    TestOrdering(StdDegRevLex(4));
  
    cout << "  ordering: MatOrd(4)" << endl;
    TestOrdering(MatOrd4);
  }

} // end of namespace CoCoA

//----------------------------------------------------------------------
// We write main() like this so we can handle uncaught CoCoA errors in
// a sensible way (i.e. by announcing them).
int main()
{
  try
  {
    CoCoA::program();
    return 0;
  }
  catch (const CoCoA::ErrorInfo& err)
  {
    cerr << "***ERROR***  UNCAUGHT CoCoA error" << endl;
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
