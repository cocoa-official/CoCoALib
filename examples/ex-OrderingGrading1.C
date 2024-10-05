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

//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-OrderingGrading1.C,v 1.11 2022/02/13 09:56:57 abbott Exp $
// $Log: ex-OrderingGrading1.C,v $
// Revision 1.11  2022/02/13 09:56:57  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.10  2020/06/17 15:49:19  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.9  2017/11/10 16:02:27  abbott
// Summary: Removed NewLexOrdering, NewStdDegLexOrdering, NewStdDegRevLexOrdering; consequential changes
//
// Revision 1.8  2015/12/01 12:26:46  abbott
// Summary: Corrected call to NewMatrixOrdering
//
// Revision 1.7  2015/06/29 15:32:38  bigatti
// *** empty log message ***
//
// Revision 1.6  2015/06/29 12:38:13  bigatti
// -- code in namespace cocoa
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
// Revision 1.2  2008/04/18 15:35:57  abbott
// (long overdue) Major revision to matrices
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.8  2007/03/08 17:43:11  cocoa
// Swapped order of args to the NewPPMonoid pseudo ctors.
//
// Revision 1.7  2007/03/08 16:55:06  cocoa
// Changed name of "range" function to "SymbolRange".
//
// Revision 1.6  2007/03/03 14:15:45  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.5  2007/02/28 15:19:38  bigatti
// -- minor: removed old comment about adjoint matrix
//
// Revision 1.4  2007/02/26 15:49:08  bigatti
// -- getting ready for unique ring Z: using NewZmod(N), NewRingQ()
//
// Revision 1.3  2007/02/12 15:31:57  bigatti
// -- added strings ShortDescription and LongDescription for indexing
//
// Revision 1.2  2007/02/10 18:44:04  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.1.1.1  2006/05/30 11:39:36  cocoa
// Imported files
//
// Revision 1.1.1.1  2005/10/17 10:46:53  cocoa
// Imported files
//
// Revision 1.2  2005/07/19 15:30:20  cocoa
// A first attempt at iterators over sparse polynomials.
// Main additions are to SparsePolyRing, DistrMPoly*.
// Some consequential changes to PPMonoid*.
//
// Revision 1.1.1.1  2005/05/03 15:47:30  cocoa
// Imported files
//
// Revision 1.4  2005/04/27 16:14:56  cocoa
// Cleaned up example programs -- added "free use" permit.
// Changed a couple of ErrorInfo object names, and added
// ERR::NotTrueGCDDomain.
//
// Revision 1.3  2005/04/21 15:12:19  cocoa
// Revised NewPolyRing as Dag Arneson suggested (perhaps just an interim
// measure).
// Brought example programs up to date (new name for CoCoA error
// information objects).
//
// Revision 1.2  2005/04/19 14:06:05  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.3  2004/11/29 16:39:21  cocoa
// -- PPOrdering no longer needs to know the adjoint of the order matrix
//
