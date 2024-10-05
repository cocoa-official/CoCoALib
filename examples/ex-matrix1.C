// Copyright (c) 2005  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "Example program illustrating the creation of matrices, and some \n"
  "operations on them. \n";

const string LongDescription =
  "Example program illustrating the creation of matrices, and some \n"
  "basic operations on them. \n";

//----------------------------------------------------------------------


namespace CoCoA
{


// convention: a function containing a "new" should be named "New.."
  matrix NewMatrixFromC(ring R, int cmat[4][4])
  {
    matrix M(NewDenseMat(R,4,4));
  
    for (int i=0; i < 4; ++i)
      for (int j=0; j < 4; ++j)
        SetEntry(M, i, j, cmat[i][j]);
    return M;
  }


  void ExBasicOps(matrix M)
  {
    cout << "M = " << M << endl;
    cout << "rk(M) = " << rk(M) << endl;
    cout << "det(M) = " << det(M) << endl;

    matrix InvM = inverse(M);  
    cout << "InvM = " << InvM << endl;
    InvM = InvM * M;  
    cout << "InvM * M = " << InvM << endl;

    matrix AdjM = adj(M);  
    cout << "AdjM = " << AdjM << endl;  
    AdjM = AdjM * M;
    cout << "AdjM * M = " << AdjM << endl;
  }


  void ExMatrixView(matrix M)
  {
    cout << "-- some examples with ConstMatrixView/MatrixView --"<< endl;
    cout << "see the documentation of matrix for the full list"<< endl;
  
    // It is easy to compute the determinant of an identity matrix, even if large.
    ConstMatrixView IdentBig = IdentityMat(RingOf(M), 1000000);
    cout << "det(IdentBig) = " << det(IdentBig) << endl << endl;

    ConstMatrixView Ident4 = IdentityMat(RingOf(M), 4);
    ConstMatrixView  A = ConcatHor(M, Ident4);
    cout << "ConcatHor(M, Ident4) = " << A << endl;
  
    ConstMatrixView Z2x8 = ZeroMat(RingOf(M), 2, 8);
    ConstMatrixView  B = ConcatVer(A, Z2x8);
    cout << "ConcatVer(A, Z2x8) = " << B << endl;

    // with MatrixView you can change the original matrix M
    //   using a "different point of view"
    // (remember: in C/C++ indices start from 0)
    cout << M << endl;
    // transpose
    MatrixView TrM = transpose(M);  SetEntry(TrM, 1,0, 11);
    cout << "TrM = transpose(M);  SetEntry(TrM, 1,0, 11);" << endl
         << "M = " << M << endl;

    // if you do not want to modify M make a copy into a proper matrix
    //  matrix TrM1 = NewDenseMat(transpose(M));

    vector<long> L12;
    L12.push_back(1); L12.push_back(2);
    // submat
    MatrixView SM = submat(M, L12, L12);  SetEntry(SM, 1,0, 22);
    cout << "SM = submat(M, L12, L12);  SetEntry(SM, 1,0, 22);" << endl
         << "M = " << M << endl;

    // BlockMat2x2
    MatrixView BM = BlockMat2x2(M, M, M, M);  SetEntry(BM, 7,0, 33);
    cout << "BM = BlockMat2x2(M, M, M, M);  SetEntry(BM, 7,0, 33);" << endl
         << "BM = " << BM << endl;
  }


  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    // the first 2 rows represent the degree matrix
    int C_matrix[4][4] = {{1, 2, 3, 4},
                          {2, 3, 4, 0},
                          {4, 0, 1, 2},
                          {0, 1, 2, 3}};

    // Create two CoCoALib matrices containing the values of M.
    // M_Z5 contains the image of M in Z/(11), and M_Q the image in Q.
    matrix M_Z11(NewMatrixFromC(NewZZmod(11), C_matrix));
    matrix M_Q(NewMatrixFromC(RingQQ(), C_matrix));

    ExBasicOps(M_Z11);
    ExBasicOps(M_Q);

    ExMatrixView(M_Q);
  }

} // end of namespace CoCoA


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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-matrix1.C,v 1.17 2022/02/13 09:57:01 abbott Exp $
// $Log: ex-matrix1.C,v $
// Revision 1.17  2022/02/13 09:57:01  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.16  2021/10/14 12:09:48  abbott
// Summary: Fixed typo in printed string
//
// Revision 1.15  2020/06/22 15:39:21  abbott
// Summary: Replaced two calls to mul(,,) with equivalent code
//
// Revision 1.14  2015/06/26 15:34:48  abbott
// Summary: Moved code into namespace CoCoA (see redmine 739)
// Author: JAA
//
// Revision 1.13  2015/04/13 15:28:59  abbott
// Summary: Changed "rank" --> "rk"; adjoint" --> "adj"
// Author: JAA
//
// Revision 1.12  2014/07/30 13:59:08  abbott
// Summary: Changed BaseRing into RingOf
// Author: JAA
//
// Revision 1.11  2013/05/31 15:05:22  bigatti
// changed BlockMat into BlockMat2x2
//
// Revision 1.10  2012/07/04 12:24:49  abbott
// Several minor improvements: better var names, clearer examples, better layout of printout.
//
// Revision 1.9  2012/06/19 14:36:50  bigatti
// -- removed "using CoCoA::rank"
//
// Revision 1.8  2012/02/08 17:52:17  bigatti
// -- changed: Z,Q -> ZZ,QQ
//
// Revision 1.7  2011/09/07 09:56:41  abbott
// Some minor improvements.
//
// Revision 1.6  2011/03/08 18:01:53  bigatti
// -- changed size_t into long
//
// Revision 1.5  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.4  2009/02/13 14:55:26  bigatti
// -- changed a lot: rearranged basic operations and added (Const)MatrixView
//
// Revision 1.3  2008/04/18 15:35:57  abbott
// (long overdue) Major revision to matrices
//
// Revision 1.2  2008/04/08 15:26:42  abbott
// Major revision to matrix implementation: added matrix views.
// Lots of changes.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.6  2007/03/03 14:15:45  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.5  2007/02/26 17:40:34  bigatti
// -- getting ready for unique ring Z: using NewZmod(N), NewRingQ()
//
// Revision 1.4  2007/02/22 17:25:40  bigatti
// -- added printing of separators: "-------"
//
// Revision 1.3  2007/02/12 16:27:43  bigatti
// -- added strings ShortDescription and LongDescription for indexing
//
// Revision 1.2  2007/02/10 18:44:03  cocoa
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
// Revision 1.2  2005/10/14 15:25:07  cocoa
// Major tidying and cleaning to small prime finite fields.
// Several consequential changes.  Improved their documentation.
//
// Added Makefile and script to include/CoCoA/ directory to
// keep library.H up to date.
//
// Revision 1.1.1.1  2005/05/03 15:47:30  cocoa
// Imported files
//
// Revision 1.6  2005/04/27 16:14:56  cocoa
// Cleaned up example programs -- added "free use" permit.
// Changed a couple of ErrorInfo object names, and added
// ERR::NotTrueGCDDomain.
//
// Revision 1.5  2005/04/21 15:12:19  cocoa
// Revised NewPolyRing as Dag Arneson suggested (perhaps just an interim
// measure).
// Brought example programs up to date (new name for CoCoA error
// information objects).
//
// Revision 1.4  2005/04/19 14:06:04  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.3  2005/03/30 15:30:43  cocoa
// -- added: computation of rank and determinant
// -- changed: declaration of I4000 as ConstMatrix
//
// Revision 1.2  2005/03/02 18:46:41  cocoa
// Added new types ConstRefMatrix, and RefMatrix following along
// the lines of ConstRefRingElem and RefRingElem.  The semantics
// should be a bit clearer now.
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.1  2004/11/29 16:50:28  cocoa
// -- first import
//
