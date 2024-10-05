// Copyright (c) 2011  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "Example program illustrating the creation of matrix views, and some  \n"
  "effects of \"aliasing\".  Examples of matrix views include transposes\n"
  "submatrices, block matrices, and concatenated matrices.              \n";

const string LongDescription =
  "Example program illustrating the creation of matrix views, and some    \n"
  "effects of \"aliasing\", i.e. where there is more than one way to refer\n"
  "to a single entry of the matrix.  We gives examples of creating various\n"
  "views: ZeroMat, IdentityMat, transpose, submat, BlockMat2x2, ConcatHor \n"
  "and ConcatVer.";


//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    // We could use any ring for these examples, but have chosen ZZ for simplicity.
    ring R = RingZZ();

    // Here is a normal 3x3 matrix (with entries 1,2,3, 4,5,6, 7,8,9)
    // we shall use it later on.
    matrix M = NewDenseMat(R,3,3);
    SetEntry(M,0,0, 1);  SetEntry(M,0,1, 2);  SetEntry(M,0,2, 3);
    SetEntry(M,1,0, 4);  SetEntry(M,1,1, 5);  SetEntry(M,1,2, 6);
    SetEntry(M,2,0, 7);  SetEntry(M,2,1, 8);  SetEntry(M,2,2, 9);


    // Zero matrix -- useful mainly in calls to BlockMat2x2 (see below)
    ConstMatrix Zero2x3 = ZeroMat(R, 2,3);
    cout << "*** ZERO MATRIX ***\n"
         << "Here is a 2x3 zero matrix: " << Zero2x3 << endl;
    cout << "It is constant; none of its entries may be assigned to." << endl
         << endl;


    // Identity matrix
    ConstMatrix Id3x3 = IdentityMat(R, 3);
    cout << "*** IDENTITY MATRIX ***\n"
         << "Here is a 3x3 identity matrix: " << Id3x3 << endl
         << "It is constant; none of its entries may be assigned to." << endl
         << endl;


    // Transpose -- entries are shared with original matrix.
    cout << "*** TRANSPOSE MATRIX ***\n"
         << "Our starting matrix is  M = " << M << endl;
    MatrixView TrM = transpose(M);
    cout << "transpose(M) = " << TrM << endl;

    cout << "If we change an entry in M, then TrM changes too:" << endl;
    SetEntry(M,0,1, 99);
    cout << "After setting M(0,1) = 99, the matrix M becomes: " << M << endl;
    cout << "and its transpose changes correspondingly: TrM = " << TrM << endl;

    cout << "Similarly, changing an entry of TrM also changes M..." << endl;
    SetEntry(TrM,1,0, 2);
    cout << "After setting TrM(0,1) = 2, the matrix TrM becomes: " << TrM << endl
         << "and the original matrix changes correspondingly: M = " << M << endl
         << endl;

    cout << "If we don't wanted shared entries, we must make an explicit copy," << endl
         << "e.g. by calling NewDenseMat(transpose(M))" << endl;
    matrix TrM_copy = NewDenseMat(transpose(M));
    SetEntry(TrM_copy,2,2, 999);
    cout << "Changing an entry in the copy does not affect the original matrix:" << endl
         << "TrM_copy has been changed to " << TrM_copy << endl
         << "but this did not change the original M = " << M << endl
         << endl;


    // Submatrices -- entries are shared with parent matrix.
    cout << "*** SUBMATRIX ***\n"
         << "We shall create the submatrix view of M which comprises just its four corners." << endl;
    vector<long> rows; rows.push_back(0); rows.push_back(2); // rows = [0,2]
    vector<long> cols; cols.push_back(0); cols.push_back(2); // cols = [0,2]
    MatrixView subm = submat(M, rows, cols);

    cout << "The submatrix is subm = " << subm << endl
         << "as for the transpose, its entries are shared with the original matrix." << endl
         << endl;


    // Block matrices -- sticking 4 matrices together to make a bigger one.
    cout << "*** BLOCK MATRIX ***\n";
    MatrixView BlockM = BlockMat2x2(M, zeroes, zeroes, TrM);
    cout << "We can stick 4 matrices together to make one larger one.  For\n"
         << "example BlockMat2x2(M,Z,Z,TrM) where Z is a 3x3 zero matrix gives:\n"
         << BlockM << endl;

    cout << "Like submat, BlockMat2x2 does not make copies of the original matrices\n"
         << "it simply refers to them.  In our specific case we must be extra\n"
         << "careful, since there is aliasing between M and TrM, so there are\n"
         << "'hidden' connections between the entries of the block matrix.\n"
         << "Changing the (0,1) entry will also change the (4,3) entry, like this:\n";
    SetEntry(BlockM,0,1, -99);
    cout << BlockM << endl;
    cout << "Naturally, also the original matrix M has changed; we now have\n"
         << "M = " << M << endl
         << endl
         << "We can avoid this phenomenon by making an explicit copy, as we did\n"
         << "for the transpose (above)." << endl
         << endl;

    // Another BlockMat2x2
    long n = 4;
    RingElem I = one(RingZZ());
    RingElem Z = zero(RingZZ());
    cout << "Another BlockMat2x2 " <<
      BlockMat2x2(
                  RowMat(vector<RingElem>(n,I)), RowMat(vector<RingElem>(1,I)),
                  StdDegLexMat(n),               ColMat(vector<RingElem>(n,Z)));
    

    // ConcatHor & ConcatVer join two matrices horizontally/vertically
    cout << "*** ConcatHor & ConcatVer ***\n";
    ConstMatrix Zero3x3 = ZeroMat(R, 3,3);
    cout << "*** CONCATENATED MATRICES ***\n"
         << "We can join two matrices horizontally or vertically:\n"
         << "ConcatHor(M,Z) = " << ConcatHor(M,Zero3x3) << endl
         << endl
         << "ConcatVer(M,Z) = " << ConcatVer(M,Zero3x3) << endl
         << "Like BlockMat2x2, these commands don't make copies of the original matrices." << endl;
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-matrix2.C,v 1.10 2022/02/13 09:57:01 abbott Exp $
// $Log: ex-matrix2.C,v $
// Revision 1.10  2022/02/13 09:57:01  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.9  2020/03/06 17:39:01  bigatti
// -- added another BlockMat2x2
//
// Revision 1.8  2015/11/30 21:53:55  abbott
// Summary: Major update to matrices for orderings (not yet complete, some tests fail)
//
// Revision 1.7  2015/06/26 15:34:48  abbott
// Summary: Moved code into namespace CoCoA (see redmine 739)
// Author: JAA
//
// Revision 1.6  2014/07/16 15:50:08  abbott
// Summary: Corrected type in printed mesg
// Author: JAA
//
// Revision 1.5  2013/05/31 15:05:22  bigatti
// changed BlockMat into BlockMat2x2
//
// Revision 1.4  2012/06/19 15:44:58  abbott
// Improved long and short descriptions.
//
// Revision 1.3  2012/02/08 17:52:17  bigatti
// -- changed: Z,Q -> ZZ,QQ
//
// Revision 1.2  2011/12/23 15:27:35  bigatti
// -- cleaning
//
// Revision 1.1  2011/12/23 15:06:53  abbott
// New example for matrix views
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
