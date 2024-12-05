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
