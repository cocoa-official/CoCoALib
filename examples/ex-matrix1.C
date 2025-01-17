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
