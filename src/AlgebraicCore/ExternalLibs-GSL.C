//   Copyright (c)  2011  Bruno Simoes
//   2022-06-07 partially cleaned by JAA: THIS CODE IS DEFINITELY NOT EXCEPTION SAFE!!!!

//   This file is part of the source of CoCoALib, the CoCoA Library.
//
//   CoCoALib is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   CoCoALib is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with CoCoALib.  If not, see <http://www.gnu.org/licenses/>.

#include "CoCoA/ExternalLibs-GSL.H"

#ifdef CoCoA_WITH_GSL

#include "CoCoA/BigRat.H"
#include "CoCoA/DenseMatrix.H"
#include "CoCoA/MatrixView.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/convert.H"
#include "CoCoA/error.H"
#include "CoCoA/matrix.H"

#include <iostream>
using std::ostream;


namespace CoCoA
{
  // GSL matrices to CoCoA matrices
  matrix NewMatrix(const ring& R, gsl_matrix* m)
  {
    const size_t M = m->size1;
    const size_t N = m->size2;

    ring myRing = R;
    if (IsZZ(R))
      myRing = RingQQ();
  
    matrix mat(NewDenseMat(myRing,M,N));
  
    for (size_t i=0; i < M; ++i)
    {
      for (size_t j=0; j < N; ++j)
      {
        SetEntry(mat, i, j, ConvertTo<BigRat>(gsl_matrix_get(m,i,j)));
      }
    }

    return mat;
  }


  matrix NewMatrix(const ring& R, gsl_matrix* m, bool transpose)
  {
    /*const*/ size_t M = m->size1;
    /*const*/ size_t N = m->size2;
    if (transpose)
    {
      std::swap(M,N);
    }

    ring myRing = R;
    if (IsZZ(R))
      myRing = RingQQ();
 
    matrix mat(NewDenseMat(myRing, M, N));

    if (transpose)
    {
      for (size_t i=0; i < M; ++i)
      {
        for (size_t j=0; j < N; ++j)
        {
          SetEntry(mat, i, j, ConvertTo<BigRat>(gsl_matrix_get(m,j,i)));
        }
      }
    }
    else {
      for (size_t i=0; i < M; ++i)
      {
        for (size_t j=0; j < N; ++j)
        {
          SetEntry(mat, i, j, ConvertTo<BigRat>(gsl_matrix_get(m,i,j)));
        }
      }
    }
    return mat;
  }


  // GSL vectors to CoCoA matrices
  matrix NewMatrix(const ring& R, gsl_vector* m)
  {
    const size_t M = m->size;

    ring myRing = R;
    if (IsZZ(R))
      myRing = RingQQ();
    
    matrix mat(NewDenseMat(myRing,1,M));

    for (size_t i = 0; i<M; ++i)
    {
      SetEntry(mat, 0, i, ConvertTo<BigRat>(gsl_vector_get(m,i)));
    }
    return mat;
  }


  // CoCoA matrices to GSL matrices
  gsl_matrix* NewMatrixGSL(ConstMatrixView N)
  {
    const size_t cols = NumCols(N);
    const size_t rows = NumRows(N);
    gsl_matrix *m = gsl_matrix_alloc(rows,cols);

    BigRat q;
    for (size_t i=0; i < rows; ++i)
    {
      for (size_t j=0; j < cols; ++j)
      {
        if (!IsRational(q, N(i,j)))
          CoCoA_THROW_ERROR2(ERR::BadArg, "Cannot convert to double");
        gsl_matrix_set(m, i,j, mpq_get_d(mpqref(q) ));
      }
    }
    return m;
  }


  void UpdateMatrix(matrix& A, gsl_matrix* m)
  {
    // todo:add exception sizes
    ring R = RingOf(A);
    if (IsZZ(R))
    {
      R = RingQQ();
      const RingHom& embed = ZZEmbeddingHom(RingQQ());
      A = embed(A);
    }
  
    const size_t M = m->size1;
    const size_t N = m->size2;  
    for (size_t i=0; i < M; ++i)
      for (size_t j=0; j < N; ++j)
      {
        SetEntry(A, i, j, ConvertTo<BigRat>(gsl_matrix_get(m,i,j)));
      }
  }


  // CoCoA vectors to GSL vectors
  gsl_vector* NewVectorGSL(std::vector<RingElem>& v)
  {
    const size_t cols = v.size();
    gsl_vector *m = gsl_vector_alloc(cols);

    BigRat q;
    for (size_t j=0; j < cols; ++j)
    {
      if (!IsRational(q, v[j]))
        CoCoA_THROW_ERROR2(ERR::BadArg, "Cannot convert to double");
      gsl_vector_set(m, j, mpq_get_d(mpqref(q) ));
    }
    return m;
  }


  std::vector<RingElem> NewVector(const ring& R, gsl_vector* m)
  {
    ring myRing = R;
    if (IsZZ(R))
      myRing = RingQQ();

    const size_t M = m->size;
    std::vector<RingElem> x(M, zero(myRing));
    for (size_t i=0; i < M; ++i)
    {
      x[i] = ConvertTo<BigRat>(gsl_vector_get(m,i));
    }
    return x;
  }


  // Print GSL vector
  void PrintVector(gsl_vector *m)
  {
    const size_t M = m->size;
    for (size_t i = 0; i<M; ++i)
      std::cout << gsl_vector_get(m,i) << " ";
    std::cout << std::endl;
  }


  // Print GSL matrix
  void PrintMatrix(gsl_matrix *m)
  {
    const size_t M = m->size1;
    const size_t N = m->size2;
    for (size_t i = 0; i<M; ++i)
    {
      for (size_t j = 0; j<N; ++j)
        std::cout << gsl_matrix_get(m,i,j) << " ";
      std::cout << std::endl;
    }
  }


  // SVD 
  std::vector<RingElem> GslSingularValues(ConstMatrixView A)
  {
    std::vector<RingElem> sv;

    const size_t cols = NumCols(A);
    const size_t rows = NumRows(A);
    if (cols == 0 || rows == 0) return sv;
    if (rows < cols) A = transpose(A);
  
    gsl_matrix *m = NewMatrixGSL(A);
    const size_t N = m->size2;
    gsl_matrix *q = gsl_matrix_alloc(N,N);
    gsl_vector *d = gsl_vector_alloc(N);
    gsl_vector *w = gsl_vector_alloc(N);
  
    int status = gsl_linalg_SV_decomp(m, q, d, w);   
    if (status != GSL_SUCCESS)
      CoCoA_THROW_ERROR2(ERR::ExternalLib, gsl_strerror(status)); // MEMORY LEAK!  
   
    const ring R = RingOf(A);
    sv = NewVector(R, d);

    gsl_matrix_free(m);
    gsl_matrix_free(q);
    gsl_vector_free(d);
    gsl_vector_free(w); 
  
    return sv;
  }


  /* Factorise a general M x N matrix A into,
   *
   *   A = U D V^T
   *
   * where U is a column-orthogonal M x N matrix (U^T U = I), 
   * D is a diagonal N x N matrix, 
   * and V is an N x N orthogonal matrix (V^T V = V V^T = I)
   *
   * U is stored in the original matrix A, which has the same size
   *
   * V is stored as a separate matrix (not V^T). You must take the
   * transpose to form the product above.
   *
   * The diagonal matrix D is stored in the vector S,  D_ii = S_i
   */

  std::vector<matrix> GslSVD(ConstMatrixView A)
  {

    const size_t cols = NumCols(A);
    const size_t rows = NumRows(A);
  
    std::vector<matrix> sv;
    if (cols == 0 || rows == 0) return sv;
    if (rows < cols) A = transpose(A);
  
    gsl_matrix *m = NewMatrixGSL(A);
 
    const size_t N = m->size2;
    gsl_matrix *q = gsl_matrix_alloc(N,N);
    gsl_vector *d = gsl_vector_alloc(N);
    gsl_vector *w = gsl_vector_alloc(N);

    int status = gsl_linalg_SV_decomp(m, q, d, w);   
    if (status != GSL_SUCCESS)
      CoCoA_THROW_ERROR2(ERR::ExternalLib, gsl_strerror(status)); // MEMORY LEAK!  
    
    const ring R = RingOf(A);

    // A^T = UDV ==> A = V^TDU^T
    matrix U = NewMatrix(R, m, rows < cols);
    matrix E = NewMatrix(R, d);
    matrix V = NewMatrix(R, q, rows >= cols);
 
    gsl_vector_free(w); 
    gsl_vector_free(d);
    gsl_matrix_free(q);
    gsl_matrix_free(m);
  
    std::vector<matrix> result;
    result.push_back(U);
    result.push_back(E);
    result.push_back(V);
    return result;
 
  }


  /*  Solves the system A x = b using the SVD factorization
   *
   *  A = U S V^T
   *
   *  to obtain x. For M x N systems it finds the solution in the least
   *  squares sense.  
   */
  std::vector<RingElem> GslSolveSVD(ConstMatrixView A, std::vector<RingElem>& b)
  {
  
    const size_t cols = NumCols(A);
    const size_t rows = NumRows(A);
    std::vector<RingElem> x;
    if (cols == 0 || rows == 0) return x;
    if (cols > rows) transpose(A);
  
    gsl_matrix* m = NewMatrixGSL(A);
    const size_t N = m->size2;
    gsl_matrix* q  = gsl_matrix_alloc(N,N);
    gsl_vector* d  = gsl_vector_alloc(N);
    gsl_vector* w  = gsl_vector_alloc(N);
  
    int status = gsl_linalg_SV_decomp(m, q, d, w);   
    if (status != GSL_SUCCESS)
      CoCoA_THROW_ERROR2(ERR::ExternalLib, gsl_strerror(status)); // MEMORY LEAK!
    
    const ring R = RingOf(A);
    matrix U(NewMatrix(R, m));
    matrix E(NewMatrix(R, d));
    matrix V(NewMatrix(R, q));

    gsl_vector* vb = NewVectorGSL(b);
    gsl_vector* vx = gsl_vector_alloc(b.size());

    status = gsl_linalg_SV_solve (m, q, d, vb, vx);
    if (status != GSL_SUCCESS)
      CoCoA_THROW_ERROR2(ERR::ExternalLib, gsl_strerror(status)); // MEMORY LEAK!
  
    x = NewVector(R, vx);
  
    gsl_vector_free(vx); 
    gsl_vector_free(vb); 
    gsl_vector_free(w); 
    gsl_vector_free(d);
    gsl_matrix_free(q);
    gsl_matrix_free(m);
  
    return x;
  }


  // LU decomposition
  matrix GslLU(ConstMatrixView A)
  {
    const size_t rows = NumRows(A);
    const size_t cols = NumCols(A);
    if (rows != cols)  CoCoA_THROW_ERROR1(ERR::ReqSquareMatrix);
    if (0 == cols)  CoCoA_THROW_ERROR1(ERR::BadMatrixSize);
  
    gsl_matrix* m = NewMatrixGSL(A);

    int s;
    gsl_permutation* p = gsl_permutation_alloc(cols);

    int status = gsl_linalg_LU_decomp (m, p, &s);
    if (status != GSL_SUCCESS)
      CoCoA_THROW_ERROR2(ERR::ExternalLib, gsl_strerror(status)); // MEMORY LEAK!
  
    const ring R = RingOf(A);
    matrix lu = NewMatrix(R, m);

    gsl_permutation_free(p);
    gsl_matrix_free(m);
 
    return lu;
  }


  // LU solver
  std::vector<RingElem> GslSolverLU(ConstMatrixView A, std::vector<RingElem>& b)
  {
    const size_t rows = NumRows(A);
    const size_t cols = NumCols(A);
    if (rows != cols)  CoCoA_THROW_ERROR1(ERR::ReqSquareMatrix);
    if (0 == cols || cols != b.size())  CoCoA_THROW_ERROR1(ERR::BadMatrixSize);
    
    gsl_matrix* m = NewMatrixGSL(A);
    gsl_vector* b1 = NewVectorGSL(b);
    gsl_vector* x = gsl_vector_alloc(b.size());
    gsl_permutation* p = gsl_permutation_alloc(cols);
    int s;
  
    int status = gsl_linalg_LU_decomp(m, p, &s);
    if (status != GSL_SUCCESS)
      CoCoA_THROW_ERROR2(ERR::ExternalLib, gsl_strerror(status));
  
    status = gsl_linalg_LU_solve(m, p, b1, x);
    if (status != GSL_SUCCESS)
      CoCoA_THROW_ERROR2(ERR::ExternalLib, gsl_strerror(status));
  
    const ring R = RingOf(A);
    std::vector<RingElem> rx = NewVector(R, x);

    gsl_permutation_free(p);
    gsl_vector_free(x); 
    gsl_vector_free(b1); 
    gsl_matrix_free(m);
  
    return rx;
  }


  // Cholesky Decomposition 
  matrix GslCholeskyDecomposition(ConstMatrixView A)
  {
    gsl_matrix *m = NewMatrixGSL(A);
  
    int status = gsl_linalg_cholesky_decomp(m);
    if (status != GSL_SUCCESS)
      CoCoA_THROW_ERROR2(ERR::ExternalLib, gsl_strerror(status)); // MEMORY LEAK!  
    const ring R = RingOf(A);
    matrix r = NewMatrix(R, m);
    gsl_matrix_free(m);
    return r;
  }


  /**
   * This function factorizes the M-by-N matrix A into the QR decomposition A = Q R. 
   * On output the diagonal and upper triangular part of the input matrix contain the matrix R.
   * The vector tau and the columns of the lower triangular part of the matrix A contain the 
   * Householder coefficients and Householder vectors which encode the orthogonal matrix Q. 
   * The vector tau must be of length k=\min(M,N). The matrix Q is related to these components 
   * by, Q = Q_k ... Q_2 Q_1 where Q_i = I - \tau_i v_i v_i^T and v_i is the Householder 
   * vector v_i = (0,...,1,A(i+1,i),A(i+2,i),...,A(m,i)). This is the same storage scheme as used by lapack.
   */

  std::vector<matrix> GslQR(ConstMatrixView A)
  {

    const size_t M = NumRows(A);
    const size_t N = NumCols(A);

    gsl_matrix* m = NewMatrixGSL(A);
    gsl_vector* gsl_tau = gsl_vector_alloc(std::min(M,N));
    gsl_matrix* q = gsl_matrix_alloc(M,M);
    gsl_matrix* r = gsl_matrix_alloc(M,N);
  
    // Q is M-by-M and R is M-by-N.
    int status = gsl_linalg_QR_decomp(m, gsl_tau);
    if (status != GSL_SUCCESS)
      CoCoA_THROW_ERROR2(ERR::ExternalLib, gsl_strerror(status)); // MEMORY LEAK!
  
    status = gsl_linalg_QR_unpack (m, gsl_tau, q, r);
    if (status != GSL_SUCCESS)
      CoCoA_THROW_ERROR2(ERR::ExternalLib, gsl_strerror(status)); // MEMORY LEAK!
  
    std::vector<matrix> result;
    const ring R = RingOf(A);
    result.push_back(NewMatrix(R, q));
    result.push_back(NewMatrix(R, r));
  
    gsl_matrix_free(r);
    gsl_matrix_free(q);
    gsl_vector_free(gsl_tau);
    gsl_matrix_free(m);

    return result;
  }


  std::vector<matrix> GslQRPT(ConstMatrixView A)
  {

    const size_t M = NumRows(A);
    const size_t N = NumCols(A);

    gsl_matrix* m = NewMatrixGSL(A);
    gsl_vector* gsl_tau = gsl_vector_alloc(std::min(M,N));
    gsl_matrix* q = gsl_matrix_alloc(M,M);
    gsl_matrix* r = gsl_matrix_alloc(M,N);
    gsl_vector* norm = gsl_vector_alloc(N);

    int s;
    gsl_permutation* p = gsl_permutation_alloc(N);
    int status = gsl_linalg_QRPT_decomp2(m, q, r, gsl_tau, p, &s, norm);
    if (status != GSL_SUCCESS)
      CoCoA_THROW_ERROR2(ERR::ExternalLib, gsl_strerror(status)); // MEMORY LEAK!  
  
    const size_t D = p->size;
    gsl_matrix* tp = gsl_matrix_alloc(1,D);
    for (size_t i = 0; i<D; ++i)
    {
      gsl_matrix_set(tp, 0,i, gsl_permutation_get(p, i));
    }
  
    std::vector<matrix> result;
    const ring R = RingOf(A);
    result.push_back(NewMatrix(R, q));
    result.push_back(NewMatrix(R, r));
    result.push_back(NewMatrix(R, tp));
  
    gsl_permutation_free(p);
    gsl_vector_free(norm);
    gsl_matrix_free(r);
    gsl_matrix_free(q);
    gsl_vector_free(gsl_tau);
    gsl_matrix_free(m);

    return result;
  }


  /**
   * This function solves the square system A x = b using the QR decomposition
   * The least-squares solution for rectangular systems can be found using gsl_linalg_QR_lssolve.
   */

  std::vector<RingElem> GslSolverQR(ConstMatrixView A, std::vector<RingElem>& b)
  {
    gsl_matrix* m = NewMatrixGSL(A);
    gsl_vector* gsl_tau = gsl_vector_alloc(std::min(NumCols(A),NumRows(A)));
    gsl_vector* gsl_b = NewVectorGSL(b);
    gsl_vector* gsl_x = gsl_vector_alloc(gsl_b->size);
  
    int status = gsl_linalg_QR_decomp(m, gsl_tau);
    if (status != GSL_SUCCESS)
      CoCoA_THROW_ERROR2(ERR::ExternalLib, gsl_strerror(status)); // MEMORY LEAK!  
  
    status = gsl_linalg_QR_solve(m, gsl_tau, gsl_b, gsl_x);
    if (status != GSL_SUCCESS)
      CoCoA_THROW_ERROR2(ERR::ExternalLib, gsl_strerror(status)); // MEMORY LEAK!  
  
    const ring R = RingOf(A); 
    std::vector<RingElem> x = NewVector(R, gsl_x);

    gsl_vector_free(gsl_x);
    gsl_vector_free(gsl_b);
    gsl_vector_free(gsl_tau);
    gsl_matrix_free(m);
 
    return x;
  }


  // Probably not exception safe!!
  void QRDecomposition(ConstMatrixView A, matrix& QR, std::vector<RingElem>& tau)
  {
    const size_t cols = NumCols(A);
    const size_t rows = NumRows(A);

    gsl_matrix* m = NewMatrixGSL(A);
    gsl_vector* gsl_tau = gsl_vector_alloc(std::min(cols,rows));
  
    int status = gsl_linalg_QR_decomp(m, gsl_tau);
    if (status != GSL_SUCCESS)
      CoCoA_THROW_ERROR2(ERR::ExternalLib, gsl_strerror(status)); // MEMORY LEAK!  
    
    UpdateMatrix(QR, m);
  
    const ring R = RingOf(A);
  
    const size_t M = gsl_tau->size;
    for (size_t i = 0; i<M; ++i)
    {
      tau.push_back(RingElem(R, ConvertTo<BigRat>(gsl_vector_get(gsl_tau,i))));
    }
  
    gsl_vector_free(gsl_tau);
    gsl_matrix_free(m);
  }


  /**
   * This function factorizes the M-by-N matrix A into bidiagonal form U B V^T. 
   * The diagonal and superdiagonal of the matrix B are stored in the diagonal 
   * and superdiagonal of A. The orthogonal matrices U and V are stored as compressed 
   * Householder vectors in the remaining elements of A. The Householder coefficients 
   * are stored in the vectors tau_U and tau_V. The length of tau_U must equal the number 
   * of elements in the diagonal of A and the length of tau_V should be one element shorter.
   */
  void GslBidiagDecomposition(ConstMatrixView A, matrix& D, std::vector<RingElem>& tau_U, std::vector<RingElem>& tau_V)
  {
    const size_t cols = NumCols(A);
  
    gsl_matrix* m = NewMatrixGSL(A);
    gsl_vector* gsl_tau_U = gsl_vector_alloc(cols);
    gsl_vector* gsl_tau_V = gsl_vector_alloc(cols-1);
  
    int status = gsl_linalg_bidiag_decomp (m, gsl_tau_U, gsl_tau_V);
    if (status != GSL_SUCCESS)
      CoCoA_THROW_ERROR2(ERR::ExternalLib, gsl_strerror(status)); // MEMORY LEAK!  
    
    const ring& R = RingOf(A);
    const size_t M = gsl_tau_U->size;
    for (size_t i = 0; i<M; ++i)
    {
      tau_U.push_back(RingElem(R, ConvertTo<BigRat>(gsl_vector_get(gsl_tau_U,i))));
    }
  
    const size_t N = gsl_tau_V->size;
    for (size_t i = 0; i<N; ++i)
    {
      tau_V.push_back(RingElem(R, ConvertTo<BigRat>(gsl_vector_get(gsl_tau_V,i))));
    }
  
    UpdateMatrix(D, m);
  
    gsl_vector_free(gsl_tau_V);
    gsl_vector_free(gsl_tau_U);
    gsl_matrix_free(m);
  }


  void GslBidiagDecompUnpack(ConstMatrixView A, matrix& U, matrix& V, std::vector<RingElem>& diagonal, std::vector<RingElem>& superdiagonal)
  {
    const size_t M = NumRows(A);
    const size_t N = NumCols(A);
  
    gsl_matrix* m = NewMatrixGSL(A);
    gsl_vector* gsl_tau_U = gsl_vector_alloc(N);
    gsl_vector* gsl_tau_V = gsl_vector_alloc(N-1);
    gsl_matrix* gsl_u = gsl_matrix_alloc(M,N);
    gsl_matrix* gsl_v = gsl_matrix_alloc(N,N);
    gsl_vector* gsl_diagonal = gsl_vector_alloc(N);
    gsl_vector* gsl_superdiagonal = gsl_vector_alloc(N-1);
  
    int status = gsl_linalg_bidiag_decomp (m, gsl_tau_U, gsl_tau_V);
    if (status != GSL_SUCCESS)
      CoCoA_THROW_ERROR2(ERR::ExternalLib, gsl_strerror(status)); // MEMORY LEAK!  
    
    status = gsl_linalg_bidiag_unpack (m, gsl_tau_U, gsl_u, gsl_tau_V, gsl_v, gsl_diagonal, gsl_superdiagonal);
    if (status != GSL_SUCCESS)
      CoCoA_THROW_ERROR2(ERR::ExternalLib, gsl_strerror(status)); // MEMORY LEAK!
    
    UpdateMatrix(U, gsl_u);
    UpdateMatrix(V, gsl_v);

    const ring& R = RingOf(A);
    size_t sz = gsl_diagonal->size;
    for (size_t i = 0; i<sz; ++i)
    {
      diagonal.push_back(RingElem(R, ConvertTo<BigRat>(gsl_vector_get(gsl_diagonal,i))));
    }
  
    sz = gsl_superdiagonal->size;
    for (size_t i = 0; i<sz; ++i)
    {
      superdiagonal.push_back(RingElem(R, ConvertTo<BigRat>(gsl_vector_get(gsl_superdiagonal,i))));
    }

    gsl_vector_free(gsl_superdiagonal);
    gsl_vector_free(gsl_diagonal);
    gsl_matrix_free(gsl_v);
    gsl_matrix_free(gsl_u);
    gsl_vector_free(gsl_tau_V);
    gsl_vector_free(gsl_tau_U);
    gsl_matrix_free(m);
  }

}

#endif
