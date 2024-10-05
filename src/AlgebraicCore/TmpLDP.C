// Copyright (c) 2013  John Abbott,  Anna M. Bigatti
// Original author: 2013  Maria-Laura Torrente

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


#include "CoCoA/ApproxPts2.H"
#include "CoCoA/BigRat.H"
#include "CoCoA/CanonicalHom.H"
#include "CoCoA/DenseMatrix.H"
#include "CoCoA/ExternalLibs-GSL.H"
#include "CoCoA/MatrixOps.H"
#include "CoCoA/MatrixView.H"
#include "CoCoA/QBGenerator.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/SparsePolyOps-RingElem.H"
#include "CoCoA/convert.H"
#include "CoCoA/matrix.H"
#include "CoCoA/ring.H"
#include "CoCoA/utils.H"

using namespace CoCoA;

#include <list>
using std::list;
#include <vector>
using std::vector;
#include <time.h>


// //----------------------------------------------------------------------
// const string ShortDescription =
//   "               Low-degree Polynomial Algorithm (LDP)            \n\n"
//   "Given a set of empirical points, LDP computes (ebar, f, QB) where \n"
//   "ebar=componentwise perturbations of the points,                   \n"
//   "f=(monic) polynomial of low degree,                               \n"
//   "QB=support of f                                                   \n";

// //----------------------------------------------------------------------

namespace CoCoA
{
  namespace ApproxPts
  {
#ifdef CoCoA_WITH_GSL
    typedef std::vector<RingElem> PointR; // all coords must be in same ring 
    //(which must be an ordered field)


// convert vector with rational entries to vector with double entries
    vector<double> convertVector(const vector<RingElem> RatV)
    {
      const int n = len(RatV);
      vector<double> DoubleV(n);
      for (int i=0; i < n; ++i)
      {
        BigRat a;
        IsRational(a,RatV[i]);
        DoubleV[i] = ConvertTo<double>(a);
      }
      return DoubleV;
    }


// convert a double to a RingElem
    void ConvertToRingElem(RingElem& x, double z)
    {
      x = ConvertTo<BigRat>(z);
    }

// convert a RingElem to a double
    double ConvertToDouble(RingElem x)
    {
      BigRat Q; IsRational(Q,x);
      return ConvertTo<double>(Q);
    }

// convert a RingElem to a double
    void ConvertToDouble(double& z, RingElem x)
    {
      z = ConvertToDouble(x);
    }


// // square 2-norm of a vector of doubles
// double SquareNormDouble(const gsl_vector *v)
// {
//   double tmpSqNorm(0);
//   const long n = v->size;
//   for (long i=0; i < n; ++i)
//     {
//       double d = gsl_vector_get(v,i);
//       tmpSqNorm += d*d;
//     }
//   return tmpSqNorm;
// }

// 2-norm of a vector of doubles
    double Norm2Double(const gsl_vector *v)
    {
      double tmpSqNorm(0);
      const long n = v->size;
      for (long i=0; i < n; ++i)
      {
        double d = gsl_vector_get(v,i);
        tmpSqNorm += d*d;
      }
      return sqrt(tmpSqNorm);
    }


// square 2-norm of a vector of RingElem
    RingElem SquareNorm(const vector<RingElem>& v)
    {
      CoCoA_ASSERT(!v.empty());
      RingElem tmpSqNorm(zero(owner(v.front())));
      for (long i=0; i < len(v); ++i)
        tmpSqNorm += power(v[i], 2);
      return tmpSqNorm;
    }

// infinity norm of a vector of RingElem
    RingElem InfinityNorm(const vector<RingElem>& V)
    {
      RingElem InfNorm(zero(owner(V.front())));
      for (long i=0; i < len(V); ++i)
        if (abs(V[i]) > InfNorm) 
          InfNorm = abs(V[i]);
      return InfNorm;
    }


// Evaluation of a poly f on a set Pts of points
    void Evaluation(vector<RingElem>& b,
                    ConstRefRingElem f,
                    const list<RingHom>& EvalAtPts)
    {
      b.clear();
      for (const RingHom& EvalMap: EvalAtPts)
        b.push_back(EvalMap(f));
    }


// Solve the LS problem A*x = b, A orthogonal matrix
// (implemented without deleting the denominators)
    void LeastSquares(vector<RingElem>& x,
                      vector<RingElem>& Rho,
                      const ConstMatrixView& A,
                      const ConstMatrixView& DiagA,
                      const vector<RingElem>& b)
    {
      ring R = owner(b.front());
      const long n = len(b);
      vector<RingElem> RhoRho(n, zero(R));

      x = vector<RingElem>(NumRows(A),zero(R));
      A->myMulByCol(x,b);
      for (long i=0; i < n; ++i)
        if (DiagA(i,i) != 0)
          x[i] /= DiagA(i,i);
        else
          x[i] = zero(R);
      A->myMulByRow(RhoRho,x);
      Rho = b;
      for (long i = 0; i < n; ++i)
        Rho[i] -= RhoRho[i];
    }


// Solve the LS problem M*x = b, M matrix
    void LeastSquares(vector<RingElem>& x,
                      vector<RingElem>& Rho,
                      ConstMatrixView A,
                      const vector<RingElem>& b)
    {
      ring R = owner(b.front());
      const long n = len(b);

      matrix PseudoA = PseudoInverse(A);
      x = vector<RingElem>(NumCols(A),zero(R));
      PseudoA->myMulByCol(x,b);
      Rho.clear();
      Rho = b;
      vector<RingElem> c(n, zero(R));
      A->myMulByCol(c,x);
      for (long i=0;  i < n; ++i)
        Rho[i] -= c[i];
    }

// multiplication of two matrices, the mult matrix is returned 
    gsl_matrix *myMul(const gsl_matrix *A, const gsl_matrix *B)
    {
      long M = A->size1, N = A->size2, K = B->size2;
      gsl_matrix *C = gsl_matrix_alloc(M, K);
      for (long i=0; i < M; ++i)
      {
	for (long j=0; j < K; ++j)
        {
          double d(0);
          for (long l=0; l < N; ++l)
            d += gsl_matrix_get(A, i, l)* gsl_matrix_get(B, l, j);
          gsl_matrix_set (C, i,j, d);
        }
      }
      return C;
    }


// algoritmo per calcolare matrice Jacobiana di Rho in un insieme di punti
    void JacRho(gsl_matrix *JRho0, 
                gsl_vector *Rho0GSL,
                gsl_vector *Alpha0GSL, 
                gsl_matrix *MGSL, 
                ConstRefRingElem tAsPoly,
                const vector<RingElem>& QBAsPoly,
                const list<RingHom>& EvalAtPts,
                const vector<int>& CoordEsatte)
    {
      const SparsePolyRing P = domain(EvalAtPts.front());
      const long dim = NumIndets(P);
      const long NumPts = len(EvalAtPts);
      //const long lenQB = len(QBAsPoly)-1;
      const long lenQB = len(QBAsPoly);
  
      //costruisco B = I-M*PseudoM
      gsl_matrix *Mt = gsl_matrix_alloc(lenQB,NumPts);
      gsl_matrix_transpose_memcpy(Mt, MGSL);
      gsl_matrix *MtM = myMul(Mt, MGSL); 
      gsl_matrix *MtMInv = gsl_matrix_alloc(lenQB,lenQB);
      gsl_permutation * p = gsl_permutation_alloc(lenQB);
      int signum;
      gsl_linalg_LU_decomp(MtM, p, &signum);
      gsl_linalg_LU_invert(MtM, p, MtMInv);
      gsl_matrix *pseudoM = myMul(MtMInv, Mt); // lenQB * NumPts
      gsl_matrix *MpseudoM = myMul(MGSL, pseudoM);
      gsl_matrix *b = gsl_matrix_alloc(NumPts,NumPts);
      gsl_matrix_set_identity(b);
      gsl_matrix_sub(b, MpseudoM);

      for (long i=0; i < dim; ++i)
      {
        int trovato=0;
        int l=0;
        while (l < len(CoordEsatte) && trovato == 0)
	{
          if (i == CoordEsatte[l]) 
            trovato = 1;
          ++l;
	}
        if (trovato == 1)
	{
	  for (int k=0; k < NumPts; ++k)
	    for (long j=0; j < NumPts; ++j)
	      gsl_matrix_set(JRho0, j, i*NumPts+k, 0);
	}
        else
	{
	  RingElem DerT = deriv(tAsPoly,i);
	  vector<RingElem> DerQB;
	  for (long j=0; j < lenQB; ++j)//cambiato!!!
	    DerQB.push_back(deriv(QBAsPoly[j], i));

	  list<RingHom>::const_iterator it = EvalAtPts.begin();
	  RingHom EvalK(*it);
	  for (long k=0; k < NumPts; ++k)
          {
            gsl_matrix *DMQB = gsl_matrix_alloc(NumPts,lenQB);
            gsl_matrix_set_zero(DMQB);
            for (long j=0; j < lenQB; ++j)
            {
              BigRat q; IsRational(q, EvalK(DerQB[j]));
              gsl_matrix_set(DMQB,k,j,ConvertTo<double>(q));
            }

            //costruisco matrice C (vettore colonna)
            gsl_matrix *DT = gsl_matrix_alloc(NumPts,1);
            gsl_matrix_set_zero(DT);
            BigRat q; IsRational(q, EvalK(DerT));
            gsl_matrix_set(DT,k,0,ConvertTo<double>(q));
	      
            gsl_matrix *C = gsl_matrix_alloc(NumPts,1);
            gsl_matrix_set_zero(C);
            gsl_matrix *Alpha = gsl_matrix_alloc(lenQB,1);
            for (long j=0; j < lenQB; ++j)
              gsl_matrix_set(Alpha, j, 0, gsl_vector_get(Alpha0GSL,j));
            gsl_matrix_sub(C, myMul(DMQB, Alpha));
            gsl_matrix_add(C, DT);

            //costruisco matrice E
            gsl_matrix *F = myMul(DMQB,pseudoM);
            gsl_matrix *Ft = gsl_matrix_alloc(NumPts,NumPts);
            gsl_matrix_transpose_memcpy(Ft, F);
            gsl_matrix *Rho = gsl_matrix_alloc(NumPts,1);
            gsl_matrix_set_col(Rho, 0, Rho0GSL);
            gsl_matrix *E = myMul(Ft, Rho);
	 
            gsl_matrix *JJ = myMul(b, C);
            gsl_matrix_sub(JJ, E);
	      
            for (long j=0; j < NumPts; ++j)
              gsl_matrix_set(JRho0, j, i*NumPts+k, gsl_matrix_get(JJ,j,0));
            ++it;
            if (it != EvalAtPts.end())
              EvalK = *it;
          }
	}
      }
    }


// Add to QB and update
    void AddToQBAndUpdate(const PPMonoidElem& t,
                          QBGenerator& QBG)
    {
      QBG.myCornerPPIntoQB(t);
    }


// Update evaluation matrices
    void UpdateMatricesNBM(const QBGenerator& QBG,
                           MatrixView& M,
                           MatrixView& A,
                           MatrixView& DiagA,
                           const vector<RingElem>& V,
                           const vector<RingElem>& Rho)
    {
      const long k = len(QBG.myQB());
      SetEntry(DiagA, k, k, SquareNorm(Rho));
      for (long j=0; j < len(Rho); ++j)
      {
        SetEntry(M, k, j, V[j]);
        SetEntry(A, k, j, Rho[j]);
      }
    }


// Add to Corner set and update
    void AddToCornerAndUpdate(const PPMonoidElem& t,
                              QBGenerator& QBG,
                              const vector<RingElem>& Alpha0,
                              //const vector<RingElem>& S,
                              const vector<RingElem>& QBAsPoly,
                              vector<RingElem>& J)
    {
      const SparsePolyRing P = owner(QBAsPoly.front());
      RingHom RToP = CoeffEmbeddingHom(P);
      const RingElem tAsPoly = monomial(P, t);

      RingElem f(tAsPoly);
      for (long i=0; i < len(QBAsPoly); ++i)
        f -= RToP(Alpha0[i])*QBAsPoly[i];
      J.push_back(f);
      QBG.myCornerPPIntoAvoidSet(t);
    }


// Build the matrix of absolute values
    matrix AbsMatrix(ConstMatrixView M)
    {
      const long nrows = NumRows(M);
      const long ncols = NumCols(M);
      matrix AbsM = NewDenseMat(RingOf(M), nrows, ncols);
      for (long i=0; i < nrows; ++i)
        for (long j=0; j < ncols; ++j)
          SetEntry(AbsM, i, j, abs(M(i,j)));
      return AbsM;
    }


// Build the vector of absolute values
    vector<RingElem> AbsVector(const vector<RingElem>& V)
    {
      vector<RingElem> AbsV;
      for (long i=0; i < len(V); ++i)
        AbsV.push_back(abs(V[i]));
      return AbsV;
    }


// It computes the NBM bound for the sufficient condition...
    void NBMBound(vector<RingElem>& Bound,
                  ConstMatrixView M,
                  const vector<RingElem>& tolerance,
                  const RingElem& g,
                  const list<RingHom>& EvalAtPts)
    {
      const long NumPts = len(EvalAtPts);
      ring R = RingOf(M);
      RingHom phi = CanonicalHom(R, owner(g));
      matrix id = NewDenseMat(IdentityMat(R, NumPts));
      matrix H = transpose(M)*PseudoInverse(transpose(M));
      matrix C = NewDenseMat(IdentityMat(R, NumPts));
      for (long i=0; i < NumPts; ++i)
        for (long j=0; j < NumPts; ++j)
          SetEntry(C, i, j, id(i,j)-H(i,j));
      matrix CAbs = AbsMatrix(C);
      // if (VerbosityLevel() >= 90)
      // {
      // VERBOSE(90) << "CAbs = " << endl;
      // for (long i=0; i < NumPts; ++i)
      //   {
      //     for (long j=0; j < NumPts; ++j)
      // 	VERBOSE(90) << " " << convertToDouble(CAbs(i,j)) << " " ;
      //     VERBOSE(90) << endl;
      //   }
      // }
      vector<RingElem> D(NumPts, zero(R));
      for (long i=0; i < len(tolerance); ++i)
      {
        vector<RingElem> V;
        Evaluation(V, phi(tolerance[i])*deriv(g,i), EvalAtPts);
        vector<RingElem> AbsV = AbsVector(V);
        for (long j=0; j < NumPts; ++j)
          D[j] += AbsV[j];
      }
      Bound.assign(NumPts, zero(R));
      CAbs->myMulByCol(Bound, D);
    }


// the value of the variable check is 0 = False, 1 = True
    void CheckNBMBound(int& Check, 
                       const vector<RingElem>& Rho,
                       ConstMatrixView M,
                       const vector<RingElem>& tolerance,
                       const RingElem& g,
                       const list<RingHom>& EvalAtPts)
    {
      const long NumPts = len(EvalAtPts);
      ring R = RingOf(M);
      vector<RingElem> Bound; 
      NBMBound(Bound, M, tolerance, g, EvalAtPts);
      //vector<double> BoundDouble = convertVector(Bound);
      //cout << "Bound = " << BoundDouble << endl;

      RingElem Diff(R);
      vector<RingElem> AbsRho = AbsVector(Rho);
      //vector<double> AbsRhoDouble = convertVector(AbsRho);
      //cout << "|Rho| = " << AbsRhoDouble << endl;

      for (long i=0; i < NumPts; ++i)
      {
        Diff = AbsRho[i]-Bound[i];
        if (Diff > 0)
          break;
      }
      if (Diff > 0) 
        Check = 1;
      else 
        Check = 0;
    }


// funzione che calcola rango numerico
// valore della funzione in [0,n]
    long NumericalRank(const gsl_matrix* m, 
                       const double& Delta, 
                       const double& k)
    {
      long M = m->size1, N = m->size2;
      gsl_matrix *v = gsl_matrix_alloc(M,N);
      if (M > N)
        gsl_matrix_memcpy(v,m);
      else
      {
        gsl_matrix_free(v);
        gsl_matrix *v = gsl_matrix_alloc(N,M);
        gsl_matrix_transpose_memcpy(v, m);
      }
      M = v->size1; N = v->size2;
      gsl_matrix *q = gsl_matrix_alloc(N,N);
      gsl_vector *sigma = gsl_vector_alloc(N);
      gsl_vector *w = gsl_vector_alloc(N);
      gsl_linalg_SV_decomp(v, q, sigma, w); // sigma contiene valori singolari
///  if (VerbosityLevel() > 90)
///  {
///    VERBOSE(90) << "Valori singolari = " << endl;
///    for (long i=0; i < N; ++i)
///      VERBOSE(90) << " " << gsl_vector_get(sigma, i) << " " ;
///    VERBOSE(90) << endl;
///  }
      long r1, r2, r;
      const long n = sigma->size;
      if (gsl_vector_get(sigma,0) < k*Delta) 
        r1 = 0;
      else
      {
        long i=1;
        while ((i<=n) && (gsl_vector_get(sigma,i-1) > k*Delta))
          ++i;
        r1 = i-1;
      }
      if (gsl_vector_get(sigma, n-1) > Delta) 
        r2 = n+1;
      else
      {
        long i=n;
        while ((i>=1) && (gsl_vector_get(sigma, i-1) < Delta))
          --i;
        r2 = i+1;	
      }
  
      if (r2 == r1+1)
      {
        if (r1 == 0) 
          r = 0; 
        if (r1 == n)
          r = n;
        r = r1;
      }
      else
      {
///      VERBOSE(90) << "non abbiamo potuto considerare la finestra (Delta, k*Delta)" << endl;
        r = r2-1;
      }
      return r;
    }


// funzione che calcola rango numerico, e parametri k e delta
    vector<double> NumericalRank(const gsl_matrix* m,
                                 const double& tol)
    {
      long M = m->size1, N = m->size2;
      gsl_matrix *v = gsl_matrix_alloc(M,N);
      if (M > N)
        gsl_matrix_memcpy(v,m);
      else
      {
        gsl_matrix_free(v);
        gsl_matrix *v = gsl_matrix_alloc(N,M);
        gsl_matrix_transpose_memcpy(v, m);
      }
      M = v->size1; N = v->size2;
      gsl_matrix *q = gsl_matrix_alloc(N,N);
      gsl_vector *sigma = gsl_vector_alloc(N);
      gsl_vector *w = gsl_vector_alloc(N);
      gsl_linalg_SV_decomp(v, q, sigma, w); // sigma contiene valori singolari
///  cout << "Valori singolari = " << endl;
///  for (long i=0; i < N; ++i)
///    cout << " " << gsl_vector_get(sigma, i) << " " ;
///  cout << endl;

      gsl_matrix_free(q); 
      gsl_vector_free(w); 
      gsl_matrix_free(v);
 
      vector<long> indexUp, index, indexDown;
      for (long i=0; i<N; ++i)
      {
        double d = gsl_vector_get(sigma,i);
        if (d >=1) indexUp.push_back(i);
        else 
	{
	  if (d > tol/100) //LAURA 03/05/13
	    index.push_back(i);
	  else
	    indexDown.push_back(i);
	}
      }

///  cout << "indexUp =" << indexUp << endl;
///  cout << "index =" << index << endl;
///  cout << "indexDown =" << indexDown << endl;
 

      double r, Delta, k;
      vector<double> res;
      if (index.empty())  // caso index vuoto
      {
        if (!indexUp.empty() && indexDown.empty())   //caso 1  
          r = N;
        if (indexUp.empty() && !indexDown.empty())   //caso 3
          r = 0;
        if (!indexUp.empty() && !indexDown.empty())  //caso 6
          r = len(indexUp);
        //Delta = 2*tol;
        Delta = tol/10;
        k=1;
      }
      else              // caso index pieno
      {
        //long indMaxDist(0);
        //double MaxDist;
        if (indexUp.empty()) // casi 2 e 5
	{
	  vector<double> dist;
	  for (int i=0; i < len(index)-1; ++i)
	    dist.push_back(gsl_vector_get(sigma,i)-gsl_vector_get(sigma,i+1));
	  dist.push_back(gsl_vector_get(sigma, len(index)-1) - tol);
	  double MaxDist = dist.front();
	  long indMaxDist(0);
	  for (int i=1; i < len(dist); ++i)
	    if (dist[i] > MaxDist)
            {
              MaxDist = dist[i];
              indMaxDist = i;
            }
	  r = indMaxDist + 1;
	  if (indexDown.empty())  //caso 2 
          {
            if (r == N)
            {
              Delta = tol + 0.05*(gsl_vector_get(sigma,N-1)-tol);
              k = (gsl_vector_get(sigma,N-1) - 0.05*(gsl_vector_get(sigma,N-1)-tol))/Delta;
            }
            else
            {
              Delta = gsl_vector_get(sigma,r) + 0.05*(gsl_vector_get(sigma,r-1)-gsl_vector_get(sigma,r));
              k = (gsl_vector_get(sigma,r-1) - 0.05*(gsl_vector_get(sigma,r-1)-gsl_vector_get(sigma,r)))/Delta;
            }
          }
	  else //caso 5
          {
            if (r == N)
            {
              Delta = tol + 0.05*(gsl_vector_get(sigma,N-1)-tol);
              k = (gsl_vector_get(sigma,N-1) - 0.05*(gsl_vector_get(sigma,N-1)-tol))/Delta;
            }
            else
            {
              Delta = gsl_vector_get(sigma,r) + 0.05*(gsl_vector_get(sigma,r-1)-gsl_vector_get(sigma,r));
              k = (gsl_vector_get(sigma,r-1) - 0.05*(gsl_vector_get(sigma,r-1)-gsl_vector_get(sigma,r)))/Delta;
            }
          }
	  if (k>2) k=2;
	}
        else                 // casi 4 e 7
	{
	  vector<double> dist;
	  dist.push_back(1 - gsl_vector_get(sigma, index[0]));
	  for (int i=0; i < len(index)-1; ++i)
	    dist.push_back(gsl_vector_get(sigma,index[i])-gsl_vector_get(sigma,index[i+1]));
	  dist.push_back(gsl_vector_get(sigma, index.back()) - tol);
	  double MaxDist = dist.front();   
	  long indMaxDist(0);
	  for (int i=1; i < len(dist); ++i)
	    if (dist[i] > MaxDist)
            {
              MaxDist = dist[i];
              indMaxDist = i;
            }
	  if (!(indMaxDist == 0)) 
	    r = index[indMaxDist-1] +1;
	  else 
	    r = indexUp.back() + 1; 
	  if (indexDown.empty())  //caso 4 
          {
            if (r == N)
            {
              Delta = tol + 0.05*(gsl_vector_get(sigma,N-1)-tol);
              k = (gsl_vector_get(sigma,N-1) - 0.05*(gsl_vector_get(sigma,N-1)-tol))/Delta;
            }
            if (r == indexUp.back() + 1)
            {
              Delta = gsl_vector_get(sigma,r) + 0.05*(1-gsl_vector_get(sigma,r));
              k = (1 - 0.05*(1-gsl_vector_get(sigma,r)))/Delta;
            }
            if ((r < N) && (r > indexUp.back() + 1))
            {
              Delta = gsl_vector_get(sigma,r) + 0.05*(gsl_vector_get(sigma,r-1)-gsl_vector_get(sigma,r));
              k = (gsl_vector_get(sigma,r-1) - 0.05*(gsl_vector_get(sigma,r-1)-gsl_vector_get(sigma,r)))/Delta;
            }
          }
	  else //caso 7
          {
            if (r == indexUp.back() + 1)
            {
              Delta = gsl_vector_get(sigma,r) + 0.05*(1-gsl_vector_get(sigma,r));
              k = (1 - 0.05*(1-gsl_vector_get(sigma,r)))/Delta;
            }
            else
            {
              Delta = gsl_vector_get(sigma,r) + 0.05*(gsl_vector_get(sigma,r-1)-gsl_vector_get(sigma,r));
              k = (gsl_vector_get(sigma,r-1) - 0.05*(gsl_vector_get(sigma,r-1)-gsl_vector_get(sigma,r)))/Delta;
            }
          }
	  if (k>2) k=2;
	}
      }
      res.push_back(r); 
      res.push_back(Delta); 
      res.push_back(k);
      return res;
    }



//soluzione di minima norma 2 di sistemi sottodeterminati Ax=b 
//(con fattorizzazione QR)
    void MinimalSolution(gsl_vector *xLS,
                         const gsl_matrix *A,
                         const gsl_vector *b)
    {
      //controllare che sia sistema sottodeterminato
      const long m = A->size1;
      const long n = A->size2;

      // QR di At=transpose(A)
      gsl_matrix *at = gsl_matrix_alloc(n,m);
      gsl_matrix_transpose_memcpy(at, A);
 
      gsl_vector *gsl_tau = gsl_vector_alloc(std::min(n,m));
      gsl_matrix *q = gsl_matrix_alloc(n,n);
      gsl_matrix *r = gsl_matrix_alloc(n,m);
      // Q is n-by-n and R is n-by-m
      gsl_matrix *atcp = gsl_matrix_alloc(n,m);
      gsl_matrix_memcpy(atcp, at);
      gsl_linalg_QR_decomp(atcp, gsl_tau);
      gsl_linalg_QR_unpack(atcp, gsl_tau, q, r);
 
      gsl_matrix *q1 = gsl_matrix_alloc(n,m);
      for (long i=0; i < n; ++i)
        for (long j=0; j < m; ++j)
          gsl_matrix_set(q1, i, j, gsl_matrix_get(q, i, j));
      gsl_matrix *r1 = gsl_matrix_alloc(m,m);
      for (long i=0; i < m; ++i)
        for (long j=0; j < m; ++j)
          gsl_matrix_set(r1, i, j, gsl_matrix_get(r, i, j));
 
      gsl_matrix *r1t = gsl_matrix_alloc(m,m);
      gsl_matrix_transpose_memcpy(r1t, r1);
      gsl_matrix *r1tinv = gsl_matrix_alloc(m,m);
      gsl_permutation *p = gsl_permutation_alloc(m);
      int signum;
      gsl_matrix *r1tcp = gsl_matrix_alloc(m,m);
      gsl_matrix_memcpy(r1tcp, r1t);
      gsl_linalg_LU_decomp (r1tcp, p, &signum);
      gsl_linalg_LU_invert(r1tcp, p, r1tinv);

      gsl_matrix *bmat = gsl_matrix_alloc(m,1);
      gsl_matrix_set_col(bmat, 0, b);
      gsl_matrix *y1=myMul(r1tinv,bmat);
      gsl_matrix *x2=myMul(q1,y1);
      gsl_matrix_get_col(xLS, x2, 0);

      gsl_matrix *xLSmat = gsl_matrix_alloc(n,1);
      gsl_matrix_set_col(xLSmat, 0, xLS);

      gsl_vector_free(gsl_tau);
      gsl_matrix_free(q); gsl_matrix_free(r);
      gsl_matrix_free(q1); gsl_matrix_free(r1);
      gsl_matrix_free(r1t); gsl_matrix_free(r1tinv);
      gsl_matrix_free(bmat); gsl_matrix_free(y1);
      gsl_matrix_free(x2); //gsl_permutation_free(p);
    }


// aggiornare i punti
    void UpdatePts(std::vector<PointR>& newpts, 
                   const gsl_vector *h)
    {
      long NumPts = len(newpts);
      long dim = len(newpts.front());
      for (long i=0; i < NumPts; ++i)
        for (long j=0; j < dim; ++j)
        {
          newpts[i][j] += ConvertTo<BigRat>(gsl_vector_get(h,j*NumPts+i));
        }
    }

// aggiornare i punti
    void UpdatePts(std::vector<PointR>& newpts, 
                   const vector<PointR>& pts, 
                   const vector<RingElem>& h)
    {
      newpts=pts;
      long NumPts = len(newpts);
      long dim = len(newpts.front());
      for (long i=0; i < NumPts; ++i)
        for (long j=0; j < dim; ++j)
          newpts[i][j] += h[j*NumPts+i];
    }


//normal flow per rho
    void NormalFlow(gsl_vector *sol,
                    const std::vector<PointR>& pts, 
                    const gsl_matrix *PiR, 
                    const gsl_matrix *JRho0, 
                    const gsl_vector *Rho0GSL,  
                    ConstRefRingElem tAsPoly, 
                    const vector<RingElem>& QBAsPoly, 
                    const double omega, 
                    const double /*tol*/, 
                    const double /*Delta*/, 
                    const double /*k*/,
                    const vector<int>& CoordEsatte)
    {
      const SparsePolyRing P = owner(tAsPoly);
      const ring R = CoeffRing(P);
      const long dim = len(pts.front());
      const long NumPts = Rho0GSL->size;
      const long lenQB = len(QBAsPoly);
      const long r = PiR->size2;

      vector<long> rows;
      for (long i=0; i < r; ++i)
      {
        rows.push_back(ConvertTo<long>(gsl_matrix_get(PiR,0,i))); 
      }
      gsl_vector *RR = gsl_vector_alloc(r);
      for (long i=0; i < r; ++i)
        gsl_vector_set(RR, i, -gsl_vector_get(Rho0GSL, rows[i]));

      gsl_matrix *JR = gsl_matrix_alloc(r, NumPts*dim);
      for (long i=0; i < r; ++i)
      {
        gsl_vector *v = gsl_vector_alloc(NumPts*dim);
        gsl_matrix_get_row(v, JRho0, ConvertTo<long>(gsl_matrix_get(PiR,0,i)));
        gsl_matrix_set_row(JR, i, v);
      }
      vector<PointR> newpts = pts;  ///// va bene????
      gsl_vector *h=gsl_vector_alloc(dim*NumPts);
      for (long i=0; i < dim*NumPts; ++i)
        gsl_vector_set(h,i,1);
      gsl_vector_set_zero(sol);
  
      long passi(0);
      //while ((SquareNormDouble(h) > omega*omega) && (gsl_vector_max(sol) < tol) && (NumericalRank(JR, Delta, k) == r))
      //while ((SquareNormDouble(h) > omega*omega) && (gsl_vector_max(sol) < tol))
      while ((Norm2Double(h) > omega) && passi < 30)
      {
        ++passi;
///      VERBOSE(95) << endl << "ITERAZIONE N^ " << passi << endl;

        MinimalSolution(h, JR, RR);
///     if (VerbosityLevel() >= 95)
///     {
///       VERBOSE(95) << "h = ";
///       for (long i=0; i < dim*NumPts; ++i)
///	     VERBOSE(95) << " " << gsl_vector_get(h,i) << " ";
///       VERBOSE(95) << endl;
///     }
        for (long i=0; i < dim*NumPts; ++i)
	{
	  double d = gsl_vector_get(sol,i);
	  gsl_vector_set(sol, i, d + gsl_vector_get(h,i));
	}
        UpdatePts(newpts, h);
        list<RingHom> EvalAtPts;
        for (long i=0; i < NumPts; ++i)
          EvalAtPts.push_back(EvalHom(P, newpts[i]));
        //creo M e b
        vector<RingElem> b, Alpha, Rho;
        Evaluation(b, tAsPoly, EvalAtPts);
        matrix M(NewDenseMat(R, NumPts, lenQB));
        for (long j=0; j < lenQB; ++j)
	{
	  vector<RingElem> V;
	  Evaluation(V, QBAsPoly[j], EvalAtPts);
	  for (long i=0; i < NumPts; ++i)
	    SetEntry(M, i, j, V[i]);
	}

        //calcolo Rho, JRho
        LeastSquares(Alpha, Rho, M, b);
        //conversioni
        gsl_vector *RhoGSL = NewVectorGSL(Rho);
        gsl_vector *AlphaGSL = NewVectorGSL(Alpha);
        MatrixView MM = transpose(transpose(M));
        gsl_matrix *MGSL = NewMatrixGSL(MM);
        gsl_matrix *JRho = gsl_matrix_alloc(NumPts, NumPts*dim);
        JacRho(JRho, RhoGSL, AlphaGSL, MGSL, tAsPoly, QBAsPoly, EvalAtPts, CoordEsatte);
        // e ora calcolo RR, JR
        for (long i=0; i < r; ++i)
	{
	  double d = gsl_vector_get(RhoGSL, rows[i]);
	  gsl_vector_set(RR, i, -d);
	}
        for (long i=0; i < len(rows); ++i)
          for (long j=0; j < NumPts*dim; ++j)
            gsl_matrix_set(JR, i, j, gsl_matrix_get(JRho, rows[i], j));
      }
    }


// algoritmo Root Finding 
    void RootFindingAlgorithm(vector<RingElem>& sol,
                              vector<RingElem>& Rho0, 
                              vector<RingElem>& Alpha0,
                              MatrixView& M0,
                              ConstRefRingElem tAsPoly,
                              const vector<RingElem>& QBAsPoly, 
                              const list<RingHom>& EvalAtPts0,
                              ConstRefRingElem tol,
                              const vector<PointR>& pts,
                              const vector<int>& CoordEsatte) 
    {
      const long NumPts = len(Rho0);
      const long dim = len(pts.front());
      ring R = owner(Rho0.front());
      //conversioni
      gsl_vector *Rho0GSL = NewVectorGSL(Rho0);
      gsl_vector *Alpha0GSL = NewVectorGSL(Alpha0);
      gsl_matrix *MGSL = NewMatrixGSL(M0);
      BigRat qq; IsRational(qq,tol);
      double tolDouble = ConvertTo<double>(qq);

      //calcolo lo Jacobiano in zero e costruisco (JRho0 | Rho0)
      gsl_matrix *JRho0 = gsl_matrix_alloc(NumPts, NumPts*dim);
      JacRho(JRho0, Rho0GSL, Alpha0GSL, MGSL, tAsPoly, QBAsPoly, EvalAtPts0, CoordEsatte);

      gsl_matrix *JRhoRho = gsl_matrix_alloc(NumPts, NumPts*dim+1);
      for (long j=0; j<NumPts*dim; ++j)
      {
        gsl_vector *v = gsl_vector_alloc(NumPts);
        gsl_matrix_get_col(v, JRho0, j);
        gsl_matrix_set_col(JRhoRho, j, v);
      }
      gsl_matrix_set_col(JRhoRho, NumPts*dim, Rho0GSL);


      //in piu //calcolo ||JRhoRho|| e poi normalizzo....
      // gsl_matrix *v = gsl_matrix_alloc(NumPts*dim+1,NumPts); //v copia di JRhoRho^t
//   gsl_matrix_transpose_memcpy(v, JRhoRho);
//   long M = v->size1; long N = v->size2;
//   gsl_matrix *q = gsl_matrix_alloc(N,N);
//   gsl_vector *sigma = gsl_vector_alloc(N);
//   gsl_vector *w = gsl_vector_alloc(N);
//   gsl_linalg_SV_decomp(v, q, sigma, w); // sigma = SingValues di JRhoRho
//   double coeff = gsl_vector_get(sigma, 0);
//   cout << "Norma JRhoRho = " << coeff << endl;
//   cout << "Valori singolari = " << endl;
//   for (long i=0; i < N; ++i)
//     cout << " " << gsl_vector_get(sigma, i) << " " ;
//   cout << endl;
//   gsl_matrix_free(q);
  
      //in piu


      //calcolo rango numerico di (JRho0 | Rho0) 
      // const double k=2;
      //   const double Delta = (2.5)*tolDouble;   //tol da convertire!!
      //   long rango = NumericalRank(JRhoRho, Delta, k);
  

      //calcolo rango numerico di (JRho0 | Rho0) 
      //con scelta automatica di k e Delta

      //in piu //normalizzo JRhoRho
      //gsl_matrix_scale(JRhoRho, 1/coeff);
      vector<double> CompleteRank = NumericalRank(JRhoRho, tolDouble); //rango, delta, k
      long rango = CompleteRank[0];
      double Delta = CompleteRank[1];
      double k = CompleteRank[2];
///  cout << endl << "Delta = " << Delta << endl;
///  cout << "k = " << k << endl;
///  cout << "k*Delta = " << k*Delta << endl;
///  cout << "Rango numerico di (JRho(0) | Rho(0)) = " << rango << endl;

      if (rango == 0) 
      {
///    cout << "valore sballato di Delta" << endl;
///    cout << "la matrice (JRho|Rho0) ha rango numerico nullo!!" << endl;
      }

      // scelta delle righe in (JRho(0) | Rho0(0))
      gsl_matrix *JRhoRhot = gsl_matrix_alloc(NumPts*dim+1, NumPts); 
      gsl_matrix_transpose_memcpy(JRhoRhot, JRhoRho);
      gsl_vector *gsl_tau = gsl_vector_alloc(std::min(NumPts*dim+1, NumPts));
      gsl_matrix *q1 = gsl_matrix_alloc(NumPts*dim+1, NumPts*dim+1);
      gsl_matrix *r = gsl_matrix_alloc(NumPts*dim+1, NumPts);
      gsl_vector *norm = gsl_vector_alloc(NumPts);
      int s;
      gsl_permutation * p = gsl_permutation_alloc(NumPts);
      gsl_linalg_QRPT_decomp2(JRhoRhot, q1, r, gsl_tau, p, &s, norm);
///  cout << "Permutazione pi greco = " << endl;
///  for (long i=0; i < NumPts; ++i)
///    cout << " " << gsl_permutation_get(p, i) << " " ;
///  cout << endl;
      //permutazione di r righe che ci interessano 
      gsl_matrix *PiR = gsl_matrix_alloc(1,rango);
      for (int i = 0; i < rango; ++i)
        gsl_matrix_set(PiR, 0,i, gsl_permutation_get(p, i));
  
      gsl_permutation_free(p);
      gsl_vector_free(gsl_tau);
      gsl_matrix_free(q1);gsl_matrix_free(r);
  
      // Normal Flow
///  cout << endl << " --- Normal Flow --- " << endl << endl;
      //double omega = 0.0000000000000001;
      constexpr double omega = 0.0000000000000000000000000001;
      gsl_vector *x = gsl_vector_alloc(NumPts*dim); 
      // in piu //normalizzo JRho0, Rho0GSL
      //gsl_matrix_scale(JRho0, 1/coeff);
      //gsl_vector_scale(Rho0GSL, 1/coeff);
      NormalFlow(x, pts, PiR, JRho0, Rho0GSL, tAsPoly, QBAsPoly, omega, tolDouble, Delta, k, CoordEsatte);
      sol.clear();
      for (long i=0; i < NumPts*dim; ++i)
      {
        double d = gsl_vector_get(x,i);
        RingElem elem(R, ConvertTo<BigRat>(d));
        sol.push_back(elem);
      }
    }



    void ComputeGamma(gsl_vector *gamma, 
                      const std::vector<RingElem>& JacF, 
                      const std::vector<PointR>& newpts,
                      const vector<RingElem>& smallR)
    {
      const SparsePolyRing P = owner(JacF.front());
      const long dim = len(newpts[0]);
      const long NumPts = len(smallR);
      srand (time(nullptr));

      for (long int i=0; i < NumPts; ++i)
      {
        for (long numIter=0; numIter < 100; ++numIter)
	{
	  vector<RingElem> pt1(newpts[i]);
	  vector<RingElem> pt2(newpts[i]);
	  vector<RingElem> diffPts;
	  for (long int j=0; j < dim; ++j)
          {
            pt1[j] += (rand() % 2000 -1000)*(smallR[i]/1000);
            pt2[j] += (rand() % 2000 -1000)*(smallR[i]/1000);
            diffPts.push_back(pt1[j] - pt2[j]);
          }
	  // Build evaluation homomorphisms at the newpts + Randspost
	  RingHom Eval1(EvalHom(P, pt1));
	  RingHom Eval2(EvalHom(P, pt2));

	  //compute the jacobians at the new random points, 
	  //compute their difference 
	  vector<RingElem> EvalJac1, EvalJac2;
	  for (long int j=0; j < dim; ++j)
          {
            EvalJac1.push_back(Eval1(JacF[j]));
            EvalJac2.push_back(Eval2(JacF[j]));
          }
	  vector<RingElem> diffJac(EvalJac1);
	  for (long int j=0; j < dim; ++j)
	    diffJac[j] -= EvalJac2[j];

	  //conversioni
	  gsl_vector *diffJacGSL = NewVectorGSL(diffJac);
	  gsl_vector *diffPtsGSL = NewVectorGSL(diffPts);

	  double NewGamma=Norm2Double(diffJacGSL)/Norm2Double(diffPtsGSL);
	  //cout << "Num = " << Norm2Double(diffJacGSL) << endl;
	  //cout << "Den = " << Norm2Double(diffPtsGSL) << endl;
	  //cout << "NewGamma = " << NewGamma << endl; 
	  // RingElem NewGamma(Norm(diffJac));
// 	  NewGamma /= Norm(diffPts);
	  if (NewGamma > gsl_vector_get(gamma,i)) 
            gsl_vector_set(gamma,i,NewGamma);
	}
      }
    }


    void ComputeCapitalR(gsl_vector *capitalR, 
                         const std::vector<RingElem>& JacF, 
                         const std::vector<PointR>& newpts,
                         const std::vector<RingElem>& smallR,
                         const gsl_vector *gamma)
    {
      const SparsePolyRing P = owner(JacF.front());
      const long dim = len(newpts[0]);
      const long NumPts = len(smallR);

      for (long int i=0; i < NumPts; ++i)
      {
        RingHom Eval(EvalHom(P, newpts[i]));
        vector<RingElem> EvalJac;
        for (long int j=0; j < dim; ++j)
          EvalJac.push_back(Eval(JacF[j]));
        //conversione
        gsl_vector *EvalJacGSL = NewVectorGSL(EvalJac);

        double tmpR=Norm2Double(EvalJacGSL);
        tmpR/=gsl_vector_get(gamma,i);
        if (ConvertToDouble(smallR[i]) < tmpR) 
          gsl_vector_set(capitalR,i,ConvertToDouble(smallR[i]));
        else 
          gsl_vector_set(capitalR,i,tmpR);
      }
    }


    void ComputeMu(gsl_vector *Mu, 
                   const std::vector<RingElem>& JacF, 
                   const std::vector<PointR>& newpts,
                   gsl_vector *capitalR)
    {
      const SparsePolyRing P = owner(JacF.front());
      const ring R = owner(newpts[0][0]);
      const long dim = len(newpts[0]);
      const long NumPts = capitalR->size;
      srand (time(nullptr));

      //conversione
      vector<RingElem> capitalRCoc = NewVector(R, capitalR);

      for (long int i=0; i < NumPts; ++i)
      {
        for (long numIter=0; numIter < 100; ++numIter)
	{
	  vector<RingElem> randPt(newpts[i]);
	  for (long int j=0; j < dim; ++j)
	    randPt[j] += (rand() % 2000 -1000)*(capitalRCoc[i]/1000);
	  // Build evaluation homomorphisms at the newpt + Randspost
	  RingHom Eval(EvalHom(P, randPt));
	  
	  //compute the jacobians at the new random point
	  vector<RingElem> EvalJac;
	  for (long int j=0; j < dim; ++j)
	    EvalJac.push_back(Eval(JacF[j]));
 
	  gsl_vector *EvalJacGSL = NewVectorGSL(EvalJac);
	  double NewMu=1/Norm2Double(EvalJacGSL);
	  if (NewMu > gsl_vector_get(Mu,i)) 
	    gsl_vector_set(Mu,i,NewMu);
	}
      }
    }


    void ComputeChi(gsl_vector *Chi, 
                    const gsl_vector *capitalR, 
                    const gsl_vector *gamma,
                    const gsl_vector *Mu)
    {
      const long NumPts = capitalR->size;
      for (long int i=0; i < NumPts; ++i)
      {
        double Ri = gsl_vector_get(capitalR,i);
        double Mui = gsl_vector_get(Mu,i);
        double Gammai = gsl_vector_get(gamma,i);
        gsl_vector_set(Chi, i, Ri/(Mui*(2+Gammai*Ri*Mui)));
      }
    }



    void KantorovichCheck(int& KantCheck,
                          ConstRefRingElem tol,
                          const std::vector<RingElem>& sol,
                          const std::vector<PointR>& newpts,
                          const std::vector<RingElem>& LowDegreePoly)
    //const std::vector<PPMonoidElem>& QB)
    {
      ring R = owner(newpts[0][0]);
      const SparsePolyRing P = owner(LowDegreePoly.front());
      const long NumPts = len(newpts);
      const long dim = len(newpts[0]);

      //costruisco Jacobiano
      vector<RingElem> JacF;
      for (long int i=0; i < dim; ++i)
        JacF.push_back(deriv(LowDegreePoly.front(),i));  
 
      //calcolo r_i
      vector<RingElem> smallR(NumPts, tol);
      //gsl_vector *smallR = gsl_vector_alloc(NumPts);
      for (long i=0; i < NumPts; ++i)
      {
        vector<RingElem> spost(dim, zero(R));
        for (long j=0; j < dim; ++j)
          spost[j] += sol[j*NumPts+i];
        //gsl_vector_set(smallR, i, tolDouble - convertToDouble(InfinityNorm(spost)));
        smallR[i] -= InfinityNorm(spost);
      }

      //calcolo gamma_i
      gsl_vector *gamma = gsl_vector_alloc(NumPts);
      gsl_vector_set_zero(gamma);
      ComputeGamma(gamma, JacF, newpts, smallR);

      //calcolo R_i
      gsl_vector *capitalR = gsl_vector_alloc(NumPts);
      gsl_vector_set_zero(capitalR);
      ComputeCapitalR(capitalR, JacF, newpts, smallR, gamma);

      //calcolo \mu_i
      gsl_vector *Mu = gsl_vector_alloc(NumPts);
      gsl_vector_set_zero(Mu);
      ComputeMu(Mu, JacF, newpts, capitalR);

      //calcolo \chi_i
      gsl_vector *Chi = gsl_vector_alloc(NumPts);
      gsl_vector_set_zero(Chi);
      ComputeChi(Chi, capitalR, gamma, Mu);

      //Stampa momentanea per verifica
///  cout << "r_i = ";
///  PrintVector(NewVectorGSL(smallR));
///  cout << "gamma = ";
///  PrintVector(gamma);
///  cout << "R_i = ";
///  PrintVector(capitalR);
///  cout << "Mu_i = ";
///  PrintVector(Mu);
///  cout << "Chi_i = ";
///  PrintVector(Chi);

      //controllo finale
      KantCheck = 1;
      //cout << "LowDegreePoly = " << LowDegreePoly.back() << endl;
///  cout << "|f_i(Ebar_i)| = ";  
      for (long int i=0; i < NumPts; ++i)
      {
        RingHom Eval(EvalHom(P, newpts[i]));
        RingElem EvalF(Eval(LowDegreePoly.front()));
        //cout << endl << "punto i = " << newpts[i] << endl;
///      cout << ConvertToDouble(abs(EvalF)) << " ";
        if (ConvertToDouble(abs(EvalF)) > gsl_vector_get(Chi,i)) 
          KantCheck = 0;
      }
///  cout << endl;

///  cout << "------------------------------------------" << endl;
      if (KantCheck == 0)
      {
///    cout << "Kantorovich non verificato" << endl;
      }
      else
      {
///      cout << "Kantorovich verificato" << endl << endl;
        double MaxR=gsl_vector_get(capitalR,0);
        for (long int k=0; k < NumPts; ++k)
          if (gsl_vector_get(capitalR,k) > MaxR) 
            MaxR = gsl_vector_get(capitalR,k);
///      cout << "Esiste una perturbazione e* dei punti su cui " << endl;
///      cout << "il polinomio si annulla. Tale perturbazione e*" << endl;
///      cout << "verifica la relazione Inf_norm(e*-EBar) < " << MaxR << endl;
      }  
    }

#endif

///programma principale
    void VanishPoly(std::vector<PointR>& newpts,
                    std::vector<RingElem>& LowDegreePoly,
                    const SparsePolyRing& P,
                    const std::vector<PointR>& pts,
                    const vector<RingElem>& tolVector,
                    ConstRefRingElem tol)
    {
#ifndef CoCoA_WITH_GSL
      CoCoA_THROW_ERROR("ExternalLib GSL not linked: VanishPoly undefined","VanishPoly");
      (void)(newpts);        // just to avoid tedious compiler warning
      (void)(LowDegreePoly); //
      (void)(P);             //
      (void)(pts);           //
      (void)(tolVector);     //
      (void)(tol);           //
#else
      // JAA cannot be bothered to check for "mixed rings"...
      if (pts.empty()) CoCoA_THROW_ERROR("No points","LowDegApproxPoly");
      const long NumPts = len(pts);
      const long dim = len(pts[0]);
      if (dim == 0) CoCoA_THROW_ERROR("dimension zero", "LowDegApproxPoly");
      for (long i=1; i < NumPts; ++i)
        if (len(pts[i]) != dim) CoCoA_THROW_ERROR("mixed dimensions", "LowDegApproxPoly");
      if (len(tolVector) != dim) CoCoA_THROW_ERROR("tolVector wrong size", "LowDegApproxPoly");
      ring R = owner(pts[0][0]);                      
      CoCoA_ASSERT(IsOrderedDomain(R) && IsField(R)); // R must be an ordered field
  
      // Build poly ring P=Q[x_1...x_n]
      // ANNA
      PPMonoid PPMon = PPM(P);
      //   PPOrdering PPO = NewStdDegRevLexOrdering(dim);
//   PPMonoid PPM = NewPPMonoid(SymbolRange("x", 0, dim-1), PPO);
//   SparsePolyRing P = NewPolyRing(R, PPM);
      RingHom coeff = CoeffEmbeddingHom(P);
  
      // Build evaluation homomorphisms at Pts
      list<RingHom> EvalAtPts;
      for (long i=0; i < NumPts; ++i)
        EvalAtPts.push_back(EvalHom(P, pts[i]));
  
      //   --------------     INITIALIZE   ---------------
      long IterNum = 1;
      vector<RingElem> J, S;     // J
  
      QBGenerator QBG(PPMon);      // to handle quotient bases
      QBG.myCornerPPIntoQB(one(PPMon)); 
      //S.push_back(one(P));
      RingElem tAsPoly(one(P));
      list<PPMonoidElem> L = QBG.myCorners();  // L = list of corner PPs
  
      // Build and initialize evaluation matrix M (by rows),
      // its orthogonal part A, and DiagA = A*A'
      // Compute eval V=t(Pts) and solve LS problem A'*b = V
      vector<RingElem> V;
      Evaluation(V, tAsPoly, EvalAtPts);

      matrix A(NewDenseMat(R, NumPts, NumPts));
      for (long i=0; i < NumPts; ++i)
        //SetEntry(A, 0, i, 1);
        SetEntry(A, 0, i, V[i]);
  
      matrix M(NewDenseMat(R, NumPts, NumPts));
      for (long i=0; i < NumPts; ++i)
        //SetEntry(M, 0, i, 1);
        SetEntry(M, 0, i, V[i]);
  
      vector<RingElem> DiagAElems(NumPts, zero(R));
      MatrixView DiagA(DiagMat(DiagAElems));//DiagA=A*A'
      //SetEntry(DiagA, 0, 0, NumPts);
      SetEntry(DiagA, 0, 0, SquareNorm(V));
  
      vector<RingElem> sol(dim*NumPts, zero(R));
  
      vector<int> CoordEsatte;
      for (int i=0; i < dim; ++i)
        if (tolVector[i]==0) 
          CoordEsatte.push_back(i);
///  cout << endl << "Coordinate esatte = " << CoordEsatte << endl << endl;
  

      // ------------   Main CYCLE   -----------------
      while (J.empty())  
        //while (!(L.empty()))
      {
        const vector<PPMonoidElem> QB = QBG.myQB();
///      clog << "---------------------------------------" << endl
///	   << "Iteration " << IterNum << endl
///	   << "Quotient Basis = " << QB << endl
///	   << "List of power products = " << L << endl;
      
        const PPMonoidElem t = L.front(); // PP considered in the current step
///      clog << "Next PP = " << t << endl;
        const RingElem tAsPoly = monomial(P, t);
        vector<RingElem> QBAsPoly;
        for (long i=0; i < len(QB); ++i)
          QBAsPoly.push_back(monomial(P, QB[i]));

        // Compute eval V=t(Pts) and solve LS problem A'*b = V
        vector<RingElem> V, Alpha, Rho;
        Evaluation(V, tAsPoly, EvalAtPts);
        //LeastSquares(Alpha, Rho, A, DiagA, V);
        vector<long> cols, rows;
        for (long i=0; i < len(QBAsPoly); ++i)  rows.push_back(i);
        // rows = LongRange(0, len(QBAsPoly)-1);
        for (long j=0; j < NumCols(M); ++j)  cols.push_back(j);
        // cols = LongRange(0, NumCols(M)-1);
        MatrixView SubM = submat(M, rows, cols);
        MatrixView SubMt = transpose(SubM);
        LeastSquares(Alpha, Rho, SubMt, V);
      
        // First trivial check
        if (IsZero(Rho))
        {
          // AddToCornerAndUpdate(t, QBG, Alpha, S, J);
	  AddToCornerAndUpdate(t, QBG, Alpha, QBAsPoly, J);
///          clog << "Truly vanishing polynomial f = " << J.back() << endl;
	  //aggiungere anche qui ..... 
          ++IterNum;
          L = QBG.myCorners();
          continue;
        }
      
        // build the poly g = t - sum(alpha_i s_i)
        RingElem g(tAsPoly);
        //cout << "Alpha = " << Alpha << endl;
        //cout << "S = " << S << endl;
        //      for (long i=0; i < len(S); ++i)
        //g -= coeff(Alpha[i])*S[i];
        for (long i=0; i < len(QBAsPoly); ++i)
          g -= coeff(Alpha[i])*QBAsPoly[i];
      
        // vector<long> cols, rows;
        // for (long i=0; i < len(S); ++i)
        // 	rows.push_back(i);
        // for (long j=0; j < NumCols(M); ++j)
        // 	cols.push_back(j);
        // MatrixView SubM = submat(M, rows, cols);
        // MatrixView SubMt = transpose(SubM);

        // check whether Rho satisfies condition NBM
        int Check;
        //cout << "g = " << g << endl;
        CheckNBMBound(Check, Rho, SubM, tolVector, g, EvalAtPts);

        // if Check = 1 then |Rho[i]| - NBMBound[i]>0, and then ...
///      cout << endl << "NBM check = " << Check << endl << endl; 
        if (Check == 1)
        {
          UpdateMatricesNBM(QBG, M, A, DiagA, V, Rho);
          // AddToQBAndUpdate(t, QBG, Alpha, S);
	  AddToQBAndUpdate(t, QBG);
	  L = QBG.myCorners();
        }
        else
        {
	  /// parte nuova ......
///	  cout << endl << " --- Root Finding --- " << endl << endl;
	  vector<RingElem> AlphaM, RhoM, xbar;
	  LeastSquares(AlphaM, RhoM, SubMt, V);
	  RootFindingAlgorithm(xbar, RhoM, AlphaM, SubMt, tAsPoly, QBAsPoly, EvalAtPts, tol, pts, CoordEsatte);
	  if ((IsZero(xbar) && !IsZero(Rho)) || (InfinityNorm(xbar) > tol))
          {
            UpdateMatricesNBM(QBG, M, A, DiagA, V, Rho);
            //AddToQBAndUpdate(t, QBG, Alpha, S);
            AddToQBAndUpdate(t, QBG);
            L = QBG.myCorners();
          }
	  else
          {
            //AddToCornerAndUpdate(t, QBG, Alpha, S, J); ///non posso più usarla!!!
            sol=xbar;
///	      cout << endl;
///	      cout << endl << "Spostamento totale= " << convertVector(sol) << endl;
            UpdatePts(newpts, pts, sol); 
///	      cout << "New Points  = " << endl;
///	      for (long i=0; i < NumPts; ++i)
///		cout << convertVector(newpts[i]) << endl;
            list<RingHom> EvalAtNewPts;
            for (long i=0; i < NumPts; ++i)
              EvalAtNewPts.push_back(EvalHom(P, newpts[i]));
	      
            //creo M e b
            vector<RingElem> Newb, NewAlpha, NewRho;
            Evaluation(Newb, tAsPoly, EvalAtNewPts);
            matrix NewM(NewDenseMat(R, NumPts, len(QB)));
            for (long j=0; j < len(QB); ++j)
            {
              vector<RingElem> V;
              Evaluation(V, QBAsPoly[j], EvalAtNewPts);
              for (long i=0; i < NumPts; ++i)
                SetEntry(NewM, i, j, V[i]);
            }

            //calcolo Rho, JRho, f
            LeastSquares(NewAlpha, NewRho, NewM, Newb);
            RingElem f(tAsPoly);
            for (long i=0; i < len(QB); ++i)
              f -= coeff(NewAlpha[i])*QBAsPoly[i];
            //stampa di f
///	      cout << endl << "Polynomial f = " << tAsPoly;
///	      for (int i=len(QB)-1; i >= 0; --i)
///		{
///		  BigRat a; double d; //RingElem c = coeff(NewAlpha[i]);
///		  IsRational(a, -NewAlpha[i]); d=ConvertTo<double>(a);
///		  cout << " +" << d << "*" << QBAsPoly[i];
///		}
///	      cout << endl << endl;
            //vedere se è almost vanishing
////??            double VV=Norm2Double(NewVectorGSL(NewRho));
            //cout << "||f(Xbar)|| = " << VV << endl;
            vector<RingElem> CoeffsF;
            CoeffsF.push_back(one(R));
            for (int i=0; i < len(NewAlpha); ++i)
              CoeffsF.push_back(coeff(NewAlpha[i]));
            //cout << "||f|| = " << Norm2Double(NewVectorGSL(CoeffsF)) << endl;
////??            VV /= Norm2Double(NewVectorGSL(CoeffsF)); 
            ///Norm2Double(NewVectorGSL(NewAlpha));
///	      cout << "||f(Xbar)||_2/||f|| = " << VV << endl;
            //J.push_back(f);


            //LAURA
            vector<RingElem> JJ;
            JJ.push_back(f);
            int KantCheck(1);
            KantorovichCheck(KantCheck,tol,sol,newpts,JJ);
            QBG.myCornerPPIntoAvoidSet(t);
            if (KantCheck == 0) 
            {
              CoCoA_THROW_ERROR("Method failed", "LowDegApproxPoly");
              //UpdateMatricesNBM(QBG, M, A, DiagA, V, Rho);
              //\AddToQBAndUpdate(t, QBG, Alpha, S);
              //AddToQBAndUpdate(t, QBG);
              L = QBG.myCorners();
              //J.push_back(f);
            }
            else
              //LAURA	
            {	
              // QBG.myCornerPPIntoAvoidSet(t);
              J.push_back(f);
              //cout << "Vanishing polynomial = " << J << endl;
            }
          }
        }     
        ++IterNum;
        L = QBG.myCorners();
      }
      std::swap(LowDegreePoly, J); // morally AlmostVanishing = J;
///  cout << "Vanishing polynomial = " << LowDegreePoly.back() << endl;

      //Kantorovich check
      //int KantCheck(1);
      //KantorovichCheck(KantCheck,tol,sol,newpts,LowDegreePoly);
#endif
    }

  } // end of namespace ApproxPts
} // end of namespace CoCoA

// void program()
// {
//   GlobalManager CoCoAFoundations;

//   cout << ShortDescription << endl;

//   std::vector<PointR> newpts;    
//   std::vector<RingElem> LowDegreePoly;
  
//   cout << "Insert the space dimension: " << endl;
//   long dim;
//   cin >> dim;
//   if ( !cin )
//     CoCoA_THROW_ERROR("Input must be a positive integer", "program");
//   if (dim < 1 || dim > 1000000)
//     CoCoA_THROW_ERROR("Ridiculous input dimension", "main program in ex-ApproxPts1");
  
//   ring Q = RingQQ();
//   vector<RingElem> tolVector(dim, zero(Q));
//   cout << "Insert the tolerances: " << endl;
//   for (long i=0; i < dim; ++i)
//     {
//       double compt;
//       cin >> compt;
//       if (!cin || compt < 0) { CoCoA_THROW_ERROR("bad input", "main program in ex-ApproxPts1"); }
//       ConvertToRingElem(tolVector[i], compt);
//     }

//   RingElem tol(zero(Q));
//   cout << "Insert the max tolerance: " << endl;
//   double compt;
//   cin >> compt;
//   if (!cin || compt < 0) { CoCoA_THROW_ERROR("bad input", "main program in ex-ApproxPts1"); }
//   ConvertToRingElem(tol, compt);

//   cout << "Insert the number of points to be processed: " << endl;
//   long NumPts;
//   cin >> NumPts;
//   if (NumPts < 1 || NumPts > 1000000)
//     CoCoA_THROW_ERROR("Ridiculous number of points", "main program in ex-ApproxPts1");

//   vector<ApproxPts::PointR> Pts(NumPts, ApproxPts::PointR(dim, zero(Q)));
//   cout << "Insert the coordinates of the points " << endl;
//   for (long i=0; i < NumPts; ++i)
//   {
//     for (long j=0; j < dim; ++j)
//     {
//       double compt;
//       cin >> compt;
//       if (!cin) { CoCoA_THROW_ERROR("bad input", "main program in ex-ApproxPts1"); }
//       ConvertToRingElem(Pts[i][j], compt);
//     }
//   }

//   cout << endl
//        << "Read " << len(Pts) << " original points." << endl
//        << endl;

//   VanishPoly(newpts, LowDegreePoly, Pts, tolVector, tol);
// }


// //----------------------------------------------------------------------
// // Use main() to handle any uncaught exceptions and warn the user about them.
// int main()
// {
//   try
//   {
//     program();
//     return 0;
//   }
//   catch (const CoCoA::ErrorInfo& err)
//   {
//     cerr << "***ERROR***  UNCAUGHT CoCoA error";
//     ANNOUNCE(cerr, err);
//   }
//   catch (const std::exception& exc)
//   {
//     cerr << "***ERROR***  UNCAUGHT std::exception: " << exc.what() << endl;
//   }
//   catch(...)
//   {
//     cerr << "***ERROR***  UNCAUGHT UNKNOWN EXCEPTION" << endl;
//   }
//   return 1;
// }

