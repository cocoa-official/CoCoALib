//   Copyright (c)  2006,2008-2009,2017  John Abbott and Anna M. Bigatti
//   Main authors: Laura Torrente (assisted by John Abbott)

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
#include "CoCoA/MatrixOps.H"
#include "CoCoA/MatrixView.H"
#include "CoCoA/PPMonoid.H"
#include "CoCoA/QBGenerator.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingTwinFloat.H"
#include "CoCoA/SparsePolyOps-RingElem.H"
#include "CoCoA/assert.H"
#include "CoCoA/convert.H"
#include "CoCoA/VectorOps.H"
#include "CoCoA/matrix.H"
#include "CoCoA/ring.H"
#include "CoCoA/symbol.H"
#include "CoCoA/utils.H"
#include "CoCoA/verbose.H" // for VerboseLog

// #include <algorithm> // already included in apply.H
// using std::swap;   in  SOI, NBM
#include <iostream>
using std::ostream;
using std::endl;
//#include <vector>  --- already included by ApproxPts2.H
using std::vector;
#include <list>
using std::list;


namespace CoCoA
{

  namespace ApproxPts
  {

    namespace // anonymous namespace for file local auxiliary functions and definitions.
    {

      //------------------------------------------------------------------
      // Auxiliary functions for SOI

      /// Convert the input data to the RingTwinFloat with current BitPrec
      // !!! PROTECT AGAINST ALIASING!!!
      void ConversionToRR(const vector<RingElem>& OrigTolerance,
                          const vector<PointR>& OrigPts,
                          ConstRefRingElem OrigGamma,
                          vector<RingElem>& tolerance,
                          vector<PointR>& pts,
                          RingElem& gamma,
                          const ring& RR)
      {
        const ring& Q = owner(OrigTolerance.front());
        RingHom QToRR = CanonicalHom(Q,RR);
        const long dim = len(OrigTolerance);
        tolerance.clear();
        for (long i=0; i < dim; ++i)
          tolerance.push_back(QToRR(OrigTolerance[i]));
        pts.clear();
        for (const PointR& OrigPt: OrigPts)
        {
          PointR NewPt;
          for (long j=0; j < dim; ++j)
            NewPt.push_back(QToRR(OrigPt[j]));
          pts.push_back(NewPt);
        }
        gamma = QToRR(OrigGamma);
      }

      void ConversionToRR(const vector<double>& OrigTolerance,
                          const vector<PointDbl>& OrigPts,
                          double OrigGamma,
                          vector<RingElem>& tolerance,
                          vector<PointR>& pts,
                          RingElem& gamma,
                          const ring& RR)
      {
        tolerance.clear();
        for (double eps: OrigTolerance)
          tolerance.push_back(RingElem(RR, ConvertTo<BigRat>(eps)));

        pts.clear();
        for (const PointDbl& pt: OrigPts)
        {
          CoCoA_ASSERT(len(pt) == len(OrigTolerance));
          vector<RingElem> tmp;
          for (double coord: pt)
            tmp.push_back(RingElem(RR, ConvertTo<BigRat>(coord)));
          pts.push_back(tmp);
        }
        gamma = RingElem(RR, ConvertTo<BigRat>(OrigGamma));
      }


//     matrix MatrixEmbedding(ConstMatrixView M, const PolyRing& P)
//     {
//       //ring R = RingOf(M);
//       //ERROR : If (CoeffRing(P) != R)   then error!!!
//       //ERROR : If they have different size   then error!!
//       matrix NewM(NewDenseMat(P, NumRows(M), NumCols(M)));
//       RingHom phi = CoeffEmbeddingHom(P);
//       for (long i=0; i < NumRows(M); ++i)
//         for (long j=0; j < NumCols(M); ++j)
//           SetEntry(NewM, i, j, phi(M(i,j)));
//       return NewM;
//     }


//     vector<RingElem> VectorEmbedding(const vector<RingElem>& V,
//                                      const PolyRing& P)
//     {
//       //ring R = owner(V.front());
//       //ERROR : If (CoeffRing(P) != R)   then error!!!
//       RingHom phi = CoeffEmbeddingHom(P);
//       vector<RingElem> NewV;
//       for (long i=0; i < len(V); ++i)
//         NewV.push_back(phi(V[i]));
//       return NewV;
//     }


      // square 2-norm of a vector
      RingElem SquareNorm(const vector<RingElem>& v)
      {
        CoCoA_ASSERT(!v.empty());
        RingElem tmpSqNorm(zero(owner(v.front())));
        for (const RingElem& coord: v)
          tmpSqNorm += power(coord, 2);
        return tmpSqNorm;
      }


      // Evaluation of a poly f on a set Pts of points
      void Evaluation(vector<RingElem>& b,
                      ConstRefRingElem f,
                      const list<RingHom>& EvalAtPts)
      {
        b.clear();
//        for (list<RingHom>::const_iterator it = EvalAtPts.begin(); it != EvalAtPts.end(); ++it)
        for (const RingHom& phi: EvalAtPts)
          b.push_back(phi(f));
      }


      void FirstOrderEval(vector<RingElem>& V,
                          ConstRefRingElem t,
                          const list<RingHom>& EvalAtPts,
                          const vector<RingElem>& tolerance,
                          const PolyRing& PErr)
      {
        RingHom phi = CoeffEmbeddingHom(PErr);
        const long dim = len(tolerance);
        const long NumPts = len(EvalAtPts);
        vector<RingElem> tmp;
        const vector<RingElem>& e = indets(PErr);

        V.assign(NumPts, zero(PErr));
        for (long i=0; i < dim; ++i)
        {
          if (!IsZero(tolerance[i]))
          {
            Evaluation(tmp, deriv(t,i), EvalAtPts);
            for (long j=0; j < NumPts; ++j)
              V[j] += phi(tmp[j]) * e[j*dim + i];
          }
        }
      }


      // Solve the LS problem A*x = b, A orthogonal matrix
      // implemented without deleting the denominators
      void LeastSquares(vector<RingElem>& x,
                        vector<RingElem>& Rho,
                        ConstMatrixView A,
                        ConstMatrixView DiagA,
                        const vector<RingElem>& b)
      {
        const ring R = owner(b.front());
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


      // First order least squares algorithm
      void FirstOrderLS(vector<RingElem>& Rho1,
                        vector<RingElem>& x1,
                        ConstMatrixView A0,
                        ConstMatrixView A1,
                        ConstMatrixView DiagA0,
                        ConstMatrixView A0A1,
                        const vector<RingElem>& b0,
                        const vector<RingElem>& b1,
                        const vector<RingElem>& x0)
      {
        SparsePolyRing P = owner(b1.front());
        RingHom phi = CoeffEmbeddingHom(P);
        matrix A0InP = phi(A0);
        const vector<RingElem> x0InP = phi(x0);//VectorEmbedding(x0, P);

        // x1 = A0*b1;
        x1.assign(NumRows(A0), zero(P));
        A0InP->myMulByCol(x1,b1);
        // B = A1*b0;
        vector<RingElem> B(NumRows(A1), zero(P));
        A1->myMulByCol(B, phi(b0));
        // C = A0A1*x0;
        vector<RingElem> C(NumRows(A0A1), zero(P));
        A0A1->myMulByCol(C,x0InP);
        for (long i=0; i < len(b1); ++i)
        {
          if (IsZero(DiagA0(i,i)))
            x1[i] = zero(P);
          else
          {
            x1[i] += (B[i] - C[i]);
            x1[i] *= phi(1/DiagA0(i,i));
          }
        }
        Rho1 = b1;
        // D = A0*x1;
        vector<RingElem> D(NumRows(A0A1), zero(P));
        A0InP->myMulByRow(D,x1);
        // E =  A1*x0;
        vector<RingElem> E(NumRows(A1), zero(P));
        A1->myMulByRow(E,x0InP);
        for (long i=0; i < len(b1); ++i)
          Rho1[i] -= (D[i] + E[i]);
      }


      // Compute the matrix C that represents the vector V
      void ErrorMatrix(MatrixView& C, const vector<RingElem>& V)
      {
        SparsePolyRing P = owner(V.front());
        const long nvars = NumIndets(P);
        RingHom PToR = EvalHom(P, vector<RingElem>(nvars, zero(CoeffRing(P))));
        for (long i=0; i < len(V); ++i)
          for (long j=0; j < nvars; ++j)
            SetEntry(C,i,j,PToR(deriv(V[i],j)));
      }


// //Compute an estimate of the squared minimal singular value of a matrix
// void LowerBoundMinSingValue(RingElem& LowerBound,
//                             ConstMatrixView A)
// {
//   list<RingElem> listOfNorms;
//   RingElem norm(RingOf(A));
//   inftyMatrixNorm(norm, A);
//   listOfNorms.push_back(1/(NumRows(A)*norm*norm));
//   GlobalOutput() << "LB1 = " << 1/(NumRows(A)*norm*norm) << endl;
//   inftyMatrixNorm(norm, transpose(A));
//   listOfNorms.push_back(1/(NumCols(A)*norm*norm));
//   GlobalOutput() << "LB2 = " << 1/(NumCols(A)*norm*norm) << endl;
//   FrobeniusNorm(norm, A);
//   listOfNorms.push_back(1/norm);
//   GlobalOutput() << "LB3 = " << 1/norm << endl;
//   listOfNorms.sort();
//   list<RingElem>::iterator it = listOfNorms.end();
//   LowerBound = *(--it);
//   GlobalOutput() << "Lower Bound is = " << LowerBound << endl;
// }



      // Compute the vector of rows of the matrix A such that the minimal
      // singular value of the corresp submatrix is greater than ||Eps||_2
      void WellIndepRows(vector<long>& rows,
                         ConstRefRingElem Eps_square,
                         ConstMatrixView A)
      {
        const ring R = RingOf(A);
        vector<long> cols;
        for (long j=0; j < NumCols(A); ++j)  cols.push_back(j);
        // cols = LongRange(0, NumCols(A)-1);  // AllCols
        rows.clear();
        for (long i=0; i < NumRows(A); ++i)
        {
          rows.push_back(i);
          ConstMatrixView B = submat(A, rows, cols);
          // Check lower bound for smallest singular value of B
          if (rk(B) < NumRows(B) || 1/FrobeniusNormSq(PseudoInverse(B)) < Eps_square)
          {
            rows.pop_back();
          }
        }
      }



      // Compute minimal 2-norm solution of a Full Rank Underdetermined System
      void SolveFRUDS(vector<RingElem>& x,
                      ConstMatrixView A,
                      const vector<RingElem>& b)
      {
        ConstMatrixView PseudoInvA = PseudoInverse(A);
        // x = PseudoInvA*b;
        x.assign(NumRows(PseudoInvA), zero(RingOf(A)));
        PseudoInvA->myMulByCol(x, b);
      }



      // Add to QB and update
      void AddToQBAndUpdate(const PPMonoidElem& t,
                            QBGenerator& QBG,
                            const vector<RingElem>& Alpha0,
                            vector<RingElem>& S)
      {
        VerboseLog VERBOSE("AddToQBAndUpdate");
        SparsePolyRing P = owner(S.front());
        RingHom RToP = CoeffEmbeddingHom(P);
        const RingElem tAsPoly = monomial(P, t);

        VERBOSE(90) << t << " added to QB" << endl;
        const long lenQB = len(QBG.myQB());
        S.push_back(tAsPoly);
        for (long i=0; i < lenQB; ++i)
          S[lenQB] -= RToP(Alpha0[i])*S[i];
        QBG.myCornerPPIntoQB(t);
      }


      // Update evaluation matrices
      void UpdateMatrices(const QBGenerator& QBG,
                          MatrixView& A0,
                          MatrixView& A1,
                          MatrixView& DiagA0,
                          MatrixView& A0A1,
                          const vector<RingElem>& Rho0,
                          const vector<RingElem>& Rho1)
      {
        RingHom coeff = CoeffEmbeddingHom(RingOf(A1));
        const long k = len(QBG.myQB());
        SetEntry(DiagA0, k, k, SquareNorm(Rho0));
        for (long j=0; j < len(Rho0); ++j)
        {
          SetEntry(A0, k, j, Rho0[j]);
          SetEntry(A1, k, j, Rho1[j]);
          SetEntry(A0A1, k, k, A0A1(k, k) + 2*Rho1[j]*coeff(Rho0[j]));
        }
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
                                const vector<RingElem>& S,
                                vector<RingElem>& J)
      {
        VerboseLog VERBOSE("AddToCornerAndUpdate");
        SparsePolyRing P = owner(S.front());
        RingHom RToP = CoeffEmbeddingHom(P);
        const RingElem tAsPoly = monomial(P, t);

        VERBOSE(90) << t << " added to Corner" << endl;
        RingElem f(tAsPoly);
        for (long i=0; i < len(S); ++i)
          f -= RToP(Alpha0[i])*S[i];
        J.push_back(f);
        QBG.myCornerPPIntoAvoidSet(t);
      }



      // Compute the (first order) border of an order ideal
      void ComputeBorderPPs(list<PPMonoidElem>& BorderPPs,
                            const QBGenerator& QBG)
      {
        const vector<PPMonoidElem> QB = QBG.myQB();
        BorderPPs.clear();
        if (QB.empty()) // can this ever happen???
        {
          BorderPPs.push_back(one(PPM(QBG)));
          return;
        }
        const vector<PPMonoidElem>& X = indets(PPM(QBG));
        const long NumIndets = len(X);
        for (long i=0; i < len(QB); ++i)
          for (long j=0; j < NumIndets; ++j)
            BorderPPs.push_back(QB[i]*X[j]);
        BorderPPs.sort();
        BorderPPs.unique();
        for (long i=0; i < len(QB); ++i)
          BorderPPs.remove(QB[i]);
      }


      // Compute the poly f = t-sum x_1 S_i
      void ComputeBorderPoly(RingElem& f,
                             const PPMonoidElem& t,
                             ConstMatrixView A,
                             ConstMatrixView DiagA,
                             const list<RingHom>& EvalAtPts,
                             const vector<RingElem>& S)
      {
        SparsePolyRing P = owner(S.front());
        RingHom phi = CoeffEmbeddingHom(P);
        f = monomial(P, t);
        vector<RingElem> b, x, Rho;
        Evaluation(b, f, EvalAtPts);
        LeastSquares(x, Rho, A, DiagA, b);
        for (long i=0; i < len(S); ++i)
          f -= phi(x[i])*S[i];
      }



      //------------------------------------------------------------------
      // Auxiliary functions for NBM

      // Build the matrix of absolute values
      matrix AbsMatrix(ConstMatrixView M)
      {
        matrix AbsM = NewDenseMat(RingOf(M), NumRows(M), NumCols(M));
        for (long i=0; i < NumRows(M); ++i)
          for (long j=0; j < NumCols(M); ++j)
            SetEntry(AbsM, i, j, abs(M(i,j)));
        return AbsM;
      }


      // Build the vector of absolute values
      vector<RingElem> AbsVector(const vector<RingElem>& V)
      {
        vector<RingElem> AbsV; AbsV.reserve(len(V));
        for (const RingElem& coord: V)
          AbsV.push_back(abs(coord));
        return AbsV;
      }


      // Computing the bound for the sufficient condition...
      void NBMBound(vector<RingElem>& Bound,
                    ConstMatrixView M,
                    const vector<RingElem>& tolerance,
                    const RingElem& g,
                    const list<RingHom>& EvalAtPts)
      {
        const long NumPts = len(EvalAtPts);
        const ring R = RingOf(M);
        RingHom phi = CanonicalHom(R, owner(g));
        matrix id = NewDenseMat(IdentityMat(R, NumPts));
        matrix H = transpose(M)*PseudoInverse(transpose(M));
        matrix C = NewDenseMat(IdentityMat(R, NumPts));
        for (long i=0; i < NumPts; ++i)
          for (long j=0; j < NumPts; ++j)
            SetEntry(C, i, j, id(i,j)-H(i,j));
        matrix CAbs = AbsMatrix(C);

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


      void ApproxPts_CheckInput(const SparsePolyRing& P,
                                const std::vector<PointR>& pts,
                                const std::vector<RingElem>& tolerance)
      {
        const ring K = CoeffRing(P); // K must be an ordered field
        CoCoA_ASSERT(IsOrderedDomain(K) && IsField(K));
        if (owner(pts[0][0]) != K)
          CoCoA_THROW_ERROR1(ERR::MixedRings);
        const long dim = len(tolerance);
        if (NumIndets(P) < dim)
          CoCoA_THROW_ERROR1(ERR::BadNumIndets);
#ifdef CoCoA_DEBUG
        for (const RingElem& eps: tolerance)
          if (eps < 0)
            CoCoA_THROW_ERROR1(ERR::ReqNonNegative);
#endif
      }
    

    } // end of anonymous namespace

//     void SOI(std::vector<PPMonoidElem>& QB,
//              std::vector<RingElem>& BBasis,
//              std::vector<RingElem>& AlmostVanishing,
//              const std::vector<PointR>& pts,
//              const std::vector<RingElem>& tolerance,
//              ConstRefRingElem gamma)

    void SOI(std::vector<PPMonoidElem>& QB,
             std::vector<RingElem>& BBasis,
             std::vector<RingElem>& AlmostVanishing,
             const SparsePolyRing& P,
             const std::vector<PointR>& pts,
             const std::vector<RingElem>& tolerance,
             ConstRefRingElem gamma)
    {
      VerboseLog VERBOSE("ApproxPts::SOI");
      ApproxPts_CheckInput(P, pts, tolerance);
      const long dim = len(tolerance);
      const long NumPts = len(pts);

      ring K = CoeffRing(P);
      PPMonoid PPMon = PPM(P);
      RingHom coeff = CoeffEmbeddingHom(P);

// P given as input AMB 2017-03
//       // Build poly rings P=QQ[x_1...x_n] and
//       // PErr = QQ[e[1,1],...,e[NumPts,dim]]
//       PPMonoid PPM = NewPPMonoid(SymbolRange("x", 0, dim-1), StdDegRevLex);
//       SparsePolyRing P = NewPolyRing(R, PPM);
      // PErr = K[e[1,1],...,e[NumPts,dim]]
      SparsePolyRing PErr = NewPolyRing(K, SymbolRange(symbol("e",1,1), symbol("e",NumPts,dim)));
      // Build evaluation homomorphisms at Pts
      list<RingHom> EvalAtPts;///  EvalAtPts.reserve(NumPts);
      for (long i=0; i < NumPts; ++i)
        EvalAtPts.push_back(EvalHom(P, pts[i]));

      //   --------------     INITIALIZE   ---------------
      long IterNum = 1;
      const RingElem Eps_square = SquareNorm(tolerance);
      const RingElem kappa = gamma+1;
      const RingElem Bound = kappa*kappa*NumPts*Eps_square;
      vector<RingElem> J, S;     //J = almost vanishing polys

      QBGenerator QBG(PPMon);      //to handle quotient bases
      QBG.myCornerPPIntoQB(one(PPMon));
      S.push_back(one(P));
      list<PPMonoidElem> L = QBG.myCorners();  //L = list of PP

      // Build and initialize matrices A0, A1 DiagA0, A0A1, CErr
      matrix A0(NewDenseMat(K, NumPts, NumPts));   // A0 = orthogonal
      for (long i=0; i < NumPts; ++i)
        SetEntry(A0, 0, i, 1);

      matrix A1(NewDenseMat(PErr, NumPts, NumPts));

      vector<RingElem> DiagA0Elems(NumPts, zero(K));
      MatrixView DiagA0 = DiagMat(DiagA0Elems); // DiagA0 = A0*A0'
      SetEntry(DiagA0, 0, 0, NumPts);

      vector<RingElem> A0A1Elems(NumPts, zero(PErr));
      MatrixView A0A1 = DiagMat(A0A1Elems); // A0A1 = 2*A0*A1'


      // ------------   MAIN CYCLE   -----------------
      while (!L.empty())
      {
        const PPMonoidElem t = L.front(); // PP considered in the current step
        const RingElem tAsPoly = monomial(P, t);
        VERBOSE(95) << "---------------------------------------" << endl;
        VERBOSE(80) << "Iteration " << IterNum << " PP = " << t << endl;
//???        VERBOSE(95) << "Quotient Basis = " << QBG.myQB() << endl;
//???        VERBOSE(95) << "List of power products = " << L << endl;

        // Compute eval V0=t(Pts) and solve LS problem A0'*b = V0
        vector<RingElem> V0, V1, Alpha0, Alpha1, Rho0, Rho1;
        Evaluation(V0, tAsPoly, EvalAtPts);
        LeastSquares(Alpha0, Rho0, A0, DiagA0, V0);
        // First trivial check
        if (IsZero(Rho0))
        {
          AddToCornerAndUpdate(t, QBG, Alpha0, S, J);
          VERBOSE(90) << "Truly vanishing polynomial " /*<< J.back()*/ << endl;
          ++IterNum;
          L = QBG.myCorners();
          continue;
        }

        // Compute f.o.eval V1, solve (f.o.error anal.)LS problem
        FirstOrderEval(V1, tAsPoly, EvalAtPts, tolerance, PErr);
        FirstOrderLS(Rho1, Alpha1, A0, A1, DiagA0, A0A1, V0, V1, Alpha0);
        // Compute CErr such that Rho1 = CErr*e
        matrix CErr(NewDenseMat(K, NumPts, NumPts*dim));//CErr=Rho1(matrix)
        ErrorMatrix(CErr, Rho1);
        // Compute the strongly independent rows of CErr
        vector<long> rows;
        WellIndepRows(rows, Eps_square, CErr);
        VERBOSE(80) << "Well Independent Rows =" << rows << endl;

        // if, for each i, ||CErr[i]|| < Eps ...
        if (rows.empty())
        {
          if (SquareNorm(Rho0) > Bound*Eps_square)
          {
            UpdateMatrices(QBG, A0, A1, DiagA0, A0A1,Rho0,Rho1);
            AddToQBAndUpdate(t, QBG, Alpha0, S);
          }
          else
          {
            AddToCornerAndUpdate(t, QBG, Alpha0, S, J);
            VERBOSE(90) << "Almost vanishing polynomial " /*<< J.back()*/ << endl;
          }
        }
        else
        {
          // Compute the submatrix of CErr with rows in rows
          vector<long> cols;
          for (long j=0; j < NumCols(CErr); ++j)  cols.push_back(j);
          // cols = LongRange(0, NumCols(CErr)-1);  // AllCols
          ConstMatrixView CErr1 = submat(CErr, rows, cols);
          vector<RingElem> MinusRho01;
          for (long i=0; i < len(rows); ++i)
            MinusRho01.push_back(-Rho0[rows[i]]);
          vector<RingElem> eLS1;
          // Compute the minimal 2-norm solution of the
          // underdetermined system CErr1 * eLS1 = MinusRho01
          SolveFRUDS(eLS1, CErr1, MinusRho01);
          VERBOSE(95) << "Square norm eLS1= " << SquareNorm(eLS1) << endl
                      << "Bound = " << Bound << endl;

          //Check and update the variables
          if (SquareNorm(eLS1) >= Bound)
          {
            UpdateMatrices(QBG, A0, A1, DiagA0, A0A1,Rho0,Rho1);
            AddToQBAndUpdate(t, QBG, Alpha0, S);
          }
          else
          {
            AddToCornerAndUpdate(t, QBG, Alpha0, S, J);
            VERBOSE(90) << "Almost vanishing polynomial " /*<< J.back()*/ << endl;
          }
        }
        ++IterNum;
        L = QBG.myCorners();
      }

      VERBOSE(95) << "==========================================" << endl
                  << "All power-products have been considered" << endl;

      // Computation of the Border Basis and printing
      BBasis.clear();
      std::swap(AlmostVanishing, J); // morally AlmostVanishing = J;
      QB.clear(); QB = QBG.myQB(); // cannot use swap because myQB is read only
      if (len(QB) < NumPts)
      {
        VERBOSE(80) << "len(QB) = " << len(QBG.myQB())
                    << " for " << len(pts) << " points:"
                    << " finer tolerance needed" << endl;
        return;
      }

      list<PPMonoidElem> BorderPPs;
      ComputeBorderPPs(BorderPPs, QBG);
//      for (list<PPMonoidElem>::const_iterator it=BorderPPs.begin(); it != BorderPPs.end(); ++it)
      for (auto&& t: BorderPPs)
      {
        RingElem BBElem(P);
        ComputeBorderPoly(BBElem, t, A0, DiagA0, EvalAtPts, S);
        BBasis.push_back(zero(P)); swap(BBasis.back(), BBElem); // morally  BBasis.push_back(f);
      }
    }


    void SOITwinFloat(std::vector<PPMonoidElem>& QB,
                      std::vector<RingElem>& BBasis,
                      std::vector<RingElem>& AlmostVanishing,
                      const SparsePolyRing& P,
                      const std::vector<PointR>& OrigPts,
                      const std::vector<RingElem>& OrigTolerance,
                      ConstRefRingElem OrigGamma)
    {
      VerboseLog VERBOSE("ApproxPts::SOITF");
      ApproxPts_CheckInput(P, OrigPts, OrigTolerance);
      const long dim = len(OrigTolerance);
///      const long NumPts = len(OrigPts);

      PPMonoid PPM = NewPPMonoid(SymbolRange("x",0,dim-1), StdDegRevLex);
      for (long BitPrec = 128; ; BitPrec *= 2)
      {
        try
        {
          const ring RR = NewRingTwinFloat(BitPrec);
          vector<RingElem> tolerance;
          vector<PointR> pts;
          RingElem gamma(RR);
          ConversionToRR(OrigTolerance, OrigPts, OrigGamma, tolerance, pts, gamma, RR);
          const SparsePolyRing P = NewPolyRing(RR, PPM);
          SOI(QB, BBasis, AlmostVanishing, P, pts, tolerance, gamma);
          VERBOSE(80) << "BBasis successfully computed with " << BitPrec << " bit precision " << endl;
          return;
        }
        catch (const RingTwinFloat::InsufficientPrecision&)
        {
          // If verbose mode, inform the user about the failure...
          VERBOSE(80) << endl << "----------------------------------"
                      << endl << "A bit precision of " << BitPrec
                      << " was not sufficient." << endl
                      << "------------------------------------------"
                      << endl << endl;
        }
      }
    }


    void SOI(std::vector<PPMonoidElem>& QB,
             std::vector<RingElem>& BBasis,
             std::vector<RingElem>& AlmostVanishing,
             const SparsePolyRing& P,
             const std::vector<PointDbl>& OrigPts,
             const std::vector<double>& OrigTolerance,
             double OrigGamma)
    {
      const ring QQ = RingQQ();
      vector<RingElem> tolerance;
      vector<PointR> pts;
      RingElem gamma(QQ);
      ConversionToRR(OrigTolerance, OrigPts, OrigGamma, tolerance, pts, gamma, QQ);
      SOITwinFloat(QB, BBasis, AlmostVanishing, P, pts, tolerance, gamma);
    }




//--------------  NBM  ---------------------------------------------

    void NBM(std::vector<PPMonoidElem>& QB,
             std::vector<RingElem>& BBasis,
             std::vector<RingElem>& AlmostVanishing,
             const std::vector<PointR>& pts,
             const std::vector<RingElem>& tolerance)
    {
      const ring R = owner(pts[0][0]); // R must be an ordered field
      CoCoA_ASSERT(IsOrderedDomain(R) && IsField(R)); //
      const long dim = len(tolerance);

      SparsePolyRing P = NewPolyRing(R, NewSymbols(dim)); // Build poly ring with dim indets
      NBM(QB, BBasis, AlmostVanishing, P, pts, tolerance);
    }

    void NBM(std::vector<PPMonoidElem>& QB,
             std::vector<RingElem>& BBasis,
             std::vector<RingElem>& AlmostVanishing,
             const SparsePolyRing& P,
             const std::vector<PointR>& pts,
             const std::vector<RingElem>& tolerance)
    {
      VerboseLog VERBOSE("ApproxPts::NBM");
      ApproxPts_CheckInput(P, pts, tolerance);
      const long NumPts = len(pts);

      ring K = CoeffRing(P);
      PPMonoid PPMon = PPM(P);
      RingHom coeff = CoeffEmbeddingHom(P);

      // Build evaluation homomorphisms at Pts
      list<RingHom> EvalAtPts;
      for (long i=0; i < NumPts; ++i)
        EvalAtPts.push_back(EvalHom(P, pts[i]));

      //   --------------     INITIALIZE   ---------------
      long IterNum = 1;
      vector<RingElem> J, S;     // J = almost vanishing polys

      QBGenerator QBG(PPMon);      // to handle quotient bases
      QBG.myCornerPPIntoQB(one(PPMon));
      S.push_back(one(P));
      list<PPMonoidElem> L = QBG.myCorners();  // L = list of corner PPs

      // Build and initialize evaluation matrix M (by rows),
      // its orthogonal part A, and DiagA = A*A'
      matrix A(NewDenseMat(K, NumPts, NumPts));
      for (long i=0; i < NumPts; ++i)
        SetEntry(A, 0, i, 1);

      matrix M(NewDenseMat(K, NumPts, NumPts));
      for (long i=0; i < NumPts; ++i)
        SetEntry(M, 0, i, 1);

      vector<RingElem> DiagAElems(NumPts, zero(K));
      MatrixView DiagA(DiagMat(DiagAElems));//DiagA=A*A'
      SetEntry(DiagA, 0, 0, NumPts);


      // ------------   MAIN CYCLE   ----------------- verb 80/90/95
      while (!L.empty())
      {
        const PPMonoidElem t = L.front(); // PP considered in the current step
        const RingElem tAsPoly = monomial(P, t);
        VERBOSE(95) << "---------------------------------------" << endl;
        VERBOSE(80) << "Iteration " << IterNum << " PP = " << t << endl;
//???        VERBOSE(95) << "Quotient Basis = " << QBG.myQB() << endl;
//???        VERBOSE(95) << "List of power products = " << L << endl;

        // Compute eval V=t(Pts) and solve LS problem A'*b = V
        vector<RingElem> V, Alpha, Rho;
        Evaluation(V, tAsPoly, EvalAtPts);
        LeastSquares(Alpha, Rho, A, DiagA, V);

        // First trivial check
        if (IsZero(Rho))
        {
          AddToCornerAndUpdate(t, QBG, Alpha, S, J);
          VERBOSE(90) << "Truly vanishing polynomial " /*<< J.back()*/ << endl;
          ++IterNum;
          L = QBG.myCorners();
          continue;
        }

        // build the poly g = t - sum(alpha_i s_i)
        RingElem g(tAsPoly);
        for (long i=0; i < len(S); ++i)
          g -= coeff(Alpha[i])*S[i];

        vector<long> cols, rows;
        for (long i=0; i < len(QBG.myQB()); ++i)
          rows.push_back(i);
        // rows = LongRange(0, len(QBG.myQB())-1);
        for (long j=0; j < NumCols(M); ++j)  cols.push_back(j);
        // cols = LongRange(0, NumCols(M)-1);  // AllCols
        ConstMatrixView SubM = submat(M, rows, cols);

        // computing the bound for NBM
        vector<RingElem> Bound;
        NBMBound(Bound, SubM, tolerance, g, EvalAtPts);

        RingElem Diff(K);
        vector<RingElem> AbsRho = AbsVector(Rho);
        for (long i=0; i < NumPts; ++i)
        {
          Diff = AbsRho[i]-Bound[i];
          if (Diff > 0)
            break;
        }

        // if Diff = |Rho[i]|-Bound[i]>0 then ...
        if (Diff > 0)
        {
          UpdateMatricesNBM(QBG, M, A, DiagA, V, Rho);
          AddToQBAndUpdate(t, QBG, Alpha, S);
        }
        else
        {
          AddToCornerAndUpdate(t, QBG, Alpha, S, J);
          VERBOSE(90) << "Almost vanishing polynomial " /*<< J.back()*/ << endl;
        }

        ++IterNum;
        L = QBG.myCorners();
      }
      VERBOSE(95) << "==========================================" << endl
                  << "All power-products have been considered" << endl;
      // Computation of the Border Basis and printing
      BBasis.clear();
      std::swap(AlmostVanishing, J); // morally AlmostVanishing = J;
      QB.clear(); QB = QBG.myQB();
      if (len(QB) < NumPts)
      {
        VERBOSE(80) << "len(QB) = " << len(QBG.myQB())
                    << " for " << len(pts) << " points:"
                    << " finer tolerance needed" << endl;
        return;
      }

      list<PPMonoidElem> BorderPPs;
      ComputeBorderPPs(BorderPPs, QBG);
//      for (list<PPMonoidElem>::const_iterator it=BorderPPs.begin(); it != BorderPPs.end(); ++it)
      for (auto&& t: BorderPPs)
      {
        RingElem BBElem(P);
        ComputeBorderPoly(BBElem, t, A, DiagA, EvalAtPts, S);
        BBasis.push_back(zero(P)); swap(BBasis.back(), BBElem); // morally  BBasis.push_back(BBElem);
      }
    }



    void NBMTwinFloat(std::vector<PPMonoidElem>& QB,
                      std::vector<RingElem>& BBasis,
                      std::vector<RingElem>& AlmostVanishing,
                      const std::vector<PointR>& OrigPts,
                      const std::vector<RingElem>& OrigTolerance)
    {
      VerboseLog VERBOSE("NBMTwinFloat");
      for (long BitPrec = 128; ; BitPrec *= 2)
      {
        try
        {
          const ring RR = NewRingTwinFloat(BitPrec);
          vector<RingElem> tolerance;
          vector<PointR> pts;
          RingElem gamma(zero(RR));
          ConversionToRR(OrigTolerance, OrigPts, zero(owner(OrigPts[0][0])), tolerance, pts, gamma, RR);
          NBM(QB, BBasis, AlmostVanishing, pts, tolerance);
          VERBOSE(80) << "BBasis successfully computed with " << BitPrec << " bit precision " << endl;
          return;
        }
        catch (const RingTwinFloat::InsufficientPrecision&)
        {
          // Inform the user about the failure...
          VERBOSE(80) << endl << "----------------------------------"
                      << endl << "A bit precision of " << BitPrec
                      << " was not sufficient." << endl
                      << "------------------------------------------"
                      << endl << endl;
        }
      }
    }


    void NBM(std::vector<PPMonoidElem>& QB,
             std::vector<RingElem>& BBasis,
             std::vector<RingElem>& AlmostVanishing,
             const std::vector<PointDbl>& OrigPts,
             const std::vector<double>& OrigTolerance)
    {
      const ring QQ = RingQQ();
      vector<RingElem> tolerance;
      vector<PointR> pts;
      RingElem gamma(QQ);
      ConversionToRR(OrigTolerance, OrigPts, 0.0, tolerance, pts, gamma, QQ);
      NBMTwinFloat(QB, BBasis, AlmostVanishing, pts, tolerance);
    }


//     void NBM(std::vector<PPMonoidElem>& QB,
//              std::vector<RingElem>& BBasis,
//              std::vector<RingElem>& AlmostVanishing,
//              const SparsePolyRing& P,
//              const std::vector<PointDbl>& OrigPts,
//              const std::vector<double>& OrigTolerance)
//     {
//       ring QQ = RingQQ();
//       vector<RingElem> tolerance;
//       vector<PointR> pts;
//       RingElem gamma(QQ);
//       ConversionToRR(OrigTolerance, OrigPts, 0.0, tolerance, pts, gamma, QQ);
//       NBMTwinFloat(QB, BBasis, AlmostVanishing, P, pts, tolerance, gamma);
//     }


  } // end of namespace ApproxPts

} // end of namespace CoCoA
