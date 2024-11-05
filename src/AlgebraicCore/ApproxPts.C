//   Copyright (c)  2006-2009,2013  John Abbott and Anna M. Bigatti
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


#include "CoCoA/ApproxPts.H"

#include "CoCoA/BigRat.H"
#include "CoCoA/RingTwinFloat.H"
#include "CoCoA/assert.H"
#include "CoCoA/convert.H"
#include "CoCoA/error.H"
#include "CoCoA/ring.H"
#include "CoCoA/utils.H"

#include <algorithm>
using std::find_if;
#include <cmath>
// using std::log; // only in UnitVolBallRadius !! Shadowed by the CoCoALib fn log(BigInt) !!
// using std::exp; // only in UnitVolBallRadius
#include <list>
using std::list;
//#include <vector>
using std::vector;


namespace CoCoA
{

  typedef vector<double> PointDbl;

  namespace // anonymous namespace for file local auxiliary functions and definitions.
  {

    // The square of x -- useful if x is a big expression.
    inline double square(double x)
    {
      return x*x;
    }
    

    // This is not really specific to approximate points, but where else should it go?
    // Compute radius of D-ball having volume 1.
    // (I got the formula from Wikipedia under "hypersphere")
    double UnitVolBallRadius(long D)
    {
      using std::log;
      using std::exp;
      CoCoA_ASSERT(D < 10000); // disallow very high dimensions
      const double pi = 4*atan(1.0);
      double logfactor = (D/2)*(-log(2*pi)); // D/2 exploits integer division!  NB -(D/2)*log(...) does NOT WORK.
      if (D & 1)
        logfactor -= log(2.0);
      for (long k=D; k > 1; k -= 2)
        logfactor += log(double(k));
      return exp(logfactor/D);
    }


    // Compute square of 2-norm between P and Q.
    double L2NormSq(const PointDbl& P, const PointDbl& Q)
    {
      CoCoA_ASSERT(len(P) == len(Q));
      double SumSq = 0;
      const long dim = len(P);
      for (long i=0; i < dim; ++i)
        SumSq += square(P[i]-Q[i]);
      return SumSq;
    }


    // Data structure for representing PreProcessed Points.
    struct PPP
    {
    public: // data members
      PointDbl myCoords;            // the coords of the preprocessed point -- average of myOrigPts
      long myNumOrigPts;            // length of list myOrigPts (used only in PreprocessSubdivAlgm)
      list<long> myOrigPts;         // indices into a master table
      vector<double> myDistanceSq;  // used only in PreprocessSubdivAlgm
      long myWorstOrigPt;         // the index of the farthest point
    };


    // Returns the average of a set of PointDbls; the set is indicated by the first arg
    // which is a list of indices into the second arg.
    PointDbl SubsetAve(const list<long>& IndexList, const vector<PointDbl>& pts)
    {
      CoCoA_ASSERT(!pts.empty());
      const long dim = len(pts[0]);
      PointDbl CofG(dim);
      for (long i: IndexList)
      {
        const PointDbl& P = pts[i];
        for (long j=0; j < dim; ++j)
          CofG[j] += P[j];
      }

      const long n = len(IndexList);
      for (long j=0; j < dim; ++j)
        CofG[j] /= n;
      return CofG;
    }


    ApproxPts::PointR SubsetAve(const list<long>& IndexList, const vector<ApproxPts::PointR>& pts)
    {
      CoCoA_ASSERT(!pts.empty());
      const long dim = len(pts[0]);
      ApproxPts::PointR ave(dim, zero(owner(pts[0][0])));
      for (long i: IndexList)
      {
        const ApproxPts::PointR& P = pts[i];
        for (long j=0; j < dim; ++j)
          ave[j] += P[j];
      }

      const long n = len(IndexList);
      for (long j=0; j < dim; ++j)
        ave[j] /= n;
      return ave;
    }


    // Choose a suitable positive value for any tolerance component which is zero,
    // and replace that component by the chosen positive value.
    void ZeroToleranceHack(const vector<ApproxPts::PointR>& OrigPts, vector<RingElem>& tolerance)
    {
      const long dim = len(tolerance);
      const long NumPts = len(OrigPts);
      for (long i=0; i < dim; ++i)
      {
        if (tolerance[i] != 0) continue;

        RingElem MinDist = tolerance[i]; // really just assignment to zero.
        for (long j=0; j < NumPts; ++j)
          for (long k=j+1; k < NumPts; ++k)
          {
            const RingElem d = abs(OrigPts[j][i]-OrigPts[k][i]);
            if (MinDist == 0 || d < MinDist)
              MinDist = d;
          }
        if (MinDist != 0)
          tolerance[i] = MinDist/4;
        else
          tolerance[i] = 1;
      }
    }


    // Rescale coords of pt by factors in array ScaleFactor.
    void ScalePt(PointDbl& ScaledPt, const ApproxPts::PointR& OrigPt, const vector<RingElem>& ScaleFactor)
    {
      CoCoA_ASSERT(len(OrigPt) > 0);
      CoCoA_ASSERT(len(OrigPt) == len(ScaleFactor));
      const long dim = len(OrigPt);
      ScaledPt.resize(dim);
      BigRat Q;
      RingElem coord(owner(OrigPt[0]));
      for (long i=0; i < dim; ++i)
      {
        coord = OrigPt[i] * ScaleFactor[i];
        if (!IsRational(Q, coord) || !IsConvertible(ScaledPt[i], Q))
          CoCoA_THROW_ERROR("Cannot convert to double", "ScalePt");
      }
    }

    // Rescale the OrigPts and convert them into PointDbl
    void NormalizePts(vector<PointDbl>& ScaledPts, const vector<ApproxPts::PointR>& OrigPts, const vector<RingElem>& tolerance)
    {
      const long NumPts = len(OrigPts);
      const long dim = len(tolerance);

      ring R = owner(tolerance[0]);
      vector<RingElem> InverseTolerance(dim, zero(R));
      for (long j=0; j < dim; ++j)
        InverseTolerance[j] = 1/tolerance[j];

      ScaledPts.resize(NumPts);
      for (long i=0; i < NumPts; ++i)
        ScalePt(ScaledPts[i], OrigPts[i], InverseTolerance);
    }


    /////////////////////////////////////////////////////////////////////////////
    // Auxiliary definitions for the grid algorithm.

    PointDbl GridNearPoint(const PointDbl& P)
    {
      const long dim = len(P);
      PointDbl GNP(dim);
      for (long i=0; i < dim; ++i)
        GNP[i] = 2*round(P[i]/2);
      return GNP;
    }


    // Only used in GridAlgm as arg to std::find.
    // WARNING data member is a reference!!
    class CoordsEq
    {
    public:
      CoordsEq(const PointDbl& P): myPointDbl(P) {};
      bool operator()(const PPP& PPPt) { return myPointDbl == PPPt.myCoords; };
    private: // data members
      const PointDbl& myPointDbl;
    };


    /////////////////////////////////////////////////////////////////////////////
    // Functions for aggregative algorithm

    // ave = (w1*P1 + w2*P2)/(w1+w2);
    void WeightedAve(PointDbl& ave, double w1, const PointDbl& P1, double w2, const PointDbl& P2)
    {
      const long dim = len(ave);
      CoCoA_ASSERT(len(P1) == dim && len(P2) == dim);
      CoCoA_ASSERT(w1+w2 != 0);
      for (long i=0; i < dim; ++i)
      {
        ave[i] = (w1*P1[i]+w2*P2[i])/(w1+w2);
      }
    }


    struct NearPair
    {
    public:
      NearPair(long i, long j, double sepsq): myLowerIndex(i), myHigherIndex(j), mySepSq(sepsq) {};
    public: // data members
      const long myLowerIndex, myHigherIndex;
      double mySepSq;
    };


    // Used only as arg to std::sort algorithm
    inline bool LessNearPair(const NearPair& P1, const NearPair& P2)
    {
      return P1.mySepSq < P2.mySepSq;
    }


    // Used in call to algorithm std::remove_if.
    class NearPairContains
    {
    public:
      NearPairContains(long j): myIndex(j) {}
      bool operator()(const NearPair& P) { return P.myLowerIndex == myIndex || P.myHigherIndex == myIndex; }
    private: // data member
      const long myIndex;
    };




    /////////////////////////////////////////////////////////////////////////////
    // Below is auxiliary code for the "subdiv" algorithm.

    // Fill in data member myDistanceSq of X
    void ComputeDistSq(PPP& X, const vector<PointDbl>& OrigPts)
    {
      for (long i=0; i < len(OrigPts); ++i)
        X.myDistanceSq[i] = L2NormSq(X.myCoords, OrigPts[i]);

      double WorstDist = 0;
      for (long i: X.myOrigPts)
        if (X.myDistanceSq[i] >= WorstDist)
        {
          WorstDist = X.myDistanceSq[i];
          X.myWorstOrigPt = i;
        }
    }


    // Make Y a copy of X but without the i-th OrigPt.
    void EraseAndUpdate(PPP& Y, const PPP& X, long i, const vector<PointDbl>& V)
    {
      Y.myOrigPts = X.myOrigPts;
      Y.myOrigPts.remove(i);
      Y.myNumOrigPts = X.myNumOrigPts-1;
      Y.myCoords = SubsetAve(Y.myOrigPts, V);
      ComputeDistSq(Y, V);
    }


    // Make Y a copy of X with the i-th OrigPt added in.
    void AddAndUpdate(PPP& Y, const PPP& X, long i, const vector<PointDbl>& V)
    {
      Y.myOrigPts = X.myOrigPts;
      Y.myOrigPts.push_back(i);
      Y.myNumOrigPts = X.myNumOrigPts+1;
      Y.myCoords = SubsetAve(Y.myOrigPts, V);
      ComputeDistSq(Y, V);
    }


    // Redistribute the assignments of OrigPts to PPPts to minimize the sum of squares of distances
    void redistribute(vector<PPP>& PPPts, const vector<PointDbl>& OrigPts)
    {
      CoCoA_ASSERT(len(OrigPts) > 0);
      const long NumPPP = len(PPPts);
      bool MadeAnImprovement = true;
      while (MadeAnImprovement)
      {
        MadeAnImprovement = false;
        double MaxGain = 0.000001; // morally zero but must allow for rounding error
        long index_from = 0;     // index of PPPt we move from (initialize to keep compiler quiet)
        long index_to = 0;       // index of PPPt we move to (initialize to keep compiler quiet)
        long IndexOrigPt = 0;    // index of orig point we want to move (initialize to keep compiler quiet)

        // For each PPPt i do ...
        for (long i=0; i < NumPPP; ++i)
        {
          const long Ni = PPPts[i].myNumOrigPts;
          if (Ni == 1) continue;
          // For each index (of OrigPt associated to PPPt i) it do...
          for (long OrigPtIndex: PPPts[i].myOrigPts)
//          for (list<long>::iterator it = PPPts[i].myOrigPts.begin(); it != PPPts[i].myOrigPts.end(); ++it)
          {
            const double diffi = (Ni * PPPts[i].myDistanceSq[OrigPtIndex]) / (Ni-1);
            // For each PPPt j (difft from PPPt i) do...
            for (long j=0; j < NumPPP; ++j)
            {
              if (j == i) continue;
              const long Nj = PPPts[j].myNumOrigPts;
              const double diffj = -(Nj * PPPts[j].myDistanceSq[OrigPtIndex]) / (Nj+1); // careful here: Nj is unsigned!
              const double gain = diffi + diffj;
              if (gain > MaxGain)
              {
                MadeAnImprovement = true;
                MaxGain = gain;
                index_from = i;
                index_to = j;
                IndexOrigPt = OrigPtIndex;
              }
            }
          }
        }

        if (MadeAnImprovement)
        {
          // Move orig point to new PPPt...
          AddAndUpdate(PPPts[index_to], PPPts[index_to], IndexOrigPt, OrigPts);
          EraseAndUpdate(PPPts[index_from], PPPts[index_from], IndexOrigPt, OrigPts);
        }
      }
    }


    // Initialize the list of preprocessed points for subdivision algm.
    void PreprocessSubdivInit(vector<PPP>& PPPts, const vector<PointDbl>& OrigPts)
    {
      CoCoA_ASSERT(PPPts.empty());
      PPPts.resize(1);
      const long NumPts = len(OrigPts);
      for (long i=0; i < NumPts; ++i)
        PPPts[0].myOrigPts.push_back(i);
      PPPts[0].myNumOrigPts = NumPts;
      PPPts[0].myCoords = SubsetAve(PPPts[0].myOrigPts, OrigPts);
      PPPts[0].myDistanceSq.resize(NumPts);
      ComputeDistSq(PPPts[0], OrigPts);
    }


    // Determine if some OrigPt is too far from its corresponding PPPt.
    // If result is true, then WorstOrigIdx and WorstPPPIdx contain the indices
    // of the most distant pair of orig and preprocessed points.
    bool ExistsBadPoint(long& WorstPPPIdx, const vector<PPP>& V)
    {
      double WorstDist = 0;
      const long NumPPP = len(V);
      for (long i=0; i < NumPPP; ++i)
      {
        const long BadIndex = V[i].myWorstOrigPt;
        if (V[i].myDistanceSq[BadIndex] > WorstDist)
        {
          WorstDist = V[i].myDistanceSq[BadIndex];
          WorstPPPIdx = i;
        }
      }
      const long dim = len(V[0].myCoords);
      const double CriticalRadiusSq = square(2*UnitVolBallRadius(dim)); // inefficient if dim is large and num of pts is small
      return WorstDist >= CriticalRadiusSq;
    }


    // Adjoin a new preprocessed point onto the end of V.
    // The new point comprises just one OrigPt, the one with index OrigPtIdx.
    void NewPPP(vector<PPP>& PPPts, const PointDbl& P, long OrigPtIdx)
    {
      const long NumPPPts = len(PPPts);
      PPPts.resize(NumPPPts+1);
      PPP& NewPt(PPPts[NumPPPts]);
      NewPt.myCoords = P;
      NewPt.myOrigPts.push_back(OrigPtIdx);
      NewPt.myNumOrigPts = 1;
      NewPt.myDistanceSq.push_back(0.0);
    }

    // Adjoin a new preprocessed point onto the end of V.
    // The new point comprises just one OrigPt, the one with index OrigPtIdx.
    void NewPPPSubdiv(vector<PPP>& PPPts,
                      const vector<PointDbl>& OrigPts,
                      long OrigPtIdx)
    {
      const long NumPPPts = len(PPPts);
      PointDbl P = OrigPts[OrigPtIdx];
      PPPts.resize(NumPPPts+1);
      PPP& NewPt = PPPts[NumPPPts];
      NewPt.myCoords = P;
      NewPt.myOrigPts.push_back(OrigPtIdx);
      NewPt.myNumOrigPts = 1;
      NewPt.myDistanceSq.resize(len(OrigPts));
      ComputeDistSq(NewPt, OrigPts);
      NewPt.myWorstOrigPt = OrigPtIdx;
    }


  } // end of anonymous namespace


  void PreprocessPts(std::vector<ApproxPts::PointR>& NewPts,
                     std::vector<long>& weights,
                     const std::vector<ApproxPts::PointR>& OrigPts,
                     std::vector<RingElem> tolerance)
  {
    // NOT EXCEPTION CLEAN!!!  But does it matter in this case??
    PreprocessPtsGrid(NewPts, weights, OrigPts, tolerance);
    // Crude heuristic for choosing between "aggr" and "subdiv"... needs improving!!!
    if (len(NewPts) > len(OrigPts)/len(NewPts)) // equiv len(NewPts) > sqrt(len(OrigPts))
      PreprocessPtsAggr(NewPts, weights, OrigPts, tolerance);
    else
      PreprocessPtsSubdiv(NewPts, weights, OrigPts, tolerance);
  }


  void PreprocessCheckArgs(const std::vector<ApproxPts::PointR>& OrigPts,
                           const std::vector<RingElem>& tolerance,
                           const char* const FnName)
  {
    // Some sanity checks on the inputs.
    const long NumPts = len(OrigPts);
    if (NumPts == 0) CoCoA_THROW_ERROR(ERR::ReqNonEmpty, FnName);
    const long dim = len(OrigPts[0]);
    if (dim == 0) CoCoA_THROW_ERROR(ERR::ReqNonEmpty, FnName);
    if (len(tolerance) != dim) CoCoA_THROW_ERROR(ERR::IncompatDims, FnName);

    // Check that all points lie in the same space over the same ring.
    ring R = owner(OrigPts[0][0]);
    if (!IsOrderedDomain(R)) CoCoA_THROW_ERROR(ERR::ReqOrdDom, FnName);
    for (long i=0; i < NumPts; ++i)
    {
      if (len(OrigPts[i]) != dim) CoCoA_THROW_ERROR(ERR::IncompatDims, FnName);
      for (long j=0; j < dim; ++j)
        if (owner(OrigPts[i][j]) != R) CoCoA_THROW_ERROR(ERR::MixedRings, FnName);
    }
    for (long j=0; j < dim; ++j)
    {
      if (owner(tolerance[j]) != R) CoCoA_THROW_ERROR(ERR::MixedRings, FnName);
      if (tolerance[j] < 0) CoCoA_THROW_ERROR(ERR::ReqNonNegative, FnName);
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  // Below are the three main algorithms.
  // First the "grid" algorithm.

  void PreprocessPtsGrid(std::vector<ApproxPts::PointR>& NewPts,
                         std::vector<long>& weights,
                         const std::vector<ApproxPts::PointR>& OrigPts,
                         std::vector<RingElem> tolerance)
  {
    PreprocessCheckArgs(OrigPts, tolerance, "PreprocessPtsGrid");
    const long NumPts = len(OrigPts);

    // Normalization effectively makes the tolerance in each direction equal to 1.
    ZeroToleranceHack(OrigPts, tolerance);
    vector<PointDbl> ScaledPts; // value is set by NormalizePts
    NormalizePts(ScaledPts, OrigPts, tolerance);

    vector<PPP> PPPts; PPPts.reserve(NumPts);
    // For each OrigPts[i] do...
    for (long i = 0; i < NumPts; ++i)
    {
      const PointDbl P = GridNearPoint(ScaledPts[i]);
      const vector<PPP>::iterator pos = find_if(PPPts.begin(), PPPts.end(), CoordsEq(P));
      if (pos != PPPts.end())
      {
        pos->myOrigPts.push_back(i);  // adjoin to existing class
        ++(pos->myNumOrigPts);        // increment myNumOrigPts
      }
      else
        NewPPP(PPPts, P, i);          // add a new PPPt
    }

    const long NumPPPts = len(PPPts);
    NewPts.resize(NumPPPts);
    weights.resize(NumPPPts);
    for (long k = 0; k < NumPPPts; ++k)
    {
      NewPts[k] = SubsetAve(PPPts[k].myOrigPts, OrigPts);
      weights[k] = PPPts[k].myNumOrigPts;
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  // The "aggregative" algorithm.

  void PreprocessPtsAggr(std::vector<ApproxPts::PointR>& NewPts,
                         std::vector<long>& weights,
                         const std::vector<ApproxPts::PointR>& OrigPts,
                         std::vector<RingElem> tolerance)
  {
    PreprocessCheckArgs(OrigPts, tolerance, "PreprocessPtsAggr");
    const long NumPts = len(OrigPts);
    const long dim = len(tolerance);

    // Normalization effectively makes the tolerance in each direction equal to 1.
    ZeroToleranceHack(OrigPts, tolerance);
    vector<PointDbl> ScaledPts; // value set by NormalizePts
    NormalizePts(ScaledPts, OrigPts, tolerance);

    vector<PPP> PPPts(NumPts);
    for (long i=0; i < NumPts; ++i)
    {
      PPPts[i].myCoords = ScaledPts[i];
      PPPts[i].myNumOrigPts = 1;
      PPPts[i].myOrigPts.push_back(i);
    }

    // nbrs contains all pairs which could conceivably collapse together
    // (i.e. only pairs separated by at most the CriticalRadius)
    const double CriticalRadiusSq = square(2*UnitVolBallRadius(dim));
    list<NearPair> nbrs;
    for (long i=0; i < NumPts; ++i)
    {
      for (long j=i+1; j < NumPts; ++j)
      {
        const double d2 = L2NormSq(PPPts[i].myCoords, PPPts[j].myCoords);
        if (d2 < 4*CriticalRadiusSq)
          nbrs.push_back(NearPair(i,j,d2));
      }
    }

    nbrs.sort(LessNearPair);

    PointDbl CofG(dim);
    while (!nbrs.empty() && nbrs.begin()->mySepSq < 4*CriticalRadiusSq)
    {
      const long i = nbrs.front().myLowerIndex;
      const long j = nbrs.front().myHigherIndex;
      nbrs.erase(nbrs.begin()); // delete first elem of nbrs

      const PPP& Pi = PPPts[i];
      const PPP& Pj = PPPts[j];
      // Compute CofG of the union of Pi an Pj...
      WeightedAve(CofG, len(Pi.myOrigPts), Pi.myCoords, len(Pj.myOrigPts), Pj.myCoords);

      // Check that all original points in the union are close to CofG.
      bool OK = true;
      for (long k: Pi.myOrigPts)
      {
        OK &= (L2NormSq(CofG, ScaledPts[k]) < CriticalRadiusSq);
      }
      for (long k: Pj.myOrigPts)
      {
        OK &= (L2NormSq(CofG, ScaledPts[k]) < CriticalRadiusSq);
      }
      if (!OK) continue; // We cannot unite these two points, so go on to next pair of points.

      // The union is fine, so accept it.  We put all data into i-th PPP
      // and empty the data of the j-th PPP.
      {
        list<long>& PPPOrigPts = PPPts[i].myOrigPts;
        PPPOrigPts.splice(PPPOrigPts.begin(), PPPts[j].myOrigPts);
        PPPts[i].myNumOrigPts += PPPts[j].myNumOrigPts;
        PPPts[j].myNumOrigPts = 0;
      }
      nbrs.remove_if(NearPairContains(j)); // since Pj has been eliminated.
      // Update all entries referring to i-th point.
//      for (list<NearPair>::iterator nbr = nbrs.begin(); nbr != nbrs.end(); ++nbr)
      for (NearPair& nbr: nbrs)
      {
        if (!NearPairContains(i)(nbr)) continue;
        long other = nbr.myLowerIndex;
        if (other == i) other = nbr.myHigherIndex;
        nbr.mySepSq = L2NormSq(CofG, PPPts[other].myCoords);
      }
      nbrs.sort(LessNearPair); // put the list back into order -- SLUG SLUG SLUG!!!
      swap(PPPts[i].myCoords, CofG); // really an assignment to PPPts[i].myCoords
    }

    // There are no further possible mergings, so copy answer into NewPts & weights.
    NewPts.clear();
    weights.clear();
    for (long i=0; i < NumPts; ++i)
    {
      const long N = PPPts[i].myNumOrigPts;
      if (N == 0) continue;
      NewPts.push_back(SubsetAve(PPPts[i].myOrigPts, OrigPts));
      weights.push_back(N);
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  // The "subdivision" algorithm.

  void PreprocessPtsSubdiv(std::vector<ApproxPts::PointR>& NewPts,
                           std::vector<long>& weights,
                           const std::vector<ApproxPts::PointR>& OrigPts,
                           std::vector<RingElem> tolerance)
  {
    PreprocessCheckArgs(OrigPts, tolerance, "PreprocessPtsSubdiv");
    const long NumPts = len(OrigPts);

    // Normalization effectively makes the tolerance in each direction equal to 1.
    ZeroToleranceHack(OrigPts, tolerance);
    vector<PointDbl> ScaledPts; // value is set by NormalizePts
    NormalizePts(ScaledPts, OrigPts, tolerance);

    // Initially all orig points are associated to a single preprocessed point.
    vector<PPP> PPPts; PPPts.reserve(NumPts);
    PreprocessSubdivInit(PPPts, ScaledPts);

    // Main loop.
    long WorstPPPIdx = 0; // Init value just to keep compiler quiet.
    // True value is set in call to ExistsBadPoint (see next line)
    while (ExistsBadPoint(WorstPPPIdx, PPPts))
    {
      const long WorstOrigIdx = PPPts[WorstPPPIdx].myWorstOrigPt;
      NewPPPSubdiv(PPPts, ScaledPts, WorstOrigIdx);
      EraseAndUpdate(PPPts[WorstPPPIdx], PPPts[WorstPPPIdx], WorstOrigIdx, ScaledPts);
      redistribute(PPPts, ScaledPts);
    }

    // Copy and denormalize result; add the weight of each PPP.
    NewPts.clear(); NewPts.reserve(len(PPPts));
    weights.clear(); weights.reserve(len(PPPts));
    for (long i=0; i < len(PPPts); ++i)
    {
      const long N = PPPts[i].myNumOrigPts;
      NewPts.push_back(SubsetAve(PPPts[i].myOrigPts, OrigPts));
      weights.push_back(N);
    }
  }


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/ApproxPts.C,v 1.36 2024/07/15 16:57:51 bigatti Exp $
// $Log: ApproxPts.C,v $
// Revision 1.36  2024/07/15 16:57:51  bigatti
// Summary: Error codes:  Renamed  ERR::Empty --> ERR::ReqNonEmpty
//
// Revision 1.35  2024/07/02 20:42:06  abbott
// Summary: Renaming of errors (redmine 308)
//
// Revision 1.34  2024/07/02 15:29:22  bigatti
// Summary: Changed error codes: LogZero into ReqNonZero
//
// and Not... into ReqNonNegative, ReqNonNegativeGrading, ReqPositive, ReqPositiveGrading
// Specialized some BadArg with explicit codes
// Tidied up some CoCoA_THROW_ERROR second argument (string)
//
// Revision 1.33  2023/03/13 19:55:54  abbott
// Summary: Changed to iterator for loop (redmine 1346)
//
// Revision 1.32  2022/02/18 14:11:52  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.31  2021/07/19 13:10:28  abbott
// Summary: Removed commented out code; changed one more for loop into new syntax.
//
// Revision 1.30  2021/07/19 11:06:16  abbott
// Summary: Indentation
//
// Revision 1.29  2021/01/07 15:07:02  abbott
// Summary: Corrected copyright
//
// Revision 1.28  2020/06/17 15:49:21  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.27  2020/02/18 11:25:12  abbott
// Summary: Redmine 1346: new for loop syntax
//
// Revision 1.26  2017/02/15 12:23:33  abbott
// Summary: Added square fn (for double)
//
// Revision 1.25  2014/05/10 20:46:51  abbott
// Summary: Corrected arg checking; also now a separate fn
// Author: JAA
//
// Revision 1.24  2014/05/06 13:13:10  abbott
// Summary: Replaced assertions by "always" arg checks
// Author: JAA
//
// Revision 1.23  2014/04/30 16:03:53  abbott
// Summary: Eliminated use of sqrt (and corrected the program logic!)
// Author: JAA
//
// Revision 1.22  2013/03/27 18:24:33  abbott
// Added approx point preprocessing to C5; also changed names of the fns, and updated doc.
//
// Revision 1.21  2011/08/24 10:24:17  bigatti
// -- renamed QQ --> BigRat
// -- sorted #include
//
// Revision 1.20  2011/08/14 15:52:17  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.19  2011/03/11 12:01:52  bigatti
// -- changed size --> len
//
// Revision 1.18  2011/03/11 11:08:27  bigatti
// -- changed size_t --> long
// -- changed size --> len
//
// Revision 1.17  2009/12/11 11:46:32  abbott
// Changed fn  convert  into  IsConvertible.
// Added template procedure  convert.
// New version because change is not backward compatible.
//
// Revision 1.16  2009/10/29 19:26:20  abbott
// Cleaned up include directives: now in alphabetical order, removed unneeded ones.
//
// Revision 1.15  2009/09/25 13:12:38  abbott
// Added initial value for WorstPPPIdx just to keep compiler quiet.
//
// Revision 1.14  2009/07/02 16:33:59  abbott
// Switched to using new IsRational interface (using a QQ value).
//
// Revision 1.13  2008/11/24 17:12:14  abbott
// Final tidying: renamed an internal fn, removed useless includes.
//
// Revision 1.12  2008/11/23 18:58:32  abbott
// Major overhaul to preprocessing and SOI/NBM code.
// Split SOI/NBM off into a separate file.
// Preprocessing is now "rational" (but internally guided by C++ doubles).
// SOI/NBM now each have 3 similar interfaces: one purely rational, one for
// input which is represented as doubles, and one which converts the input
// to RingTwinFloat values and produces a result which is over some RingTwinFloat
// (the precision is increased automatically until an answer is obtained).
//
// Revision 1.11  2008/11/20 10:47:30  abbott
// Now the code actually computes the list of weights.
//
// Revision 1.10  2008/11/20 10:08:15  abbott
// The preprocessing fns now return also the weights associated with the representatives.
//
// Revision 1.9  2008/11/06 12:50:44  abbott
// Moved definitions of square and round to utils.H from ApproxPts.H
//
// Revision 1.8  2008/10/07 15:45:22  abbott
// Changed ErrorInfo objects so they include the name of their own error ID.
// Changed catch statements to catch const objects.
// Removed calls to the member fn which accessed the error ID member of an
// ErrorInfo; now you simply compare directly with the error ID (makes the
// code easier to read).
//
// Revision 1.7  2008/09/12 13:28:43  bigatti
// -- new: NBM implementation
//
// Revision 1.6  2008/06/04 18:27:37  abbott
// Modified the server interface for "SOI": it now accepts a 3rd arg (gamma).
//
// Revision 1.5  2008/05/30 14:20:43  abbott
// SOI now returns also the "almost vanishing" polynomials.
//
// Revision 1.4  2008/05/30 12:51:08  abbott
// Added a comment.
//
// Revision 1.3  2008/05/29 15:46:29  bigatti
// -- added Approximate Border Basis (by Abbott,Torrente)
//
// Revision 1.2  2007/10/30 17:14:08  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.8  2007/03/08 18:22:30  cocoa
// Just whitespace cleaning.
//
// Revision 1.7  2006/12/21 13:48:33  cocoa
// Made all increment/decrement calls prefix (except where the must be postfix).
//
// Revision 1.6  2006/11/22 14:43:32  cocoa
// -- minor cleaning (indicated by Intel compiler)
//
// Revision 1.5  2006/11/20 15:53:38  cocoa
// Minor cleaning.  Improved comments.
//
// Revision 1.4  2006/10/31 14:26:03  cocoa
// Added some missing consts; should make the code a little easier to comprehend.
//
// Revision 1.3  2006/10/06 14:04:15  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.2  2006/06/21 17:05:47  cocoa
// Major overhaul of approx point preprocessing algms.
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.3  2006/05/29 16:16:47  cocoa
// Added "disj" preprocessing algorithm.
//
// Revision 1.2  2006/05/22 15:52:16  cocoa
// Added preprocess-disg algorithm to ApproxPts.
// Sundry minor improvements.
//
// Revision 1.1  2006/05/12 13:16:30  cocoa
// Added functions for preprocessing approximate points.
//
//
