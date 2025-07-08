//   Copyright (c)  2017  John Abbott,  Anna M. Bigatti

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

#include "CoCoA/RootBound.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/BigRatOps.H"
#include "CoCoA/NumTheory-gcd.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/SparsePolyIter.H"
#include "CoCoA/SparsePolyOps-RingElem.H"
//#include "CoCoA/SparsePolyRing.H" // from SparsePolyOps-RingElem.H
#include "CoCoA/VectorOps.H"
#include "CoCoA/error.H"
#include "CoCoA/verbose.H"

#include<algorithm>
using std::min;
#include <cmath>
using std::abs;
using std::log;
#include <iostream>
using std::endl;
#include <vector>
using std::vector;

namespace CoCoA
{

  namespace // anonymous
  {

    // Throws if f is not univariate in char 0 with deg >= 1.
    // WARNING: does not check that coeffs are rationals!
    void RootBound_CheckArg(ConstRefRingElem f, const ErrorContext& FnName)
    {
      const ring& P = owner(f);
      if (!IsPolyRing(P) || !IsZero(characteristic(P)))
        CoCoA_THROW_ERROR_WITH_CONTEXT3(ERR::BadArg, "expected univariate poly in characteristic 0", FnName);
      if (IsZero(f) || deg(f) == 0)
        CoCoA_THROW_ERROR_WITH_CONTEXT3(ERR::BadArg, "must have deg >= 1", FnName);
      if (UnivariateIndetIndex(f) < 0)
        CoCoA_THROW_ERROR_WITH_CONTEXT2(ERR::ReqUnivariate, FnName); // not univariate
    }


    // This fn is now obsolete; fn below is better (but does use exp(...))
    //
    // // upper bound for 1/(2^(1/d)-1) for d integer >= 2
    // // For d > 4 the formula (185*d-64)/128 requires proof.
    // BigRat BirkhoffScaleFactor(long d)
    // {
    //   if (d < 2)  CoCoA_THROW_ERROR2(ERR::BadArg,"degree must be at least 2");
    //   if (d > 4)  return BigRat(185*d-64, 128); // error less than 0.2%
    //   switch (d)
    //   {
    //   case 2: return BigRat(619,256);
    //   case 3: return BigRat(493,128);
    //   case 4: return BigRat(677,128);
    //   }
    //   CoCoA_THROW_ERROR1(ERR::ShouldNeverGetHere);
    //   return BigRat(-1,1); // just to keep compiler quiet
    // }


    // Upper bound for Birkhoff's scale factor (excess is about 0.01%-0.1%)
    double BirkhoffScaleFactor2(long d)
    {
      if (d < 2)  CoCoA_THROW_ERROR2(ERR::BadArg,"degre must be at least 2");
      const double log2 = std::log(2.0);
      if (d >= 7)  return (1+1.0/1024)*(d/log2-0.5);
      return (1+1.0/1024)/(exp(log2/d)-1);
    }


    long CoeffSizeEstimate(ConstRefRingElem f)
    {
      CoCoA_ASSERT(IsZZ(CoeffRing(owner(f))));
      long ans = 0;
      for (SparsePolyIter it=BeginIter(f); !IsEnded(it); ++it)
      {
        ans += 1+FloorLog2(ConvertTo<BigInt>(coeff(it)));
      }
      return ans;
    }


    // Make f "primitive": all coeffs are integer, content is 1, LC is positive.
    // Divide f by highest power x^k so that ConstantCoeff(f) is non-zero.
    // Then if f(x) = g(x^k) some k then it returns g(x)
    RingElem RootBound_preprocess(ConstRefRingElem f)
    {
      RootBound_CheckArg(f, CoCoA_ERROR_CONTEXT);
      
      VerboseLog VERBOSE("RootBound_preprocess");
      const PolyRing ZZx = NewPolyRing(RingZZ(), symbols("x"));
      const PPMonoidElem x = indet(PPM(ZZx),0);
      long GcdExp = 0;
      long LastExp = -1;
      BigInt CommonDenom(1);
      BigInt CommonNumer;
      for (SparsePolyIter it=BeginIter(f); !IsEnded(it); ++it)
      {
        const long d = deg(PP(it));
        if (LastExp != -1 && GcdExp != 1)
          GcdExp = gcd(GcdExp, LastExp-d);
        LastExp = d;
        const BigRat c = ConvertTo<BigRat>(coeff(it));
        CommonNumer = gcd(CommonNumer, num(c));
        CommonDenom = lcm(CommonDenom, den(c));
      }
      // Trick to make sure result has positive leading coeff.
      if (sign(LC(f)) < 0)
        CommonNumer = -CommonNumer;
    
      RingElem ans(ZZx);
      for (SparsePolyIter it=BeginIter(f); !IsEnded(it); ++it)
      {
        const long d = deg(PP(it));
        const BigRat c = ConvertTo<BigRat>(coeff(it));
        const RingElem IntCoeff(RingZZ(), (num(c)/CommonNumer)*(CommonDenom/den(c)));
        PushBack(ans, IntCoeff, power(x, (d-LastExp)/GcdExp));
      }
      if (VerbosityLevel() >= 50)
      {
        VERBOSE(50) << "Removed power of x: " << LastExp << endl;
        VERBOSE(50) << "Poly was in x^" << GcdExp << endl;
      }
      return ans;
    }

  
    // Result is approximately exp(d) expressed as a rational of the form
    // n*2^k where n is an integer less than 256 (and k is integer, of course!)
    BigRat ApproxExp(double d)
    {
      const double ln2 = std::log(2.0);
      const long pwr2 = (-8) + static_cast<long>(std::floor(d/ln2));

      const double excess = d - pwr2*ln2;
      const int numer = static_cast<int>(std::ceil(exp(excess)));
      if (pwr2 < 0)  return BigRat(numer,power(2,-pwr2));
      return BigRat(numer*power(2,pwr2), 1);
    }
  

    BigRat ApproxRoot(const BigRat& q, int n)
    {
      if (n == 1)  return q;
      return ApproxExp(LogAbs(q)/n);
    }


    constexpr double LogZero = -1.0e15; // "minus infinity"


    // Returns vector<double> such that entry k contains
    // log(abs(a_k/a_d))  where  a_d  is the leading coeff.
    // if a_k is 0 then entry is LogZero (see defn above)
    vector<double> LogCoeffVec(ConstRefRingElem f)
    {
      CoCoA_ASSERT(!IsZero(f) && deg(f) > 0);
      const int d = deg(f);
      vector<double> LogCoeff(d, LogZero);
      const double LogLCF = LogAbs(ConvertTo<BigRat>(LC(f)));
      for (SparsePolyIter it=++BeginIter(f); !IsEnded(it); ++it)
      {
        const int i = deg(PP(it));
        LogCoeff[i] = LogAbs(ConvertTo<BigRat>(coeff(it))) - LogLCF;
      }
      return LogCoeff;
    }


    // (exact) Root bound for polys of degree 1
    BigRat RootBound_deg1(ConstRefRingElem f)
    {
      CoCoA_ASSERT(!IsZero(f) && deg(f) == 1);
      if (NumTerms(f) == 1)  return BigRat(0);
      const BigRat lcf = ConvertTo<BigRat>(LC(f));
      return abs(ConvertTo<BigRat>(coeff(++BeginIter(f)))/lcf);
    }

  } // end of namespace anonymous
  

  //------------------------------------------------------------------
  // Classical Cauchy bound.
  // Ref: Wikipedia, "properties of polynomial roots"
  
  double LogRootBound_Cauchy(const vector<double>& LogCoeff)
  {
    CoCoA_ASSERT(!LogCoeff.empty());
//??    VerboseLog VERBOSE("LogRootBound_Cauchy");
    const int d = len(LogCoeff);  // must have d > 0
    if (d == 1)  return LogCoeff[0];

    double MaxVal = LogZero;
    for (int i=0; i < d; ++i)
    {
      if (LogCoeff[i] == LogZero)  continue;
      if (LogCoeff[i] > MaxVal)  MaxVal = LogCoeff[i];
//      MaxVal = std::max(LogCoeff[i], MaxVal);
    }

    const double ans = log(1 + ApproxExp(MaxVal));  // any cleverer way?
//??    VERBOSE(51) << "RETURN " << ans << endl;
    return ans;
  }


  BigRat RootBound_Cauchy(ConstRefRingElem f)
  {
    RootBound_CheckArg(f, CoCoA_ERROR_CONTEXT);

    const int d = deg(f);
    if (d == 1)  return RootBound_deg1(f);
    if (IsMonomial(f))  return BigRat(0);
    const double LogBound = LogRootBound_Cauchy(LogCoeffVec(f));
    return ApproxExp(LogBound);
  }

  //------------------------------------------------------------------
  // Classical Lagrange bound -- see also the LMS bound below!
  // Ref: Wikipedia, "properties of polynomial roots"

  double LogRootBound_Lagrange(const vector<double>& LogCoeff)
  {
    CoCoA_ASSERT(!LogCoeff.empty());
    VerboseLog VERBOSE("LogRootBound_Lagrange");
    const int d = len(LogCoeff);  // must have d > 0
    if (d == 1)  return LogCoeff[0];

    double MaxVal = LogZero;
    for (int i=0; i < d; ++i)
    {
      if (LogCoeff[i] == LogZero)  continue;
      if (LogCoeff[i] > MaxVal)  MaxVal = LogCoeff[i];
    }

//    if (MaxVal == LogZero) return LogZero;
    if (std::log(d) + MaxVal < 0)  return 0; // 0 = log(1)
    BigRat SumOfBigCoeffs;
    const double LWB = MaxVal - std::log(d) - 7;  // we shall ignore "negligible" values
    VERBOSE(55) << "ignoring coeffs with logs below " << LWB << endl;
    for (int i=0; i < d; ++i)
    {
      if (LogCoeff[i] < LWB)  continue;
      SumOfBigCoeffs += ApproxExp(LogCoeff[i]);
    }

    const double ans = std::max(0.0, log(SumOfBigCoeffs));
    VERBOSE(52) << "RETURN " << ans << endl;
    return ans;
  }


  BigRat RootBound_Lagrange(ConstRefRingElem f)
  {
    RootBound_CheckArg(f, CoCoA_ERROR_CONTEXT);

    const int d = deg(f);
    if (d == 1)  return RootBound_deg1(f);
    if (IsMonomial(f))  return BigRat(0);
    const double LogBound = LogRootBound_Lagrange(LogCoeffVec(f));
    return ApproxExp(LogBound);
  }

  //------------------------------------------------------------------
  // This was the root bound used by Zassenhaus (1969)

  double LogRootBound_Birkhoff(const vector<double>& LogCoeff)
  {
    CoCoA_ASSERT(!LogCoeff.empty());
    VerboseLog VERBOSE("LogRootBound_Birkhoff");
    const int d = len(LogCoeff);  // must have d > 0
    if (d == 1)  return LogCoeff[0];

    double LogBinomial = 0.0;
    double MaxVal = LogZero;
    VERBOSE(55) << "MaxVal = " << MaxVal << endl;
    for (int i=0; i < d; ++i)
    {
      if (LogCoeff[i] != LogZero)
      {
        const double val = (LogCoeff[i]-LogBinomial)/(d-i);
        if (val > MaxVal)  MaxVal = val;
//      MaxVal = std::max(val, MaxVal);
      }
      LogBinomial += std::log(d-i) - std::log(i+1);
    }

    const double ScaleFactor = BirkhoffScaleFactor2(d);
    VERBOSE(55) << "ScaleFactor = " << ScaleFactor << endl;
    VERBOSE(55) << "Max log val = " << MaxVal << endl;
    VERBOSE(52) << "RETURN " << std::log(ScaleFactor) + MaxVal << endl;
    return std::log(ScaleFactor) + MaxVal;
  }


  BigRat RootBound_Birkhoff(ConstRefRingElem f)
  {
    RootBound_CheckArg(f, CoCoA_ERROR_CONTEXT);

//    VerboseLog VERBOSE("RootBound_Birkhoff");
    const int d = deg(f);
    if (d == 1)  return RootBound_deg1(f);
    if (IsMonomial(f))  return BigRat(0);
    const double LogBound = LogRootBound_Birkhoff(LogCoeffVec(f));
    return ApproxExp(LogBound);
    
//     vector<double> LogCoeff(d, LogZero);
//     const double LogLCF = LogAbs(ConvertTo<BigRat>(LC(f)));
//     for (SparsePolyIter it=++BeginIter(f); !IsEnded(it); ++it)
//     {
//       const int i = deg(PP(it));
//       LogCoeff[i] = LogAbs(ConvertTo<BigRat>(coeff(it)))-LogLCF;
//     }
//     VERBOSE(51) << "LogCoeff = " << LogCoeff << endl;
//     double LogBinomial = 0.0;
//     double MaxVal = LogZero;
//     for (int i=0; i < d; ++i)
//     {
//       if (LogCoeff[i] == LogZero)  continue;
//       const double val = (LogCoeff[i]-LogBinomial)/(d-i);
//       MaxVal = std::max(val, MaxVal);
//       LogBinomial += std::log(d-i) - std::log(i+1);
//     }

//     const BigRat ScaleFactor = BirkhoffScaleFactor(d);
//     VERBOSE(50) << "ScaleFactor = " << ScaleFactor << endl;
//     VERBOSE(50) << "Max log val = " << MaxVal << endl;
//     return ScaleFactor *  ApproxExp(MaxVal);
  }


//   //------------------------------------------------------------------
//   // 2019-04-26: found this fn is an old (2011) test program...
//   // Originally from a 1969 version of a Knuth book (Seminumerical Algorithms)
//   // LMS bound (below) is always <= Knuth bound.  So Knuth bound is superseded by LMS.

//   double LogRootBound_Knuth(const vector<double>& LogCoeff)
//   {
//     CoCoA_ASSERT(!LogCoeff.empty());
//     const int d = len(LogCoeff);  // must have d > 0
//     if (d == 1)  return LogCoeff[0];

//     double MaxVal = LogZero;
//     for (int i=0; i < d; ++i)
//     {
//       if (LogCoeff[i] == LogZero)  continue;
//       const double val = LogCoeff[i]/(d-i);
//       if (val > MaxVal)  MaxVal = val;
// //      MaxVal = std::max(val, MaxVal);
//     }

//     return std::log(2.0) + MaxVal;
//   }

//   BigRat RootBound_Knuth(ConstRefRingElem f)
//   {
//     RootBound_CheckArg(f, CoCoA_ERROR_CONTEXT);

//     const int d = deg(f);
//     if (d == 1)  return RootBound_deg1(f);
//     if (IsMonomial(f))  return BigRat(0);
//     const double LogBound = LogRootBound_Knuth(LogCoeffVec(f));
//     return ApproxExp(LogBound);
//   }
  

  //------------------------------------------------------------------
  // LMS = Lagrange-Mignotte-Stefanescu, and slightly better than Fujiwara (1916)
  // Supersedes bounds from Fujiwara (1916) and Knuth (Seminumerical Algms, 1969)

  double LogRootBound_LMS(const vector<double>& LogCoeff)
  {
    CoCoA_ASSERT(!LogCoeff.empty());
    VerboseLog VERBOSE("RootBound_LMS");
    const int d = len(LogCoeff);  // must have d > 0
    if (d == 1)  return LogCoeff[0];

    double max1 = LogZero;
    double max2 = LogZero;
    for (int i=0; i < d; ++i)
    {
      if (LogCoeff[i] == LogZero)  continue;
      const double val = LogCoeff[i]/(d-i);
      if (val < max2)  continue;
      if (val < max1)  { max2 = val; continue; }
      max2 = max1; max1 = val;
    }
    VERBOSE(55) << "max1 = " << max1 << "   max2 = " << max2 << endl;
    if (max2 == LogZero)  return max1;
    if (max2 < max1 - 36)  return max1; // NB 36 = 51*ln(2)
    return max1 + log(1+ApproxExp(max2-max1));
  }


  // LMS = Lagrange-Mignotte-Stefanescu
  BigRat RootBound_LMS(ConstRefRingElem f)
  {
    RootBound_CheckArg(f, CoCoA_ERROR_CONTEXT);

    const int d = deg(f);
    if (d == 1)  return RootBound_deg1(f);
    if (IsMonomial(f))  return BigRat(0);
    const double LogBound = LogRootBound_LMS(LogCoeffVec(f));
    return ApproxExp(LogBound);

//     VerboseLog VERBOSE("RootBound_LMS");
//     const int d = deg(f);
//     if (d == 1)  return RootBound_deg1(f);
//     if (IsMonomial(f))  return BigRat(0);
//     vector<double> LogCoeff(d, LogZero);
//     const double LogLCF = LogAbs(ConvertTo<BigRat>(LC(f)));
//     for (SparsePolyIter it=++BeginIter(f); !IsEnded(it); ++it)
//     {
//       const int i = deg(PP(it));
//       LogCoeff[i] = LogAbs(ConvertTo<BigRat>(coeff(it))) - LogLCF;
//     }
//     VERBOSE(51) << "LogCoeff = " << LogCoeff << endl;

//     double max1 = LogZero;
//     double max2 = LogZero;
//     for (int i=0; i < d; ++i)
//     {
//       if (LogCoeff[i] == LogZero)  continue;
//       double val = LogCoeff[i]/(d-i);
//       if (val < max2)  continue;
//       if (val < max1)  { max2 = val; continue; }
//       max2 = max1; max1 = val;
//     }
//     VERBOSE(50) << "max1 = " << max1 << "   max2 = " << max2 << endl;
//     if (max2 == LogZero)  return ApproxExp(max1);
//     return ApproxExp(max1) + ApproxExp(max2);
  }


  //------------------------------------------------------------------

  // LogRootBound_simple is the best of the "simple" bounds
  double LogRootBound_simple(const vector<double>& LogCoeff)
  {
    const double LogBound_Cauchy = LogRootBound_Cauchy(LogCoeff);
    const double LogBound_Lagrange = LogRootBound_Lagrange(LogCoeff);
    const double LogBound_Birkhoff = LogRootBound_Birkhoff(LogCoeff);
    const double LogBound_LMS = LogRootBound_LMS(LogCoeff);
    VerboseLog VERBOSE("LogRootBound_simple");
    VERBOSE(51) << "LogBound_Cauchy = "  << LogBound_Cauchy << endl;
    VERBOSE(51) << "LogBound_Lagrange = "  << LogBound_Lagrange << endl;
    VERBOSE(51) << "LogBound_Birkhoff = "  << LogBound_Birkhoff << endl;
    VERBOSE(51) << "LogBound_LMS = "  << LogBound_LMS << endl;
    const double LogBound_min =  min(min(LogBound_Cauchy, LogBound_Lagrange),
                                     min(LogBound_Birkhoff, LogBound_LMS));
    VERBOSE(51) << "LogBound_min = "  << LogBound_min << endl;
    return LogBound_min;
  }

  BigRat RootBound_simple(ConstRefRingElem f)
  {
    RootBound_CheckArg(f, CoCoA_ERROR_CONTEXT);

    const int d = deg(f);  // must have d > 0
    if (d == 1)  return RootBound_deg1(f);
    if (IsMonomial(f))  return BigRat(0);
    const double LogBound = LogRootBound_simple(LogCoeffVec(f));
    return ApproxExp(LogBound);
  }

  //------------------------------------------------------------------

  BigRat RootBound(ConstRefRingElem f, long NumGraeffeIters)
  {
    RootBound_CheckArg(f, CoCoA_ERROR_CONTEXT);
    if (NumGraeffeIters < -1 || NumGraeffeIters > 25)
      CoCoA_THROW_ERROR2(ERR::OutOfRange, "number of iters must be in 0..25");
    
    VerboseLog VERBOSE("RootBound");
    if (IsMonomial(f))
    {
      VERBOSE(50) << "Input is monomial --> ans is 0" << endl;
      return BigRat(0);
    }
    RingElem g = RootBound_preprocess(f);
    const long f_step = deg(f) - deg(PP(++BeginIter(f)));
    const long d = deg(g);
    VERBOSE(40) << "Preprocessed poly has deg: " << d << endl;
    const long g_step = d - deg(PP(++BeginIter(g)));
    VERBOSE(40) << "Input was poly in x^" << (f_step/g_step) << endl;
    if (d == 1)
    {
      VERBOSE(40) << "Preprocessed poly is linear" << endl;
      const BigRat lcg = ConvertTo<BigRat>(LC(g));
      const BigRat ConstCoeff = ConvertTo<BigRat>(coeff(++BeginIter(g)));
      return ApproxRoot(abs(ConstCoeff/lcg), f_step/g_step);
    }
    BigRat B = RootBound_simple(g);
//    BigRat RB_LMS = RootBound_LMS(g);
//    BigRat RB_B = RootBound_Birkhoff(g);
    // VERBOSE(40) << "LMS bound: " << RB_LMS << endl;
    // VERBOSE(40) << "Birkhoff bound: " << RB_B << endl;
    // BigRat B = min(RB_LMS, RB_B);
    VERBOSE(40) << "Overall bound: " << B << endl;

    int pwr = 1;
    const long CoeffSize = CoeffSizeEstimate(g);
    VERBOSE(40) << "CoeffSizeEstimate = " << CoeffSize << endl;
    long MaxIters = NumGraeffeIters;
    if (NumGraeffeIters < 0)
      MaxIters = min(5L, 22-FloorLog2(CoeffSize));

    for (long NumIters = 1; NumIters <= MaxIters; ++NumIters)
    {
      VERBOSE(40) << "Graeffe loop [" << NumIters << "]  start" << endl;
      pwr *= 2;
      g = graeffe(g);
      const BigRat NewB = ApproxRoot(RootBound_simple(g), pwr);
      VERBOSE(40) << "Graeffe loop [" << NumIters << "]  Bound for transformed poly: " << NewB << endl << endl;
      B = min(B, NewB);
      VERBOSE(40) << "Graeffe loop [" << NumIters << "]  Overall bound: " << B << endl << endl;
    }
    VERBOSE(40) << "After graeffe loop B = " << B << endl;
    if (deg(f) != deg(g))
    {
      B = ApproxRoot(B, f_step/g_step);
      VERBOSE(40) << "Take " << f_step/g_step << " root..  B= " << B << endl;
    }
    return B;
  }


  double LogRootBound(ConstRefRingElem f, long NumGraeffeIters)
  {
    CoCoA_ASSERT(NumGraeffeIters >= -1 && NumGraeffeIters <= 25);
    VerboseLog VERBOSE("LogRootBound");
    if (IsMonomial(f))
    {
      VERBOSE(50) << "Input is monomial --> ans is 0" << endl;
      return LogZero;
    }
    RingElem g = RootBound_preprocess(f);
    const long f_step = deg(f) - deg(PP(++BeginIter(f)));
    const long d = deg(g);
    VERBOSE(40) << "Preprocessed poly has deg: " << d << endl;
    const long g_step = d - deg(PP(++BeginIter(g)));
    VERBOSE(40) << "Input was poly in x^" << (f_step/g_step) << endl;
     
    vector<double> LogCoeffs = LogCoeffVec(g); //not const; updated in graeffe loop
    if (d == 1)
    {
      VERBOSE(40) << "Preprocessed poly was linear" << endl;
      return LogCoeffs[0]/(f_step/g_step);
    }
    double LogBound = LogRootBound_simple(LogCoeffs);
    VERBOSE(40) << "Initial log bound: " << LogBound << endl;

    VERBOSE(40) << "CoeffSizeEstimate = " << CoeffSizeEstimate(g) << endl;
    long MaxIters = NumGraeffeIters;
    if (NumGraeffeIters < 0)
      MaxIters = min(5L, 22-FloorLog2(CoeffSizeEstimate(g)));

    int pwr = 1;
    for (long NumIters = 1; NumIters <= MaxIters; ++NumIters)
    {
      VERBOSE(40) << "Graeffe loop [" << NumIters << "]  start" << endl;
      pwr *= 2;
      g = graeffe(g);
      LogCoeffs = LogCoeffVec(g);
      const double NewLogBound = LogRootBound_simple(LogCoeffs)/pwr;
      VERBOSE(45) << "Graeffe loop [" << NumIters << "]  new log bound: " << NewLogBound << endl;
      LogBound = std::min(LogBound, NewLogBound);
      VERBOSE(45) << "Graeffe loop [" << NumIters << "]  Overall log bound: " << LogBound << endl << endl;
    }
    VERBOSE(40) << "After graeffe loop Log bound = " << LogBound << endl;
    if (deg(f) != deg(g))
    {
      LogBound = LogBound/(f_step/g_step);
    VERBOSE(40) << "Take " << f_step/g_step << " root..  LogBound= " << LogBound << endl;
    }
    return LogBound;
  }


  BigRat RootBound2(ConstRefRingElem f, long NumGraeffeIters)
  {
    RootBound_CheckArg(f, CoCoA_ERROR_CONTEXT);
    if (NumGraeffeIters < -1 || NumGraeffeIters > 25)
      CoCoA_THROW_ERROR2(ERR::OutOfRange, "number of iters must be in 0..25");

    const double LogBound = LogRootBound(f, NumGraeffeIters);
    if (LogBound == LogZero)  return BigRat(0);
    return ApproxExp(LogBound);
  }


  // Input is f(x) = a_0 + a_1*x + ... + a_d*x^d with deg >= 1
  // Output is -|a_0| - |a_1|*x - ... -|a_{d-1}|*x^{d-1} + |a_d|*x^d
  // Output might also be scaled(???) or just call prim(ans)???
  RingElem RootBoundTransform(ConstRefRingElem f)
  {
    static const ErrorInfo ErrMesg(ERR::BadConvert, "RootBoundTransform");
    RootBound_CheckArg(f, CoCoA_ERROR_CONTEXT);
    const ring& P = owner(f);
    RingElem ans = monomial(P, abs(ConvertTo<BigRat>(LC(f),ErrMesg)), LPP(f));
    bool FirstTime = true;
    for (SparsePolyIter it=BeginIter(f); !IsEnded(it); ++it)
    {
      if (FirstTime)  { FirstTime = false; continue; }
      ans += monomial(P, -abs(ConvertTo<BigRat>(coeff(it),ErrMesg)), PP(it));
    }
    return ans;
  }
  

} // end of namespace CoCoA
