// Copyright (c) 2016  John Abbott,  Anna Bigatti
// Orig authors: 2016  Anna Maria Bigatti, Eduardo Saenz-de-Cabezon
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

#include <map>     // std::multimap

using namespace std;
//using back_inserter;

//----------------------------------------------------------------------
const string ShortDescription =
  "MVT for square free monomial ideals  \n";

const string LongDescription =
  "MVT for square free monomial ideals  \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  bool DEBUG = false;

  typedef std::vector<DynamicBitset> PPVectorSqFr;


  class MVTSquareFree
  {
  public:
    explicit MVTSquareFree(const ideal& I);
    void myComputeBetti(std::vector<int>& bettis);
    std::ostream& myOutput(std::ostream& out) const;

  private:
    explicit MVTSquareFree(const PPVectorSqFr& I, int n);
    void myComputeRec(PPVectorSqFr& I, const DynamicBitset& AccPivot, const int dim, const BigInt& pos);
    void myComputeRec2(PPVectorSqFr& PrimeI, PPVectorSqFr& EdgeI, const DynamicBitset& AccPivot, const int dim, const BigInt& pos_orig, const MVTSquareFree& EdgeMVT);
    void myInsert(const DynamicBitset&, const int&, const BigInt&);
    void myInsertGens(const PPVectorSqFr& I, const DynamicBitset& AccPivot, const int dim, const BigInt& pos);
    void myInsertTree(const MVTSquareFree& EdgeMVT, const DynamicBitset& AccPivot, const int dim, const BigInt& pos);
    typedef std::multimap<DynamicBitset, BigInt> MultimapDynBS_BigInt;
    
  private: // data members
    std::vector<MultimapDynBS_BigInt> myDimVec;
    bool IamEmpty;
  };


  //*******************************************************//
  void InterreduceSort(PPVectorSqFr& g);
  void tilde(PPVectorSqFr& Prime_g, PPVectorSqFr& Edge_g, DynamicBitset& AccPivot, PPVectorSqFr& I);
  

  //********************************************************//
  // FUNCTIONS and MEMBER FUNCTIONS for MVTSquareFree class //
  //********************************************************//
  
  std::ostream& operator<<(std::ostream& out, const MVTSquareFree& M)
  { return M.myOutput(out); }

  //-- public member functions -------------------------------------

  //-- public constructor
  MVTSquareFree::MVTSquareFree(const ideal& I):
    myDimVec(len(gens(I)))
  {
    VerboseLog VERBOSE("MVTSquareFree(ideal)");
    IamEmpty = true;
    PPMonoid PPMI = PPM(RingOf(I));
    const std::vector<RingElem>& gensI = gens(I);
    PPVectorSqFr J;
    J.reserve(len(gensI));
    for (long i=0; i<len(gensI); ++i)
    {
      // ANNA: should be a CoCoA_ASSERT
      VERBOSE(9) << gensI[i] << endl;
      if (!IsMonomial(gensI[i])) CoCoA_THROW_ERROR("is not monomial","MVTSquareFree");
      if (!IsOne(LC(gensI[i]))) CoCoA_THROW_ERROR("coeff is not 1","MVTSquareFree");
      if (deg(LPP(gensI[i]))>2) CoCoA_THROW_ERROR("deg>2","MVTSquareFree");
      if (!IsSqFree(LPP(gensI[i]))) CoCoA_THROW_ERROR("not squarefree","MVTSquareFree");
      J.push_back(DynamicBitset(LPP(gensI[i])));
    }
    InterreduceSort(J);
    DynamicBitset AccumulatingPivot(NumIndets(PPMI));
    BigInt pos(1);
    myInsertGens(J, AccumulatingPivot, 0, pos);
    myComputeRec(J, AccumulatingPivot, 0, pos);
  }


  std::ostream& MVTSquareFree::myOutput(std::ostream& out) const
  {
    if (IamEmpty) out << " [empty]"; else out << endl;
    for (long i=0; i<len(myDimVec); ++i)
    {
      const MultimapDynBS_BigInt& m=myDimVec[i];
      if (!m.empty())
      {
        out << "tree in dim = " << i << endl;
        for (MultimapDynBS_BigInt::const_iterator it=m.begin();it!=m.end();++it)
        {
          out<<"  Multidegree "<< it->first << " (deg="<<count(it->first)<<"):";
          out<<"  pos "<< it->second << endl;
        }
      }
    }
    return out;
  }


  void MVTSquareFree::myComputeBetti(vector<int>& bettis)
  {
    if (len(bettis) < len(myDimVec))
    {
      for (long i=len(bettis);  i<len(myDimVec); ++i)
        if ((myDimVec[i]).size() != 0)
          CoCoA_THROW_ERROR("impossible Betti number","myComputeBetti");
      for (long i=0; i<len(bettis); ++i)   bettis[i] = (myDimVec[i]).size();
    }
    else // len(bettis) >= len(myDimVec)
    {
      for (long i=0; i<len(myDimVec); ++i) bettis[i] = (myDimVec[i]).size();
      for (long i=len(myDimVec); i<len(bettis); ++i) bettis[i] = 0;
    }
  }


  //-- private member functions -------------------------------------

  //-- private constructor (intermediate subtree)
  MVTSquareFree::MVTSquareFree(const PPVectorSqFr& I, int n):
    myDimVec(n)
  {
    IamEmpty = true;
    if (I.empty()) return;
    //if (DEBUG) cout << " ------MVT-subtree-begin --\\" << endl;
    PPVectorSqFr copyI;
    copyI.assign(I.begin(), I.end());
    DynamicBitset AccumulatingPivot(n);
    BigInt pos(1);
    // --DON'T INSERT GENERATORS--
    myComputeRec(copyI, AccumulatingPivot, 0, pos);
    //if (DEBUG) cout << " ------MVT-subtree-end ----//" << endl;
  }
  

  void MVTSquareFree::myInsert(const DynamicBitset& multidegree, const int& dim, const BigInt& pos)
  {
    IamEmpty = false;
    (myDimVec[dim]).insert(std::make_pair(multidegree, pos));
  }    


  void MVTSquareFree::myInsertGens(const PPVectorSqFr& I, const DynamicBitset& AccPivot, const int dim, const BigInt& pos)
  {
    for (PPVectorSqFr::const_iterator it=I.begin(); it!=I.end(); ++it)
      myInsert((*it)|AccPivot, dim, pos);
  }
  

  void MVTSquareFree::myInsertTree(const MVTSquareFree& EdgeMVT, const DynamicBitset& AccPivot, const int dim, const BigInt& pos)
  {
    if (DEBUG) cout << "InsertTree: EdgeMVT is" << EdgeMVT << endl
                    << "  AccPivot=" << AccPivot
                    << "  dim=" << dim
                    << "  pos=" << pos << endl;
    for (long i=0; i<len(EdgeMVT.myDimVec); ++i)
    {
      const MultimapDynBS_BigInt& DimMap = EdgeMVT.myDimVec[i];
      for (MultimapDynBS_BigInt::const_iterator it2=DimMap.begin(); it2!=DimMap.end(); ++it2)
        myInsert((it2->first)|AccPivot, i+dim, (it2->second)*pos);
    }
    if (DEBUG) cout << "-->InsertTree: tree now is";
    if (DEBUG) myOutput(cout);
  }


  void MVTSquareFree::myComputeRec2(PPVectorSqFr& PrimeI, PPVectorSqFr& EdgeI, const DynamicBitset& AccPivot, const int dim, const BigInt& pos_orig, const MVTSquareFree& EdgeMVT)
  {
    BigInt pos = pos_orig;
    myInsertGens(PrimeI, AccPivot, dim, pos);
    myInsertGens(EdgeI,  AccPivot, dim, pos);
    while (!PrimeI.empty())
    {
      DynamicBitset pivot = PrimeI.back();
      PrimeI.pop_back();
      // vertical branch
      DynamicBitset TildeAccPivot = AccPivot|pivot;
      PPVectorSqFr TmpPrimeI;
      TmpPrimeI.assign(PrimeI.begin(), PrimeI.end());
      myComputeRec2(TmpPrimeI, EdgeI, TildeAccPivot, 1+dim, 2*pos, EdgeMVT);
      // horizontal branch
      pos = 2*pos+1;
    }
    if (!EdgeMVT.IamEmpty) myInsertTree(EdgeMVT, AccPivot, dim, pos);
  }
  

  void MVTSquareFree::myComputeRec(PPVectorSqFr& I, const DynamicBitset& AccPivot, const int dim, const BigInt& pos_orig)
  {
    BigInt pos = pos_orig;
    if (dim!=0) myInsertGens(I, AccPivot, dim, pos);
    while (len(I)>2)
    {
      //      cout << " -- MVT     pos = " << pos   << "  I = " << I << endl;
      PPVectorSqFr PrimeTildeI, EdgeTildeI;
      DynamicBitset TildeAccPivot = AccPivot;
      tilde(PrimeTildeI, EdgeTildeI, TildeAccPivot, I); // modifies I
      // vertical branch
      if (PrimeTildeI.empty())
        myComputeRec(EdgeTildeI, TildeAccPivot, 1+dim, 2*pos);
      else 
      {
//         if (len(PrimeTildeI)==1 && len(EdgeTildeI)<2)
//         {
//           myInsert(PrimeTildeI[0]|TildeAccPivot, 1+dim, 2*pos);
//           if (!EdgeTildeI.empty())
//           {
//             myInsert(EdgeTildeI[0]|TildeAccPivot, 1+dim, 2*pos);
//             myInsert(PrimeTildeI[0]|EdgeTildeI[0]|TildeAccPivot, 2+dim, 4*pos);
//           }
//         }
//         else
        {
          MVTSquareFree TmpMdM(EdgeTildeI, len(TildeAccPivot));
          myComputeRec2(PrimeTildeI, EdgeTildeI, TildeAccPivot, 1+dim, 2*pos, TmpMdM);
        }
      }
      // horizontal branch
      pos = 2*pos+1;
    }  // while
    if (len(I)==2) myInsert(I[0]|I[1]|AccPivot, 1+dim, 2*pos);
  }

  //********************************************************//
  //********************************************************//

  //  namespace // anonymous for functions local to this file/compilation unit.
  //  {
   
    bool DBLessThan(const DynamicBitset& f, const DynamicBitset& g)
    { return f<g; }


    bool divides00(const PPVectorSqFr& g, DynamicBitset b)
    {
      for (PPVectorSqFr::const_iterator it=g.begin(); it!=g.end(); ++it)
        if (IsSubset(*it,b)) return true;
      return false;
    }
      

    void InterreduceSort(PPVectorSqFr& g)
    {
      PPVectorSqFr newg;
      newg.reserve(len(g));
      sort(g.begin(), g.end(), DBLessThan);
      for (PPVectorSqFr::const_iterator it=g.begin(); it!=g.end(); ++it)
        if (!divides00(newg, *it)) newg.push_back(*it);
      swap(newg,g);
    }


//     void tilde1(PPVectorSqFr& Prime_g, PPVectorSqFr& Edge_g, const PPVectorSqFr& I, const DynamicBitset& pivot)
//     {
//       long count=0;
//       for (long i=0; i<len(I); ++i)
//         if (IsDisjoint(I[i], pivot))
//           g.push_back(I[i]);
//         else
//         {
//           g.push_back(I[i]-pivot);
//           ++count;
//         }
//       if (count!=0)
//         InterreduceSort(g);
//     }
  

    void tilde2(PPVectorSqFr& Prime_g, PPVectorSqFr& Edge_g, const PPVectorSqFr& I, const DynamicBitset& pivot)
    {
      PPVectorSqFr Edge_reduced;
      for (long i=0; i<len(I); ++i)
        if (IsDisjoint(I[i], pivot))
          Edge_g.push_back(I[i]);
        else
          Prime_g.push_back(I[i]-pivot);
      InterreduceSort(Prime_g);
      for (long i=0; i<len(Edge_g); ++i)
        if (!divides00(Prime_g, Edge_g[i]))
          Edge_reduced.push_back(Edge_g[i]);
      swap(Edge_g, Edge_reduced);
    }


    void tilde(PPVectorSqFr& Prime_g, PPVectorSqFr& Edge_g, DynamicBitset& AccPivot, PPVectorSqFr& I)
    {
      if (!Prime_g.empty()) CoCoA_THROW_ERROR("Prime_g: not empty", "tilde");
      if (!Edge_g.empty())  CoCoA_THROW_ERROR("Edge_g: not empty", "tilde");
      Prime_g.reserve(len(AccPivot));
      Edge_g.reserve(len(I));
      DynamicBitset pivot = I.back();
      I.pop_back();
      if (DEBUG) cout << "tilde: pivot = " << pivot << endl;
      AccPivot |= pivot;
      //      if (len(I) <= 20)  tilde1(Prime_g, Edge_g, I, pivot);  else
      tilde2(Prime_g, Edge_g, I, pivot);
      if (DEBUG) cout << "tilde: Prime_g = " << Prime_g << endl;
      if (DEBUG) cout << "tilde: Edge_g  = " << Edge_g << endl;
    }
  

    //     bool HasHomology(const PPWithMask& pp, const PPVector& J)
    //     {
    //       PPMonoid M=owner(PP(pp));
    //       int n=NumIndets(M);
    //       PPMonoidElem ppi=PP(pp)/product(indets(M));
    //       if (IsDivisible(ppi, J))
    //         return false;
    //       else
    //         for (int i=0; i<n; ++i)
    //           if (!IsDivisible(indet(M,i)*ppi, J))
    //             return false;
    //       return true;
    //     }

    // bool HasHomology_gen(const PPWithMask& pp, const PPVector& J, const PPMonoidElem& m)
    //     {
    //       PPMonoid M=owner(PP(pp));
    //       int n=NumIndets(M);
    //       PPMonoidElem ppi=(PP(pp)*m)/product(indets(M));
    //       if (IsDivisible(ppi, J))
    //         return false;
    //       else
    //         for (int i=0; i<n; ++i)
    // 	if(exponent(m,i)==0)
    // 	{
    //           if (!IsDivisible(indet(M,i)*ppi, J))
    //             {return false;}
    // 	}
    // 	 return true;
    //     }

  //  }  // end of anonymous namespace



  //////////////////////////////////////////////////////////////////////

  ideal ConsecutiveLinear(long n, long k)
  {
    SparsePolyRing P = NewPolyRing(RingQQ(), SymbolRange("x",0,n-1));
    PPMonoidElem pp(PPM(P));
    std::vector<RingElem> g;
    for (long i=0; i<n-k+1; ++i)
    {
      pp = indet(PPM(P), i);;
      for (long j=1; j<k; ++j)  pp *= indet(PPM(P), i+j);
      g.push_back(monomial(P,pp));
    }
    return ideal(g);
  }


  ideal wheel(long n)
  {
    VerboseLog VERBOSE("wheel");
    SparsePolyRing P = NewPolyRing(RingQQ(), SymbolRange("x",0,n-1));
    PPMonoidElem pp(PPM(P));
    std::vector<RingElem> g;
    for (long i=1; i<n-1; ++i)
    {
      g.push_back(monomial(P, indet(PPM(P), i) * indet(PPM(P), 0)));
      VERBOSE(9) << g.back() << "   " << flush;
      g.push_back(monomial(P, indet(PPM(P), i) * indet(PPM(P), i+1)));
      VERBOSE(9) << g.back() << endl;
    }
    g.push_back(monomial(P, indet(PPM(P), n-1) * indet(PPM(P), 0)));
    VERBOSE(9) << g.back() << "   " << flush;
    g.push_back(monomial(P, indet(PPM(P), n-1) * indet(PPM(P), 1)));
    VERBOSE(9) << g.back() << endl;

    return ideal(g);
  }


  ideal star(long n)
  {
    SparsePolyRing P = NewPolyRing(RingQQ(), SymbolRange("x",0,n-1));
    PPMonoidElem pp(PPM(P));
    std::vector<RingElem> g;
    for (long i=1; i<n; ++i)
      g.push_back(monomial(P, indet(PPM(P), i) * indet(PPM(P), 0)));
    return ideal(g);
  }


  //////////////////////////////////////////////////////////////////////
  void test_program(const ideal& I)
  {
    double t0;
    SparsePolyRing P = RingOf(I);
    DivMaskRule DMR = NewDivMaskSingleBitWrap();
    VerboseLog VERBOSE("test_program");

    cout << I << endl << endl;

    VERBOSE(5) << "MVTSquareFree  " << endl;
    t0 = CpuTime();
    MVTSquareFree MVTSqFr(I);
    cout << "Construction of MV Tree SqFr: sec " << CpuTime() - t0 << endl;

    VERBOSE(5) << "general MVT  " << endl;
    t0 = CpuTime();
    MultidegreeMap MVT;
    PPVector PPs(PPM(P), DMR);
    VERBOSE(9) << " -before convert" << endl;
    convert(PPs, gens(I));
    VERBOSE(9) << " -before MayerVietorisTree" << endl;
    MayerVietorisTree(MVT, PPs);
    //ReducedMayerVietorisTree(MVT, PPs);
    cout << "Construction of MV Tree:      sec " << CpuTime() - t0 << endl;

    if (DEBUG) cout << " --MVTSqFr is" << MVTSqFr << endl;
    if (DEBUG){cout << " --MVT is" << endl; PrintMultidegreeMap(MVT);}

    VERBOSE(5) << "--==compare-bettis==--" << endl;
    std::vector<int> bettis2(NumIndets(P),0), bettis3(NumIndets(P),0);

    VERBOSE(5) << "MVTSquareFree  " << endl;
    t0 = CpuTime();
    VERBOSE(9) << " -before myComputeBetti" << endl;
    MVTSqFr.myComputeBetti(bettis2);
    cout << "Reading the Betti numbers SqFr: sec " << CpuTime() - t0 << endl;
    OutputRange(cout, bettis2.begin(), bettis2.end());  cout << endl;

    VERBOSE(5) << "general MVT  " << endl;
    t0 = CpuTime();
    Bettis(bettis3, MVT);
    cout << "Reading the Betti numbers:      sec " << CpuTime() - t0 << endl;
    OutputRange(cout, bettis3.begin(), bettis3.end());  cout << endl;

    if (len(bettis3) != len(bettis2)) cout << " -different len!!";
    for (long i=0; i<len(bettis2); ++i)
      if (bettis2[i]!=bettis3[i]) cout << " -different(" << i << ")";
    
    cout << " ----" << endl;
  }
  
  //////////////////////////////////////////////////////////////////////
  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;
    //    cout << boolalpha; // so that bools print out as true/false

    //\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
    //\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
    //\/  ***PUT YOUR CODE HERE***  \/\/

    SetVerbosityLevel(5);
    
    DynamicBitset::ourOutputStyle = DynamicBitset::WithSeparators;

    if (false)
    {
      cout << " ----small-------------------------------------" << endl;
      long N = 9;
      SparsePolyRing P = NewPolyRing(RingQQ(), SymbolRange("x",0,N-1));
    
      std::vector<RingElem> f;
      f.push_back(RingElem(P,"x[1]*x[2]"));
      f.push_back(RingElem(P,"x[3]*x[2]"));
      f.push_back(RingElem(P,"x[4]*x[2]"));

      test_program(ideal(f));
    }    

    //    if (true)
    if (false)
    {
      cout << " ---ConsecutiveLinear--------------------------" << endl;
      long N = 28;
      DEBUG = false;
      test_program(ConsecutiveLinear(N, 2));
    }

    if (true)
    {
      cout << " ---wheel--------------------------------------" << endl;
      long N = 5;
      DEBUG = false;
      test_program(wheel(N));
    }

    //    if (true)
    if (false)
    {
      cout << " ----star------------------" << endl;
      long N = 22;
      DEBUG = false;
      test_program(star(N));
    }

    if (false)
    {
      cout << " ----Test-for-subtree-recursion2------------------" << endl;
      long N = 10;
      DEBUG = false;
      SparsePolyRing P = NewPolyRing(RingQQ(), SymbolRange("x",0,N-1));
    
      std::vector<RingElem> f;
      for (long i=0; i<N-2; ++i) f.push_back(indet(P,i)*indet(P,N-2));
      f.push_back(indet(P,0)*indet(P,N-1));
      test_program(ideal(f));
    }

    if (false)
    {
      cout << " ---indets--------------------------------------" << endl;
      long N = 22;
      DEBUG = false;
      SparsePolyRing P = NewPolyRing(RingQQ(), SymbolRange("x",0,N-1));
      test_program(ideal(indets(P)));
    }    

    //\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
    //\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/

  }


} // end of namespace CoCoA

//----------------------------------------------------------------------
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-MVT-edge.C,v 1.11 2022/03/04 11:29:57 bigatti Exp $
// $Log: ex-MVT-edge.C,v $
// Revision 1.11  2022/03/04 11:29:57  bigatti
// Summary: added dates for Orig author
//
// Revision 1.10  2022/02/13 09:56:57  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.9  2020/06/17 15:49:19  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.8  2017/07/24 14:55:37  bigatti
// -- added verbosity
//
// Revision 1.7  2017/04/27 15:24:52  bigatti
// -- changed ReadExpr --> RingElem
//
// Revision 1.6  2016/06/06 11:45:02  bigatti
// -- aligned code with other MVT code examples
//
// Revision 1.5  2016/05/27 15:52:52  bigatti
// -- new implementation using class
//
// Revision 1.4  2016/05/25 15:19:41  bigatti
// -- new data type for tree
//
// Revision 1.3  2016/05/25 13:29:15  bigatti
// -- now it works (slow, but works)
//
// Revision 1.2  2016/05/24 12:34:23  bigatti
// Summary: changed names for redeclaration (temporary, for deleting anonymous namespace)
//
// Revision 1.1  2016/05/24 12:22:39  bigatti
// -- first import
