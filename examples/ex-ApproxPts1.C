// Copyright (c) 2006  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "A short example showing how to use the approximate point preprocessing algorithms.\n"
  "See the file ex-ApproxPt1.in for a sample input.\n";

const string LongDescription =
  "A short example which reads a set of approximate points (with common epsilon)\n"
  "and then applies the various preprocessing algorithms to them.  It prints out\n"
  "the resulting preprocessed set with the corresponding weights, and the time taken.";
//----------------------------------------------------------------------


namespace CoCoA
{

  RingElem ConvertToRingElem(const ring& R, double z)
  {
    return RingElem(R, ConvertTo<BigRat>(z));
  }

  vector<RingElem> ConvertToVectorRingElem(const ring& R, const vector<double>& v)
  {
    const long n = len(v);
    vector<RingElem> ans; ans.reserve(n);
    for (long i=0; i < n; ++i)
      ans.push_back(ConvertToRingElem(R, v[i]));
    return ans;
  }

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;
  
    cout << "Insert the space dimension: " << endl;
    long dim;
    cin >> dim;
    if ( !cin )
      CoCoA_THROW_ERROR("Input must be a positive integer", "main program in ex-ApproxPts1");
    if (dim < 1 || dim > 1000000)
      CoCoA_THROW_ERROR("Ridiculous input dimension", "main program in ex-ApproxPts1");
  
    ring QQ = RingQQ();
    vector<double> InputEps(dim);
    cout << "Insert the tolerance in each dimension: " << endl;
    for (long i=0; i < dim; ++i)
    {
      cin >> InputEps[i];
      if (!cin || InputEps[i] < 0) { CoCoA_THROW_ERROR("bad input", "main program in ex-ApproxPts1"); }
    }
    const vector<RingElem> epsilon = ConvertToVectorRingElem(QQ, InputEps);

    cout << "Insert the number of points to be preprocessed: " << endl;
    long NumPts;
    cin >> NumPts;
    if (NumPts < 1 || NumPts > 1000000)
      CoCoA_THROW_ERROR("Ridiculous number of points", "main program in ex-ApproxPts1");

    vector<ApproxPts::PointR> OrigPts;  OrigPts.reserve(NumPts);
    cout << "Insert the coordinates of the points " << endl;
    for (long i=0; i < NumPts; ++i)
    {
      vector<double> InputPt(dim);
      for (long j=0; j < dim; ++j)
      {
        cin >> InputPt[j];
        if (!cin) { CoCoA_THROW_ERROR("bad input", "main program in ex-ApproxPts1"); }
      }
      OrigPts.push_back(ConvertToVectorRingElem(QQ, InputPt));
    }

    cout << endl
         << "Read " << len(OrigPts) << " original points." << endl
         << endl;

    double StartTime, EndTime;
    cout << "-------------------------------------------------------" << endl;
    vector<ApproxPts::PointR> PreprocessedPts;
    vector<long> weights;
    StartTime = CpuTime();
    PreprocessPtsGrid(PreprocessedPts, weights, OrigPts, epsilon);
    EndTime = CpuTime();
    cout << "Grid algm produces " << len(PreprocessedPts) << " preprocessed points." << endl
         << endl
         << "The preprocessed points are: " << PreprocessedPts << endl
         << "and their weights are: " << weights << endl
         << endl
         << "CPU time for grid algm: " << EndTime-StartTime << endl;

    cout << "-------------------------------------------------------" << endl;
    StartTime = CpuTime();
    PreprocessPtsAggr(PreprocessedPts, weights, OrigPts, epsilon);
    EndTime = CpuTime();
    cout << "Aggr algm produces " << len(PreprocessedPts) << " preprocessed points." << endl
         << endl
         << "The preprocessed points are: " << PreprocessedPts << endl
         << "and their weights are: " << weights << endl
         << endl
         << "CPU time for aggr algm: " << EndTime-StartTime << endl;

    cout << "-------------------------------------------------------" << endl;
    StartTime = CpuTime();
    PreprocessPtsSubdiv(PreprocessedPts, weights, OrigPts, epsilon);
    EndTime = CpuTime();
    cout << "Subdiv algm produces " << len(PreprocessedPts) << " preprocessed points." << endl
         << endl
         << "The preprocessed points are: " << PreprocessedPts << endl
         << "and their weights are: " << weights << endl
         << endl
         << "CPU time for subdiv algm: " << EndTime-StartTime << endl;
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
