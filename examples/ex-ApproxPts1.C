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

//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-ApproxPts1.C,v 1.17 2022/02/13 09:56:55 abbott Exp $
// $Log: ex-ApproxPts1.C,v $
// Revision 1.17  2022/02/13 09:56:55  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.16  2020/06/17 15:49:19  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.15  2015/06/25 16:40:27  abbott
// Summary: Moved code into namespace CoCoA (redmine 739)
// Author: JAA
//
// Revision 1.14  2014/05/10 20:46:08  abbott
// Summary: 3 very minor improvements
// Author: JAA
//
// Revision 1.13  2013/03/27 18:24:33  abbott
// Added approx point preprocessing to C5; also changed names of the fns, and updated doc.
//
// Revision 1.12  2012/09/21 13:38:12  abbott
// Replaced t0 by StartTime and EndTime.
//
// Revision 1.11  2012/02/08 17:40:23  bigatti
// -- changed: Z,Q -> ZZ,QQ
//
// Revision 1.10  2011/12/23 15:48:54  bigatti
// -- fixed typo
//
// Revision 1.9  2011/12/13 14:20:35  abbott
// Cleaned code following Caleo's comments.
// Added LongDescription.
//
// Revision 1.8  2011/08/24 10:46:32  bigatti
// -- renamed QQ --> BigRat
// -- sorted #include
//
// Revision 1.7  2011/03/11 12:04:00  bigatti
// -- changed size_t --> long
// -- changed size --> len
//
// Revision 1.6  2010/12/26 13:02:03  abbott
// Finished changing "GlobalXXXput()" into the corresponding standard C++ stream.
//
// Revision 1.5  2010/12/17 16:07:55  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.4  2009/12/11 11:46:32  abbott
// Changed fn  convert  into  IsConvertible.
// Added template procedure  convert.
// New version because change is not backward compatible.
//
// Revision 1.3  2008/11/23 18:23:22  abbott
// Now uses points with coords in RingQ rather than represedented as C++ doubles.
//
// Revision 1.2  2008/11/21 21:12:19  abbott
// Update inline with the changes to the preprocessing code.
// Changed semiwidth into tolerance.
// Program now prints out the weights too.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.5  2007/03/03 14:15:45  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.4  2007/02/12 15:29:07  bigatti
// -- added strings ShortDescription and LongDescription for indexing
//
// Revision 1.3  2007/02/10 18:44:04  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.2  2006/06/21 17:05:47  cocoa
// Major overhaul of approx point preprocessing algms.
//
// Revision 1.1.1.1  2006/05/30 11:39:36  cocoa
// Imported files
//
// Revision 1.3  2006/05/29 16:17:10  cocoa
// Changed a variable name.
//
// Revision 1.2  2006/05/22 15:52:16  cocoa
// Added preprocess-disg algorithm to ApproxPts.
// Sundry minor improvements.
//
// Revision 1.1  2006/05/12 13:16:30  cocoa
// Added functions for preprocessing approximate points.
//
