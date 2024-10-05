#ifdef CoCoA_WITH_BOOST_temporarily_disabled // needs boost 1.41 or later

//   Copyright (c)  2013  Bruno Simoes

//   This file is part of the source of CoCoALib, the CoCoA Library.

//   CoCoALib is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.

//   CoCoALib is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.

//   You should have received a copy of the GNU General Public License
//   along with CoCoALib.  If not, see <http://www.gnu.org/licenses/>.

#include "CoCoA/library.H"
//#include "CoCoA/RBGWV.H"

#include <iostream>
#include <vector>
#include <iomanip>

#ifdef CoCoA_WITH_ZToolkit
	#include "ZToolkit.C"
#endif

#define TEST_ASSERT(cond) do { if (!(cond)) CoCoA::AssertionFailed(#cond, __FILE__, __LINE__); } while(0)

using namespace CoCoA;
				
std::vector<SBStats> myStatsList;

// Setup table headers
void SetupHeadersFn(CoCoA::BenchmarkToolkit &toolkit)
{
	toolkit.header.push_back("Algorithm"); 
	toolkit.header.push_back("Ring"); 
	toolkit.header.push_back("Char"); 
	toolkit.header.push_back("Ordering"); 
	toolkit.header.push_back("Nr. Generators"); 
	toolkit.header.push_back("Timings");			
	toolkit.header.push_back("Pairs Created");	
		
	// Create tables (title, array index)
	toolkit.tableKV.push_back(std::make_pair("Number of generators", 4));
	toolkit.tableKV.push_back(std::make_pair("Timings (ms)", 5));
	toolkit.tableKV.push_back(std::make_pair("Pairs Created", 6));	
}

// Function used to add values to the table
std::vector<std::string> UpdateStatsFn( CoCoA::BenchmarkExample &example, 
										const std::vector<RingElem> &gbasis, 
										const double time, 
										const int k, // test case
										const std::string algorithm)
{

	SBStats myStats = myStatsList[myStatsList.size()-1];
	std::vector<std::string> stat;
	std::stringstream ss;

	// Algorithm
	stat.push_back(example.id);
	
	// Ring
	stat.push_back(example.TestRings[k].name);

	// Char
	ss.str("");	
	ss << example.TestRings[k].characteristic;
	stat.push_back(ss.str());

	// Ordering
	stat.push_back(example.TestRings[k].ordering);
	
	// Number of generators
	ss.str("");	
	ss << gbasis.size();	
	stat.push_back(ss.str());
	
	// Timings
	if( time == -1 )
	{
		stat.push_back("-");
	}
	else 
	{
		ss.str("");	
		ss << std::setprecision(2) << time;
		stat.push_back(ss.str());
	}
	
	if(gbasis.size())
	{
		if(algorithm == "RBGWV")
		{
			ss.str("");	
			ss << myStats.myStatsNewPairs;
			stat.push_back(ss.str());				
		}
		else
		{
			ss.str("");	
			ss << myStats.myPInserted;
			stat.push_back(ss.str());	
		}
	}
	else
	{
		stat.push_back("-");
	}
	
	ss.str("");	
	return stat;				
}

// Function wrapper for GBasis
void myGBasis(const ideal& in, std::vector<RingElem>& output)
{
	std::vector<RingElem> generators = gens(in);
	SBStats stats(generators.size(), GReductor::ourDefaultStatLevel);	
	ComputeGBasis( output, stats, generators, GReductor::ourDefaultStatLevel);
	myStatsList.push_back(stats);
}

// Function wrapper for SATGBasis
void mySATGBasis(const ideal& in, std::vector<RingElem>& output)
{
	std::vector<RingElem> generators = gens(in);
	SBStats stats(generators.size(), GReductor::ourDefaultStatLevel);
	ComputeSATGBasis( output, stats, generators, GReductor::ourDefaultStatLevel);
	myStatsList.push_back(stats);
}

// Function wrapper for SATMixGBasis
void mySATMixGBasis(const ideal& in, std::vector<RingElem>& output)
{
	std::vector<RingElem> generators = gens(in);
	SBStats stats(generators.size(), GReductor::ourDefaultStatLevel);
	ComputeSATMixGBasis( output, stats, generators, GReductor::ourDefaultStatLevel);
	myStatsList.push_back(stats);
}

#ifdef CoCoA_WITH_THRUST
	void myRBGWVGBasis( const ideal& in, std::vector<RingElem>& theGB )
	{
		std::vector<RingElem> generators = gens(in);
		SBStats stats(generators.size(), GReductor::ourDefaultStatLevel);	
		SparsePolyRing SPR(owner(generators));
		GRingInfo GRI(SPR,true,false,NewDivMaskEvenPowers());
		RBGWV g(GRI, generators);	
		g.computeGBasis(theGB);
		myStatsList.push_back(g.myStats);
	}
#endif

int main()
{
	try
	{
		GlobalManager CoCoAFoundations;

		BenchmarkToolkit db;
		
		db.init("ex-BenchmarkToolkit-xml.in");

		// Declare algorithms label
		const std::string LabelGBasis("GBasis");
		const std::string LabelSATGBasis("SAT GBasis");
		const std::string LabelSATMixGBasis("SAT Mix GBasis");

#ifdef CoCoA_WITH_THRUST
		const std::string LabelRBGWVGBasis("RBGWV");
#endif

		// Add functions to the benchmark list	
		std::map<std::string, MyTestFn> BenchmarkFns;
		BenchmarkFns[LabelGBasis] = &myGBasis;
   		BenchmarkFns[LabelSATGBasis] = &mySATGBasis;
   		BenchmarkFns[LabelSATMixGBasis] = &mySATMixGBasis;

#ifdef CoCoA_WITH_THRUST
   		BenchmarkFns[LabelRBGWVGBasis] = &myRBGWVGBasis;   		
#endif
   		
		// Create a filter to select specific tags
   		std::vector<std::string> filterTags;
   		
  		filterTags.push_back("katsura7");
   		filterTags.push_back("chandra6");
   		filterTags.push_back("butcher");
   		filterTags.push_back("chandra4");	
   		
   		// Create a filter to select specific rings
   		std::vector<std::string> filterRings;
		filterRings.push_back("NewZZmod-32003");
		
		// Create a filter to select specific orderings
		// An empty array should be passed to accept all
   		std::vector<std::string> filterOrds;   		
   		filterOrds.push_back("DegRevLex");
	
   		// Setup table headers
		SetupHeadersFn(db);

   		std::cout << std::endl << "Computing ... " << std::endl;

   		// Execute 
		db.execute(BenchmarkFns, &UpdateStatsFn, filterTags, filterRings, filterOrds, 1);

		// Export to LaTeX
   		std::cout << std::endl << std::endl << "Exporting LaTeX file ... ";
		db.ExportAsLatex();

		// Print to ctdout		
   		std::cout << std::endl << "Printing to PlainText ... " << std::endl;
		db.ExportAsPlainText(); 
				
		return 0;
	}
	catch (const CoCoA::ErrorInfo& err)
	{
		std::cerr << "***ERROR***  UNCAUGHT CoCoA Error";
		ANNOUNCE(std::cerr, err);
	}
	catch (const std::exception& exc)
	{
		std::cerr << "***ERROR***  UNCAUGHT std::exception: " << exc.what() << std::endl;
	}
	catch(...)
	{
		std::cerr << "***ERROR***  UNCAUGHT UNKNOWN EXCEPTION" << std::endl;
	}

	BuildInfo::PrintAll(std::cerr);
	return 1;
}
#endif

