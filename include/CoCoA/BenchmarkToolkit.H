#ifdef CoCoA_WITH_BOOST_temporarily_disabled // needs boost 1.41 or later
#ifndef CoCoA_BENCHMARKTOOLKIT_H
#define CoCoA_BENCHMARKTOOLKIT_H

//   Copyright (c)  2013  Bruno Simoes

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

//#include "CoCoA/TmpGRStats.H"

#include <map>
#include <string>
#include <vector>
#include <iostream>

#include <boost/date_time/posix_time/ptime.hpp>
#include <boost/date_time/microsec_time_clock.hpp>

namespace CoCoA
{

	class RingElem;     // fwd decl -- defined in BigRat.H
	class ideal;        // fwd decl -- defined in ideal.H
	class ring;         // fwd decl -- defined in BigRat.H

	struct PolyRingInfo
	{
		std::string id;
		std::string name;
		int characteristic;
		std::string ordering;
		std::vector<CoCoA::ring> ring;
	};	
		
	struct BenchmarkExample
	{
		std::string id;
		std::string name;
		std::string description;
		std::string references;
		std::string notes;

		std::vector<std::string> 	TagsList;	
		std::vector<PolyRingInfo>	TestRings;
		std::vector<ideal> 			TestIdeals;
	};
	
	struct BenchmarkInfo
	{
		// Groebner basis algorithm being tested.
		std::string algorithm;
	
		// Statistics
		std::vector<std::vector<std::string> > stats;	
	};

	typedef void (*MyTestFn)(const CoCoA::ideal &, std::vector<RingElem> &);
	typedef std::vector<std::string> (*UpdateStatsMyTestFn)(  CoCoA::BenchmarkExample &,
									  						const std::vector<RingElem> &, 
															const double, 
															const int,
															const std::string);

	struct BenchmarkToolkit
	{

		// Contains all the information about the benchmark examples
		std::vector<BenchmarkExample> database;

		// Map between name and ring implementation
		std::map<std::string, std::vector<PolyRingInfo> > PrMap;
		
		// Pointer to benchmark examples used
		std::vector<std::pair<int,int> > testset;

		// Table headers.
		std::vector< std::string > header;

		// Table caption and associated benchmark index.
		std::vector< std::pair<std::string,int> > tableKV;
				
		// List of benchmark results
		std::vector<BenchmarkInfo> benchmarks;

		void reset();
		void init(const std::string &filename);
		void execute(const std::map<std::string, MyTestFn> &functions,
					 const UpdateStatsMyTestFn &UpdateStatsFn,
					 const std::vector<std::string> &tags,
					 const std::vector<std::string> &rings,
					 const std::vector<std::string> &orderings,
					 const bool MakeIdealHomog = false);

        void ExportAsPlainText();
		void ExportAsLatex( const std::string &filename = "benchmarks.tex",
							const bool ShowComparisonCharts = true,
							const bool ShowPerfBenchmarksByAlgo = true,
							const bool ShowTimeComparisonTable = true,
							const bool ShowExamplesDesc = true);

		ideal Homogenized( const ideal& j );

		private: 
			void timestamp();
			std::string GetBenchmarkExampleLabel( int, int );		
			void PrintLongTable(const std::string caption,
								const std::string title,
								const bool IsComparison,
								const int AlgorithmIndex,
								std::stringstream &out);
	};

	class TimeTag
	{
		public:
			TimeTag(const std::string & name) : 
				name(name),
				start(boost::date_time::microsec_clock<boost::posix_time::ptime>::local_time())
			{}
	
			~TimeTag()
			{
                          boost::posix_time::ptime now(date_time::microsec_clock<boost::posix_time::ptime>::local_time());
                          boost::posix_time::time_duration d = now - start;
		
				std::cout << name;
				std::cout << " completed in " << d.total_milliseconds() / 1000.0;
				std::cout << " seconds" << std::endl;
			}
	
		private:
			std::string name;
			boost::posix_time::ptime start;
	};
	
}

#endif
#endif
