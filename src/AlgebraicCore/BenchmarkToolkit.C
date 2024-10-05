#ifdef CoCoA_WITH_BOOST_temporarily_disabled // needs boost 1.41 or later

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

#include "CoCoA/BenchmarkToolkit.H"
#include "CoCoA/BigInt.H"
#include "CoCoA/CanonicalHom.H"
#include "CoCoA/DenseMatrix.H"
#include "CoCoA/PolyRing.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/TmpHilbert.H"
#include "CoCoA/VectorOps.H"
#include "CoCoA/matrix.H"
#include "CoCoA/symbol.H"
#include "CoCoA/time.H"
//#include "CoCoA/library.H"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp> // split function

#include <boost/thread.hpp>
#include <boost/date_time/posix_time/posix_time_duration.hpp>
#include <boost/bind.hpp>
#include <boost/ref.hpp>

#include <iostream>
#include <vector>
#include <iomanip>
#include <sstream>
#include <algorithm>

#define BENCHMARK_DEBUGGER 0

using namespace CoCoA;

namespace CoCoA
{

	/**
	 * Creates an unique label for the benchmark example.
	 */

	std::string BenchmarkToolkit::GetBenchmarkExampleLabel(int i, int j)
	{
		BenchmarkExample b = database[i];
		std::stringstream ss;
		ss << b.id << "-" << b.TestRings[j].ordering << "-" << b.TestRings[j].id;
		return ss.str();
	}

	/**
	 * Reset stored results.
	 */

	void BenchmarkToolkit::reset()
	{
		testset.clear();
		benchmarks.clear();
	}

	/**
	 * Print benchmark results to std::cout.
	 */

	void BenchmarkToolkit::ExportAsPlainText()
	{
		int i, j, size;

		if(testset.size() == 0 || header.size() == 0) return;

		/**
		 * Stats for each Groebner Basis Algorithm
		 */

		size = 20 * header.size() + 1;

		BOOST_FOREACH( BenchmarkInfo const &b, benchmarks )
		{
			std::cout << std::endl << "Algorithm: " << b.algorithm;
			std::cout << std::endl << std::setfill ('-') << std::setw (size) << " ";
			std::cout << std::endl << "|" << std::setfill (' ');

			BOOST_FOREACH( std::string const &e, header )
				std::cout << std::setw(17) << e << std::setw(3) << " | ";

			std::cout << std::endl << std::setfill ('-') << std::setw (size) << " ";

			BOOST_FOREACH( std::vector<std::string> const &r, b.stats )
			{
				std::cout << std::endl << "|" << std::setfill (' ');
				for( i = 0; i < (int)r.size(); i++ )
				{
					std::string const e = r[i];
					std::cout  << std::setw(17) << e << std::setw(3) << " | ";
				}
			}
			std::cout << std::endl << std::setfill ('-');
			std::cout << std::setw (size) << " " << std::endl;
		}
		std::cout << std::endl;

		/**
		 * Comparison stats for each attribute
		 */

		size = 46 + 20 * benchmarks.size();

		for( i = 0; i < (int)tableKV.size(); i++ )
		{
			std::pair<std::string, int> p = tableKV[i];

			std::cout << std::endl << "Table: " << p.first;

			std::cout << std::endl << std::setfill ('-') << std::setw (size) << " ";
			std::cout << std::endl << "|" << std::setfill (' ');
			std::cout << std::setw(42) << "Algorithm" << std::setw(3) << " | ";

			BOOST_FOREACH( BenchmarkInfo const &b, benchmarks )
				std::cout << std::setw(17) << b.algorithm << std::setw(3) << " | ";

			std::cout << std::endl << std::setfill ('-') << std::setw (size) << " ";

			for( j = 0; j < (int)testset.size(); j++ )
			{
				std::pair<int,int> ts = testset[j];
				std::cout << std::endl << "|" << std::setfill (' ');
				std::cout << std::setw(42) << GetBenchmarkExampleLabel(ts.first, ts.second);
				std::cout << std::setw(3) << " | ";

				BOOST_FOREACH( BenchmarkInfo const &b, benchmarks )
				{
					std::cout << std::setw(17) << b.stats[j][p.second] << std::setw(3) << " | ";
				}
			}

			std::cout << std::endl << std::setfill ('-');
			std::cout << std::setw (size) << " " << std::endl;
			std::cout << std::endl;
		}
	}

	void BenchmarkToolkit::PrintLongTable(const std::string caption,
											const std::string title,
											const bool IsComparison,
											const int DataIndex,
											std::stringstream &out)
	{
		int j, k, i;

		k = (IsComparison) ? (int)benchmarks.size() : (int) header.size();

		/**
		 * Table Header
		 */

		std::stringstream myHeader;
		myHeader << "\\toprule & \\multicolumn{" <<k<< "}{c}{";
		myHeader << title << "}\\\\";
		myHeader << std::endl << "\\cmidrule(l){2-" << (k+1) << "}";

		if(IsComparison)
		{
			myHeader << std::endl << "Algorithm ";
			for(j = 0; j < k; j++)
			{
				myHeader << "&" << benchmarks[j].algorithm;
			}
		}
		else
		{
			myHeader << header[0];
			for( j = 1; j < k; j++ )
			{
				myHeader << " & " << header[j];
			}
		}

		myHeader << "\\\\" << std::endl << "\\midrule " << std::endl;

		/**
		 * LaTeX Table
		 */

		out << std::endl << std::endl;
		out << "\\begin{longtable}{l";

		for( j = 0; j < (int) header.size(); j++ )
		{
			out << "c";
		}

		out << "}" << std::endl;
		out << "\\caption{"<< caption << "} \\\\" << std::endl;

		out << myHeader.str();
		out << "\\endfirsthead" << std::endl << std::endl;

		out << "\\multicolumn{" << (k+1) << "}{c}{\\tablename\\ \\thetable\\";
		out << "-- \\textit{Continued from previous page}} \\\\" << std::endl;
		out << myHeader.str();
		out << "\\endhead" << std::endl << std::endl;

		out << "\\hline \\multicolumn{" << (k+1);
		out << "}{r}{\\textit{Continued on next page}} \\\\" << std::endl;
		out << "\\endfoot" << std::endl << std::endl;

		out << "\\hline" << std::endl;
		out << "\\endlastfoot" << std::endl;

		if(IsComparison)
		{
			for( j = 0; j < (int)testset.size(); j++ )
			{
				std::pair<int,int> ts = testset[j];
				out << GetBenchmarkExampleLabel(ts.first, ts.second);
				BOOST_FOREACH( BenchmarkInfo const &b, benchmarks )
				{
					out << "&" << b.stats[j][DataIndex];
				}
				out << "\\\\" << std::endl;
				j++;
			}
		}
		else
		{
			BenchmarkInfo const &b = benchmarks[DataIndex];
			for(j=0; j < (int)b.stats.size(); j++)
			{
				for(i=0; i < (int)b.stats[j].size()-1; i++)
				{
					out << ((i > 0) ? " & " : "") << b.stats[j][i];
				}
				out << " & " << b.stats[j][i];
				out << "\\\\" << std::endl;
			}
		}

		out << "\\hline" << std::endl;
		out << "\\end{longtable}" << std::endl;
	}

	void BenchmarkToolkit::ExportAsLatex( const std::string &filename,
										  const bool ShowComparisonCharts,
										  const bool ShowPerfBenchmarksByAlgo,
										  const bool ShowTimeComparisonTable,
										  const bool ShowExamplesDesc )
	{
		if(testset.size() == 0) return;
		int j, i, k, f;

		k = (int)benchmarks.size();

		std::stringstream LatexCode;

		LatexCode << "\\documentclass{article}" << std::endl;
		LatexCode << "\\usepackage{booktabs}" << std::endl;
		LatexCode << "\\usepackage{graphicx}" << std::endl;
		LatexCode << "\\usepackage{longtable}" << std::endl;
		LatexCode << "\\usepackage{pgfplots}" << std::endl;
		LatexCode << "\\usepackage{multirow}" << std::endl << std::endl;

		LatexCode << "\\usepackage{anysize}" << std::endl;
		LatexCode << "\\marginsize{2cm}{2cm}{2cm}{2cm}" << std::endl;

		LatexCode << "\\begin{document}" << std::endl;
		LatexCode << "\\tableofcontents " << std::endl;

		// Charts - Time comparison
		if( ShowComparisonCharts )
		{
			LatexCode << std::endl << "\\newpage" << std::endl;
			LatexCode << "\\section{Benchmark Comparison Charts}" << std::endl;

			LatexCode << "\\pgfplotsset{compat=newest}" << std::endl;

			std::vector<std::vector<std::string> > rows = benchmarks[ 0 ].stats;

			int NumRows = (int)rows.size();
			double ChartHeight = 11;

			int NumBars = 55 / k;

			for( f = 0; f < (int)tableKV.size(); f++ )
			{
				std::pair<std::string, int> kv = tableKV[f];

				for(i = 0; i < NumRows; i += NumBars)
				{
					int p = std::min(i + NumBars, NumRows);
					double BarWidth = std::min(ChartHeight / (p * k) / k, 1.1);

					LatexCode << "\\newpage" << std::endl;
					LatexCode << "\\begin{figure}[h]" << std::endl;
					LatexCode << "\\caption{Benchmark Comparison Chart: "<< kv.first;
					LatexCode << "}\\vspace{0.5cm}" << std::endl << "\\begin{tikzpicture}" << std::endl;
					LatexCode << "\\begin{axis}[";
					LatexCode << "legend style={legend columns=" << k;
					LatexCode <<", at={(xticklabel cs:0.5,40)}, anchor=south,draw=none },";
					LatexCode << "height=" << ChartHeight << "cm, symbolic y coords={";

					std::pair<int,int> ts = testset[i];
					LatexCode << GetBenchmarkExampleLabel(ts.first, ts.second);

					for ( j = i+1; j < p; j++ )
					{
						ts = testset[j];
						LatexCode << ", " <<GetBenchmarkExampleLabel(ts.first, ts.second);
					}

					LatexCode << "},";
					LatexCode << "y tick label style={text width=5cm,font=\\footnotesize}, ytick=data, ";
					LatexCode << "xlabel=" << kv.first << ", xmin=0, xbar, bar width=" << BarWidth;
					LatexCode << "cm]" << std::endl;

					BOOST_FOREACH( BenchmarkInfo const &b, benchmarks )
					{
						LatexCode << "\\addplot coordinates {" << std::endl;
						LatexCode << "\t";

						for ( j = i; j < p; j++ )
						{
							std::vector<std::string> const &r = b.stats[j];
							ts = testset[j];
							LatexCode << "(" << r[kv.second] << "," << GetBenchmarkExampleLabel(ts.first, ts.second) << ") ";
						}

						LatexCode << std::endl << "};" << std::endl;
						LatexCode << "\\addlegendentry{"<< b.algorithm << "}" << std::endl;
					}

					LatexCode << "\\end{axis}" << std::endl;
					LatexCode << "\\end{tikzpicture}" << std::endl;
					LatexCode << "\\end{figure}" << std::endl;
				}
			}
		}

		/**
		 * Table - Benchmark by algorithm
		 */

		if(ShowPerfBenchmarksByAlgo)
		{
			LatexCode << std::endl;
			LatexCode << "\\newpage" << std::endl;
			LatexCode << "\\section{Benchmark Algorithm Statistics}" << std::endl;

			int r = 0;
			BOOST_FOREACH( BenchmarkInfo const &b, benchmarks )
			{
				std::string title(b.algorithm);
				std::string caption(b.algorithm);
				caption.append(" - Performance Benchmarks");
				PrintLongTable(caption, title, 0, r++, LatexCode);
			}
		}

		/**
		 * Table - Benchmark comparison
		 */

		if( ShowTimeComparisonTable && k > 1 )
		{
			LatexCode << std::endl;
			LatexCode << "\\newpage" << std::endl;
			LatexCode << "\\section{Benchmark Comparison Tables}" << std::endl;

			std::string title("Implementation");
			for( f = 0; f < (int)tableKV.size(); f++ )
			{
				std::pair<std::string, int> kv = tableKV[f];
				std::string caption(kv.first);
				caption.append(" - Comparison");
				PrintLongTable(caption, title, 1, kv.second, LatexCode);
			}
		}

		/**
		 * Table - Description of benchmark examples
		 */

		if(ShowExamplesDesc)
		{
			LatexCode << std::endl << "\\newpage" << std::endl;
			LatexCode << "\\appendix" << std::endl;
			LatexCode << "\\section{Benchmark Examples} \\label{App:BenchmarkToolkit}" << std::endl;

			LatexCode << "\\begin{longtable}{|l|l|p{10cm}|} " << std::endl;
			LatexCode << "%\\caption{Benchmark Examples}\\label{App:BenchmarkToolkit} \\\\" << std::endl;

			LatexCode << "\\hline" << std::endl;
			LatexCode << "\\multicolumn{3}{|c|}{Example} \\\\ \\hline";
			LatexCode << std::endl << "\\endfirsthead" << std::endl << std::endl;

			LatexCode << "\\hline" << std::endl;
			LatexCode << "\\multicolumn{3}{|c|}%" << std::endl;
			LatexCode << "{\\tablename\\ \\thetable\\ -- \\textit{Continued from previous page}} \\\\" << std::endl;
			LatexCode << "\\hline" << std::endl;
			LatexCode << "\\multicolumn{3}{|c|}{Example} \\\\" << std::endl;
			LatexCode << "\\hline" << std::endl;
			LatexCode << "\\endhead" << std::endl << std::endl;

			LatexCode << "\\hline \\multicolumn{3}{r}{\\textit{Continued on next page}} \\\\" << std::endl;
			LatexCode << "\\endfoot" << std::endl << std::endl;

			LatexCode << "\\hline" << std::endl;
			LatexCode << "\\endlastfoot" << std::endl << std::endl;

			LatexCode << "\\hline \\hline" << std::endl;

			BOOST_FOREACH( BenchmarkExample const &t, database )
			{
				int i = 1;
				i += ((!t.references.empty()) ? 1 : 0);
				i += ((!t.notes.empty()) ? 1 : 0);

				// No information about this example
				if(i == 1 && t.description.empty()) continue;

				if(i == 1)
				{
					LatexCode << "\\vspace*{0.2cm}" << t.name << "\\vspace*{0.2cm} &";
				}
				else
				{
					LatexCode << "\\multirow{" << i << "}{*}{\\vspace*{0.2cm}";
					LatexCode << t.name << "\\vspace*{0.2cm}} &";
				}

				LatexCode << "Description &" ;
				if(!t.description.empty())
				{
					LatexCode << "\\parbox{10cm}{\\vspace*{0.2cm}" << t.description;
					LatexCode << "\\vspace*{0.2cm}}\\\\" << std::endl;
				}
				else
					LatexCode << "N.A." << "\\\\" << std::endl;

				if(!t.references.empty())
				{
					LatexCode << "\\cline{2-3} & References &" ;
					LatexCode << "\\parbox{10cm}{\\vspace*{0.2cm}" << t.references;
					LatexCode << "\\vspace*{0.2cm}}\\\\" << std::endl;
				}

				if(!t.notes.empty())
				{
					LatexCode << "\\cline{2-3} & Notes &" ;
					LatexCode << "\\parbox{10cm}{\\vspace*{0.2cm}" << t.notes;
					LatexCode << "\\vspace*{0.2cm}}\\\\" << std::endl;
				}

				LatexCode << "\\hline \\hline" << std::endl << std::endl;
			}
			LatexCode << "\\hline" << std::endl;
			LatexCode << "\\end{longtable}" << std::endl;
		}

		LatexCode << "\\end{document}" << std::endl;

		std::string buff(LatexCode.str());
		std::replace( buff.begin(), buff.end(), '_', '-');

		std::string FilePath = "bin/"; FilePath.append( filename );
		std::ofstream LatexCodeContent ( FilePath.c_str() );
		LatexCodeContent << buff.c_str();
	}

	ring GetRingByName(const std::string &RingName, const int RingChar)
	{
		if(RingName.compare("NewZZmod") == 0)
		{
			return NewZZmod(RingChar);
		}
		else if(RingName.compare("RingZZ") == 0)
		{
			return RingZZ();
		}
		else if(RingName.compare("RingQQ") == 0)
		{
			return RingQQ();
		}
		return RingQQ();
	}

	RingElem ReadPoly(std::istringstream &in, const SparsePolyRing &P)
	{
		ring R = CoeffRing(P);
		long NumParams = 0;
		const long NumInds = NumIndets(P);

		if (IsFractionField(R) && IsPolyRing(BaseRing(R)))
		{
			NumParams = NumIndets(BaseRing(R));
		}

		std::vector<long> v(NumInds);
		std::vector<long> ParV(NumParams);

		long NumSummands;
		in >> NumSummands;
		if (!in) CoCoA_THROW_ERROR(ERR::InputFail, "ReadPoly -- NumSummands");
		BigInt IntCoeff; // to avoid ctor-dtor inside loop
		RingElem c(R), ans(P), mon(P);
		for (long NS=0; NS<NumSummands; ++NS)
		{
			in >> IntCoeff;
			for (long i=0 ; i<NumParams ; ++i)  in >> ParV[i];
			for (long i=0 ; i<NumInds ; ++i)  in >> v[i];
                        if (!in) CoCoA_THROW_ERROR(ERR::InputFail, "ReadPoly -- summand");
			if (NumParams==0)
				c = IntCoeff;
			else
			{
                                const SparsePolyRing BR = BaseRing(R);
				c = CanonicalHom(BR,R)(monomial(BR, IntCoeff, ParV));
			}
			mon = monomial(P, c, v);
			P->myAddClear(raw(ans), raw(mon)); // ANNA: use geobucket
		}
		return ans;
	}

	/**
	 * Loads a database of examples from a XML file.
	 */

	void BenchmarkToolkit::init(const std::string &filename){

		/** Create an empty boost property tree object */
		using boost::property_tree::iptree;
		iptree pt;

		/**
		 * Load the XML file into the boost property tree.
		 * If an error happens then an exception is thrown.
		 */
		read_xml(filename, pt);

		/** Ordering functions */
		std::map<std::string, PPOrdering (*)(long)> OrderingMap;
		OrderingMap["Lex"] = &NewLexOrdering;
   		OrderingMap["DegLex"] = &NewStdDegLexOrdering;
   		OrderingMap["DegRevLex"] = &NewStdDegRevLexOrdering;

		BOOST_FOREACH( iptree::value_type const &root, pt.get_child("Benchmarks") )
		{
			if( root.first == "Environments" )
			{
				BOOST_FOREACH( iptree::value_type const &envs, root.second )
				{
					if( envs.first == "Environment" )
					{
						const std::string TestSet = envs.second.get<std::string>("Name");

						#if(BENCHMARK_DEBUGGER)
							std::cout << TestSet << std::endl;
						#endif

						const std::string NrInd = envs.second.get<std::string>("NrIndeterminates");
						int NumIndets = atoi(NrInd.c_str());

						const iptree& rts = envs.second.get_child("Rings");
						const iptree& ords = envs.second.get_child("Orderings");

						BOOST_FOREACH( iptree::value_type const &rt, rts )
						{
							const std::string RingName = rt.second.get<std::string>("Name");
							const int RingChar = rt.second.get<int>("Characteristic", 0);
							std::stringstream ss;
							ss << RingName << "-" << RingChar;

							#if(BENCHMARK_DEBUGGER)
								std::cout << "\t" << ss.str() << std::endl;
							#endif

							ring rk = GetRingByName(RingName, RingChar);

							BOOST_FOREACH( iptree::value_type const &ord, ords )
							{
								const std::string RingOrd = ord.second.get<std::string>("OrderingType");

								#if(BENCHMARK_DEBUGGER)
									std::cout << "\t\t" << RingOrd << std::endl;
								#endif

                                                                //std::map<std::string,PPOrdering (*)(long)>::iterator it = OrderingMap.find(RingOrd);
                                                                const auto it = OrderingMap.find(RingOrd);
								PPOrdering o = (it != OrderingMap.end()) ? it->second(NumIndets) : NewStdDegRevLexOrdering(NumIndets);

								/** User-defined ordering and grading */
								if(RingOrd == "Matrix")
								{
									int mi, mj, i, j, val;

									const int GradingDim = ord.second.get("GradingDim", 0);
									const std::string s = ord.second.get<std::string>("IntegerMatrix");

									#if(BENCHMARK_DEBUGGER)
										std::cout << s << std::endl;
									#endif

									std::istringstream mat(s);

									mat >> mi;
									mat >> mj;
									if(mi != NumIndets)
										std::cerr << "Invalid Ordering matrix." << std::endl; // check this

									matrix M(NewDenseMat(RingZZ(), mi, mj));
									for(i=0; i<mi; i++)
									{
										for(j=0; j<mj; j++)
										{
											mat >> val;
											SetEntry(M, i, j, val);
										}
									}
									o = NewMatrixOrdering(NumIndets, GradingDim, M);
								}

								boost::optional< const iptree& > child = envs.second.get_child_optional( "Indeterminates" );
								if( child )
								{
									const std::string RingIndets = envs.second.get<std::string>("Indeterminates");
									std::vector<std::string> vars;
									boost::split(vars, RingIndets, boost::is_any_of(","));
									if(NumIndets == (int)vars.size())
									{
										ring r = NewPolyRing(rk, NewPPMonoid(symbols(vars),o));
										PolyRingInfo myRingInfo;
										myRingInfo.id = ss.str();
										myRingInfo.characteristic = RingChar;
										myRingInfo.name = RingName;
										myRingInfo.ordering = RingOrd;
										myRingInfo.ring.push_back(r);
										PrMap[TestSet].push_back(myRingInfo);
									}
									else
									{
										std::cerr << "Error 12" << std::endl;
									}
								}
								else
								{
									ring r = NewPolyRing(rk, NewPPMonoid(SymbolRange("x", 0, NumIndets-1),o));
									PolyRingInfo myRingInfo;
									myRingInfo.id = ss.str();
									myRingInfo.characteristic = RingChar;
									myRingInfo.name = RingName;
									myRingInfo.ordering = RingOrd;
									myRingInfo.ring.push_back(r);
									PrMap[TestSet].push_back(myRingInfo);
								}
							}
						}
					}
				}
			}

			else if( root.first == "Examples" )
			{
				BOOST_FOREACH( iptree::value_type const &node, root.second )
				{
					if( node.first == "Example" )
					{
						BenchmarkExample myBenchmarkExample;
						myBenchmarkExample.id = node.second.get<std::string>("ID");
						myBenchmarkExample.name = node.second.get<std::string>("Name");

						myBenchmarkExample.description = node.second.get<std::string>("Description", "");
						myBenchmarkExample.references = node.second.get<std::string>("References", "");
						myBenchmarkExample.notes = node.second.get<std::string>("Notes", "");

						const std::string tags = node.second.get<std::string>("TagsList");
						boost::split(myBenchmarkExample.TagsList, tags, boost::is_any_of(","));
						std::for_each(myBenchmarkExample.TagsList.begin(), myBenchmarkExample.TagsList.end(),
									boost::bind(&boost::trim<std::string>, _1, std::locale() ));
						std::sort(myBenchmarkExample.TagsList.begin(), myBenchmarkExample.TagsList.end());

						const std::string EnvironmentID = node.second.get<std::string>("TestEnvironment");
						const std::string polynomials = node.second.get<std::string>("Polynomials");
						int NumPolys;

						#if(BENCHMARK_DEBUGGER)
							std::cout << std::endl << "Example:\t" << myBenchmarkExample.name << std::endl;
							if(!myBenchmarkExample.description.empty()) std::cout << "Description:\t" << myBenchmarkExample.description << std::endl;
							if(!myBenchmarkExample.references.empty()) std::cout << "References:\t" << myBenchmarkExample.references << std::endl;
							if(!myBenchmarkExample.notes.empty()) std::cout << "Notes:\t\t" << myBenchmarkExample.notes << std::endl;
							std::cout << "Environment ID:\t" << EnvironmentID << std::endl;
						#endif

						if(PrMap.find( EnvironmentID ) == PrMap.end())
						{
							std::cerr << "Test environment is missing: " << EnvironmentID << std::endl;
						}
						else
						{
							std::vector<PolyRingInfo> RingInfo = PrMap[EnvironmentID];
							BOOST_FOREACH( PolyRingInfo myRingInfo, RingInfo )
							{
								SparsePolyRing myRing = myRingInfo.ring[0];
								std::vector<RingElem> j;

								std::istringstream in(polynomials);
								in >> NumPolys;
                                                                if (!in) CoCoA_THROW_ERROR(ERR::InputFail, "ReadPolyList -- NumPolys");
								j.reserve(NumPolys);

								for (int i=0 ; i<NumPolys ; ++i)
								{
									j.push_back(ReadPoly(in, myRing));
								}

								ideal I(myRing,j);
								myBenchmarkExample.TestIdeals.push_back(I);
							}

							myBenchmarkExample.TestRings = RingInfo;
							database.push_back(myBenchmarkExample);
						}
					}
				}
			}
    	}
	}

	// Do some pre-processing to the ideal before computing the GB
	// In this case we are going to homog it.
	ideal BenchmarkToolkit::Homogenized( const ideal& j )
	{
		std::vector<RingElem> g = gens(j);
		SparsePolyRing P = owner(g[0]);

		// Update PPOrdering
		PPOrdering PPO = ordering(PPM(P));
		int NewNumIndets = NumIndets(PPO) + 1;

		if (IsLex(PPO)) PPO = NewStdDegLexOrdering(NewNumIndets);
		else if (IsStdDegLex(PPO)) PPO = NewStdDegLexOrdering(NewNumIndets);
		else if (IsStdDegRevLex(PPO)) PPO = NewStdDegRevLexOrdering(NewNumIndets);
		else { /** Implement */}

		// Create a new ring
		SparsePolyRing K = NewPolyRing(CoeffRing(P), NewPPMonoid(NewSymbols(NewNumIndets),PPO));

		// Create a RingHom
		std::vector<RingElem> images = indets(K); images.resize(NewNumIndets-1);
		RingHom phi = PolyAlgebraHom(P, K, images);

		// Perform the homogenization
		RingElem h = indet(K, NewNumIndets-1);
		return homog(ideal(phi(g)), h);
	}

	void BenchmarkToolkit::timestamp()
	{
		time_t ltime; /* calendar time */
		ltime = time(NULL); /* get current cal time */
		std::cout << std::endl;
		char *t = asctime(localtime(&ltime));
		t[strlen(t)-1] = '\0';
		std::cout << "[" << t << "] ";
	}

	void BenchmarkToolkit::execute(const std::map<std::string, MyTestFn> &functions,
									const UpdateStatsMyTestFn &UpdateStatsFn,
									const std::vector<std::string> &tags,
									const std::vector<std::string> &rings,
									const std::vector<std::string> &orderings,
									const bool MakeIdealHomog)
	{

		int i, j;

		// Clear previous results
		reset();

		if( (int)functions.size() == 0 ) return;

		// Apply filters
		for( i = 0; i < (int)database.size(); i++ )
		{
			BenchmarkExample t = database[i];

			if(tags.size() > 0)
			{
				bool AcceptExample = false;
				BOOST_FOREACH( std::string tag, tags )
				{
					if(std::binary_search(t.TagsList.begin(), t.TagsList.end(), tag))
					{
						AcceptExample = true;
						break;
					}
				}
				if(!AcceptExample) continue;
			}

			for(j = 0; j < (int)t.TestRings.size(); j++)
			{
				PolyRingInfo ri = t.TestRings[j];

				bool AcceptOrd = false;
				if (orderings.size() == 0 ||
					 std::binary_search(orderings.begin(), orderings.end(), ri.ordering))
				{
					AcceptOrd = true;
				}

				bool AcceptRing = false;
				if (rings.size() == 0 ||
					 std::binary_search(rings.begin(), rings.end(), ri.id))
				{
					AcceptRing = true;
				}

				if( AcceptOrd && AcceptRing )
				{
					std::pair <int, int> p (i,j);
					testset.push_back(p);
				}
			}
		}

		std::pair<std::string,MyTestFn> p;
		BOOST_FOREACH(p, functions)
		{
			BenchmarkInfo info;
			info.algorithm = p.first;
			for( i = 0; i < (int)testset.size(); i++ )
			{
				std::pair <int,int> pair = testset[i];
				BenchmarkExample t = database[pair.first];
				j = pair.second;

				timestamp();
				std::cout << p.first << "(" << t.name << ", " << t.TestRings[j].id << ");" << std::flush;

				std::vector<RingElem> gb;
				ideal myIdeal( (MakeIdealHomog) ? Homogenized(t.TestIdeals[j]) : t.TestIdeals[j] );	HilbertNumQuot(myIdeal);
        double timer = CpuTime();
				p.second(myIdeal, gb);
				std::vector<std::string> stat;
        timer = CpuTime()-timer;
				stat = UpdateStatsFn(t, gb, timer, j, p.first);
				info.stats.push_back(stat);
				gb.clear();

			}
			benchmarks.push_back(info);
		}
	}
}
#endif
