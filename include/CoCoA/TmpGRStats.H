#ifndef CoCoA_GRStats_H
#define CoCoA_GRStats_H

//   Copyright (c)  2005  John Abbott,  Anna M. Bigatti
//   Original author: 2005 Massimo Caboara

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


#include "CoCoA/degree.H"

#include <iosfwd>
// using std::ostream;
#include <list>
// using std::list;

namespace CoCoA
{

  class DegStats
  {
  public:

    DegStats(const degree& Deg);
    DegStats(const degree& Deg,
             unsigned int PInserted,
             unsigned int GMLilled,
             unsigned int CopKilled,
             unsigned int BKilled,
             unsigned int Useful,
             unsigned int Useless,
             unsigned int Pairsno);
    ~DegStats(){};
    // fields
    degree myDeg;
    unsigned int myPInserted;
    unsigned int myGMKilled;
    unsigned int myCopKilled;
    unsigned int myBKilled;
    unsigned int myUseful;
    unsigned int myUseless;
    unsigned int myPairsNo;
  };


//------------------------Stats--------------------------------

  class Stats
  {
  public:
    Stats(unsigned int NumGens);// standard ctor
    Stats(unsigned int NumGens,// detailed ctor
          unsigned int PInserted,
          unsigned int GMKilled,
          unsigned int CopKilled,
          unsigned int Useful,
          unsigned int Useless,
          unsigned int BKilled,
          unsigned int BTouched,
          unsigned int GMTouched,
          unsigned int PolyDeleted,
          unsigned int PolyDHed,
          unsigned int DegDH,
 	  unsigned int NReductions,
          double ReductionTime,
          double TotalTime,
          long PrintReduction,
          long PrintDeg,
          long PrintGM,
          long PrintCop,
          long PrintBC,
          long PrintNumPair,
          long PrintFinal,
          long PrintFinalFull,
          long PrintFinalSimple,
          long PrintNewPairs,
          long PrintPolyDeleted,
          long PrintPolyDH,
          long PrintPolyLen,
          long PrintKill);
    ~Stats(){};
	
    void mySetLevels();
    void myUpgradeDegStats(const degree& new_deg, unsigned int pairs_no);
    void myStampa(std::ostream&)const;
    Stats& operator=(const Stats& rhs);
    

    // These are the statistics data. They are computed in any case
    unsigned int myPInserted;///< Number of Pairs Inserted
    unsigned int myGMKilled;///< Number of Pairs killed by the GM criteria
    unsigned int myCopKilled;///< Number of Pairs killed by the coprime criteria
    unsigned int myUseful;///< Number of Pairs that doesn't reduce to 0
    unsigned int myUseless;///< Number of Pairs that reduce to 0
    unsigned int myBKilled;///< Number of Pairs killed by the Back criterion
    unsigned int myBTouched;///< Number of pairs examined by the Back criterion
    unsigned int myGMTouched;///< Number of Pairs examined by the GM criteria
    unsigned int myPolyDeleted;///< Nuymber of polys deleted by new polys - 0 if the example is homog
    unsigned int myPolyDHed;///< Nuymber of polys dehomog'ed by dehomogalgorithm - pseudo homog only
    unsigned int myDegDH;///< Total of degrees dropped by the dehomog algorithm
    unsigned int myNumGens;///< Nuymber of the initial generators of the GBasis
    unsigned int myNReductions;///< Nuymber of reductions during computations
    double myReductionTime;// Used for statistics.
    double myTotalTime;// Used for statistics.

    std::list<DegStats> myDegByDeg;///< A list of deg-by deg stats computed during execution
    std::list<std::pair<unsigned int,unsigned int> > myPolyLens;///< List of <len Spoly,len reduced spoly>

    // These are the fine statistics settings. They control what is printed.
    // Default values: you can manipulate it by the
    // set_stat_level procedure. Fine tuning can be done by changing the values
    // one by one in the GReductor ctor
    long myReductionLevel;///< Print reductiob stats during Buchberger
    long myDegLevel;///< Print new deg when deg change during Buchberger
    long myCopLevel;///< Print coprime killed during Buchberger
    long myGMLevel;///< Print GM criteria stats during Buchberger
    long myBCLevel;///< Print Back criterion stats during Buchberger
    long myKillLevel;///< Print KIll criterion stats during Buchberger
    long myNumPairLevel;///< Print number of pairs after every reduction
    long myNewPairsLevel;///< Print data about new pairs creation during Buchberger
    long myPolyDeletedLevel;///< Prints the number of polys deleted by new polys -
    long myPolyDHLevel;///< Prints the number of polys dhmog'ed - new dehomog alg
    long myFinalLevel;///< Print stat data after Buchberger
    long myFinalFullLevel;///< Print full stat data degree by degree after Buchberger
    long myFinalSimpleLevel;///< Print minimal stat data after Buchberger
    long myPolyLenLevel;///< Print poly len stats after Buchberger


    //  private:
    //    int myLevel;

  };


}// end namespace cocoa

#endif
