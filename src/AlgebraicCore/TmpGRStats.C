//   Copyright (c)  2005-2017  John Abbott and Anna M. Bigatti
//   Author: 2005  Massimo Caboara, 2017  Anna M. Bigatti

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

#include "CoCoA/TmpGRStats.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/VectorOps.H"  // for templates which print lists/vectors (in myStampaPairs)
#include "CoCoA/verbose.H"

#include <iostream>
using std::ostream;
using std::endl;
#include <utility>
using std::pair;
#include <list>
using std::list;

namespace CoCoA
{


/*----------------- class DegStats functions ----------------*/

  DegStats::DegStats(const degree& Deg):
    myDeg(Deg)
  {
    myPInserted=0;
    myGMKilled=0;
    myCopKilled=0;
    myBKilled=0;
    myUseful=0;
    myUseless=0;
    myPairsNo=0;
  }//ctor
  
  DegStats::DegStats(const degree& Deg,
           unsigned int PInserted,
	   unsigned int GMKilled,
           unsigned int CopKilled,
	   unsigned int BKilled,
	   unsigned int Useful,
	   unsigned int Useless,
           unsigned int PairsNo):
    myDeg(Deg)
  {
    myPInserted=PInserted;
    myGMKilled=GMKilled;
    myCopKilled=CopKilled;
    myBKilled=BKilled;
    myUseful=Useful;
    myUseless=Useless;
    myPairsNo=PairsNo;
  }//ctor

/*----------------- class Stats functions ----------------*/


  Stats::Stats(unsigned int NumGens)
  {
    mySetLevels();// Sets the flags. Has to be there.
    myNumGens=NumGens;
    myPInserted=0;
    myGMKilled=0;
    myCopKilled=0;
    myBKilled=0;
    myBTouched=0;
    myGMTouched=0;
    myUseful=0;
    myUseless=0;
    myPolyDeleted=0;
    myPolyDHed=0;
    myDegDH=0;
    myNReductions=0;
    myReductionTime=0.0;
    myTotalTime=0.0;
   }

// WARN: copy all the data, level may be out of sync with the flags
 Stats::Stats(unsigned int NumGens,
              //        int level,
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
              long PrintKill)
{
      myNumGens=NumGens;
      myPInserted=PInserted;
      myGMKilled=GMKilled;
      myCopKilled=CopKilled;
      myBKilled=BKilled;
      myBTouched=BTouched;
      myGMTouched=GMTouched;
      myUseful=Useful;
      myUseless=Useless;
      myPolyDeleted=PolyDeleted;
      myPolyDHed=PolyDHed;
      myDegDH=DegDH;
      myNReductions=NReductions,
      myReductionTime=ReductionTime;
      myTotalTime=TotalTime;
      myReductionLevel=PrintReduction;
      myDegLevel=PrintDeg;
      myGMLevel=PrintGM;
      myCopLevel=PrintCop;
      myBCLevel=PrintBC;
      myNumPairLevel=PrintNumPair;
      myFinalLevel=PrintFinal;
      myFinalFullLevel=PrintFinalFull;
      myFinalSimpleLevel=PrintFinalSimple;
      myNewPairsLevel=PrintNewPairs;
      myPolyDeletedLevel=PrintPolyDeleted;
      myPolyDHLevel=PrintPolyDH;
      myPolyLenLevel=PrintPolyLen;
      myKillLevel=PrintKill;
}

  void Stats::mySetLevels()
  {
    //        case 6:// poly reduction history
          myPolyLenLevel =140;
          myFinalFullLevel =140;
          //        case 5://
          myGMLevel =130;
          myBCLevel =130;
          myPolyDeletedLevel =130;
          myPolyDHLevel =130;
          myNewPairsLevel =130;
          myCopLevel =130;
          //        case 4:// pair by pair stats
          myKillLevel =120;
          myReductionLevel =120;
          //        case 3:// deg by deg and number of pairs
          myNumPairLevel =120;
          myDegLevel =100;
          //        case 2:// disable at the moment
          //myFinalFullLevel=true;
          //        case 1:// some final stats
          myFinalLevel =120;
          //          break;
          //	case 0:
          myFinalSimpleLevel =120;
//           break;
//         default:
//         case -1:
//           break;
//    };
  }

  void Stats::myUpgradeDegStats(const degree& new_deg, unsigned int Pairs_no)
  {
    DegStats New(new_deg);
    //  New.Deg=new_deg;
// for (const auto& entry: myDegByDeg)
    for (list<DegStats>::const_iterator it=myDegByDeg.begin();
         it!=myDegByDeg.end();
         ++it)
    {
      New.myPInserted+=it->myPInserted;
      New.myGMKilled+=it->myGMKilled;
      New.myCopKilled+=it->myCopKilled;
      New.myUseful+=it->myUseful;
      New.myUseless+=it->myUseless;
      New.myBKilled+=it->myBKilled;
    }
    New.myPInserted=myPInserted-New.myPInserted;
    New.myGMKilled=myGMKilled-New.myGMKilled;
    New.myCopKilled=myCopKilled-New.myCopKilled;
    New.myBKilled=myBKilled-New.myBKilled;
    New.myUseful=myUseful-New.myUseful;
    New.myUseless=myUseless-New.myUseless;
    New.myPairsNo=Pairs_no;

    myDegByDeg.push_back(New);
  }


  void Stats::myStampa(ostream& out)const
  {
    if (!out) return;  // short-cut for bad ostreams

    out<<"-- GBasis stats -----------------\n";
    if (VerbosityLevel() >= myFinalFullLevel)
    {

      if (!myDegByDeg.empty()) {
        out<<"The degree by degree Stats\n";
        for (list<DegStats>::const_iterator it=myDegByDeg.begin();
             it!=myDegByDeg.end();++it) {
          out<<"Degree "<<it->myDeg<<endl;
          out<<"** Pairs Reduced         "<<it->myUseful+it->myUseless
              <<" of which "<<it->myUseful<<" useful and "
              <<it->myUseless<<" useless\n";
          out<<"** Pairs Inserted        "<<it->myPInserted;
          if (it->myPInserted!=0)	
	  out<<" of which "<<it->myPInserted-it->myGMKilled-it->myCopKilled<<" survived";
          out<<endl;
          if (it->myPInserted!=0) {
            out<<"****  of which GMKilled    "<<it->myGMKilled<<endl;
            out<<"****  of which CopKilled   "<<it->myCopKilled<<endl;
            out<<"****  of which BKilled     "<<it->myBKilled<<endl;	
          };
          out<<"Pairs at end of deg      "<<it->myPairsNo<<endl;
          out<<"\n";
        };//for
      };//if
      if (VerbosityLevel() >= myPolyLenLevel)
      {
        out<< "Poly Lens\n";
        out<<"[[";
        list< pair<unsigned int, unsigned int> >::const_iterator it1;
        for (it1=myPolyLens.begin();it1!=myPolyLens.end();++it1)
          out<<",["<<it1->first<<","<<it1->second<<"]";
        out<<"]"<<endl;
      }
    };// if (myFinalFullLevel)

    if (VerbosityLevel() >= myFinalLevel)
    {
      out << " Poly in Basis  "<<myUseful-myPolyDeleted<<endl;
      out << " Pairs Reduced        "<<myUseful+myUseless
          << " = ("<<myUseful+myUseless-myNumGens<<"+"<<myNumGens<<" Gens)"
          << " of which "<<myUseful<<" useful and "
          << myUseless<<" useless\n";
      if (myBKilled==0) out<<" Minimal Pairs"<<myUseful+myUseless+myCopKilled-myNumGens
                           <<" = (Reduced "<<myUseful+myUseless<<" + Coprime "
                           <<myCopKilled<<" - Gens "<<myNumGens<<")"<<endl;
      if (myPolyDeleted!=0)
        out<<" Poly Deleted         "<<myPolyDeleted<<endl;
      if (myPolyDHed!=0)
      {
        out<<" Poly Dehomog'ed      "<<myPolyDHed<<endl;// number of polynomials dehomogenized
        out<<" Degs dropped           "<<myDegDH<<endl;// sum of the degrees dropped
      }
      out<<" Pairs Inserted        "<<myPInserted
          <<" + "<<myNumGens<<" generators"
          <<" of which "<<myUseful+myUseless<<" survived"<<endl;
      out<<"   of which GMKilled     "<<myGMKilled<<endl;
      out<<"   of which CopKilled    "<<myCopKilled<<endl;
      out<<"   of which BKilled      "<<myBKilled<<endl;
      out<<" GM considered pairs   "<<myGMTouched<<endl;

      out<<" Back considered pairs "<<myBTouched<<endl;
     }//if (myFinalLevel)
     if (VerbosityLevel() >= myFinalSimpleLevel)
      {
        out << "[log]  reductions=" << myNReductions;
        out << "\tTotalTime="<<-myTotalTime<< endl;
      }//if (myFinalSimpleLevel)
     out<<"-- GBasis stats end -----------------" << endl;
  }//myStampa


  Stats& Stats::operator=(const Stats& rhs)
  {
    if (this == &rhs) return *this;
    myNumGens = rhs.myNumGens;
    myPInserted = rhs.myPInserted;
    myGMKilled = rhs.myGMKilled;
    myCopKilled = rhs.myCopKilled;
    myBKilled = rhs.myBKilled;
    myBTouched = rhs.myBTouched;
    myGMTouched = rhs.myGMTouched;
    myUseful = rhs.myUseful;
    myUseless = rhs.myUseless;
    myPolyDeleted = rhs.myPolyDeleted;
    myPolyDHed = rhs.myPolyDHed;
    myDegDH = rhs.myDegDH;
    myNReductions = rhs.myNReductions;
    myReductionTime = rhs.myReductionTime;
    myTotalTime = rhs.myTotalTime;
    myReductionLevel = rhs.myReductionLevel;
    myDegLevel = rhs.myDegLevel;
    myGMLevel = rhs.myGMLevel;
    myCopLevel = rhs.myCopLevel;
    myBCLevel = rhs.myBCLevel;
    myNumPairLevel = rhs.myNumPairLevel;
    myFinalLevel = rhs.myFinalLevel;
    myFinalFullLevel = rhs.myFinalFullLevel;
    myFinalSimpleLevel = rhs.myFinalSimpleLevel;
    myNewPairsLevel = rhs.myNewPairsLevel;
    myPolyDeletedLevel = rhs.myPolyDeletedLevel;
    myPolyDHLevel = rhs.myPolyDHLevel;
    myPolyLenLevel = rhs.myPolyLenLevel;
    myKillLevel = rhs.myKillLevel;
    
    return *this;
  }
  


}// end namespace cocoa
		
