//   Copyright (c)  2017  John Abbott, Anna M. Bigatti
//   Author:  2017  John Abbott

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


#include "CoCoA/ExternalLibs.H"
#include "CoCoA/PREPROCESSOR_DEFNS.H"

#include "gmp.h"
#include "CoCoA/ExternalLibs-Frobby.H"
#include "CoCoA/ExternalLibs-GFan.H"
#include "CoCoA/ExternalLibs-GSL.H"
#include "CoCoA/ExternalLibs-MathSAT.H"
#include "CoCoA/ExternalLibs-Normaliz.H"

#include "CoCoA/ToString.H"

#include <iostream>
//#include <vector>
using std::vector;


namespace CoCoA
{

  ExternalLibInfo::ExternalLibInfo(const std::string& LibName,
                                   const std::string& LibFile,
                                   const std::string& LibVersion,
                                   const std::string& LibWebsite):
        myName(LibName),
        myFile(LibFile),
        myVersion(LibVersion),
        myWebsite(LibWebsite)
  {}

  
  std::ostream& operator<<(std::ostream& out, const ExternalLibInfo& info)
  {
    if (!out) return out;  // short-cut for bad ostreams
    out << "ExternalLibInfo(LibName=\"" << info.myName << "\",  "
        << "LibFile=\"" << info.myFile << "\",  "
        << "LibVersion=\"" << info.myVersion << "\",  "
        << "LibInfo=\"" << info.myWebsite << "\")";
    return out;
  }


  const std::vector<ExternalLibInfo>& ExternalLibs()
  {
    static std::vector<ExternalLibInfo> info;
    if (!info.empty()) return info;
#ifdef __GNU_MP_RELEASE
    info.push_back(ExternalLibInfo("GMP","UNKNOWN libgmp",ToString(__GNU_MP_RELEASE),"https://gmplib.org/"));
#else
    info.push_back(ExternalLibInfo("GMP","UNKNOWN libgmp","!!OLD VERSION!!","https://gmplib.org/"));
#endif
    #ifdef CoCoA_WITH_NORMALIZ
    info.push_back(ExternalLibInfo("Normaliz","UNKNOWN libnormaliz",ToString(NMZ_RELEASE),"https://www.normaliz.uni-osnabrueck.de/"));
    #endif
    #ifdef CoCoA_WITH_FROBBY
    info.push_back(ExternalLibInfo("Frobby","UNKNOWN libfrobby", FrbVersion() /* >= 0.9.3*/,"http://www.broune.com/frobby/"));
    #endif
    #ifdef CoCoA_WITH_GFAN
    info.push_back(ExternalLibInfo("Gfan","UNKNOWN libgfan","UNKNOWN" /* >= 0.6*/,"http://home.math.au.dk/jensen/software/gfan/gfan.html"));
    info.push_back(ExternalLibInfo("CDD","UNKNOWN libcddgmp","UNKNOWN" /* tested only with 094h*/,"https://www.inf.ethz.ch/personal/fukudak/cdd_home/index.html"));
    #endif
    #ifdef CoCoA_WITH_GSL
    //    info.push_back(ExternalLibInfo("GSL","UNKNOWN libgsl","UNKNOWN" /* >= ?*/,"https://www.gnu.org/software/gsl/"));
    #endif
    #ifdef CoCoA_WITH_MATHSAT
    char* const MSAT_VERSION = msat_get_version();
    info.push_back(ExternalLibInfo("MathSAT","UNKNOWN libmathsat",ToString(MSAT_VERSION) ,"http://mathsat.fbk.eu"));
    msat_free(MSAT_VERSION);
    #endif
    return info;
  }


} // end of namespace CoCoA
