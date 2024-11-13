//   Copyright (c)  2005,2007,2012  John Abbott and Anna M. Bigatti

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


#include "CoCoA/BuildInfo.H"
#include "CoCoA/config.H"
#include "CoCoA/long32or64.H"
#include "CoCoA/PREPROCESSOR_DEFNS.H"

#include <limits>
#include <iostream>
using std::ostream;
using std::endl;
#include <string>
using std::string;

namespace CoCoA
{

  namespace BuildInfo
  {
    // Build info is actually passed in via three preprocessor variables:
    // COCOA_VERSION, COCOA_CXX, COCOA_CXXFLAGS -- all three are quoted strings.

    const std::string& version()
    {
      static const string info(COCOA_VERSION);
      return info;
    }

    const std::string& compiler()
    {
      static const string info(COCOA_CXX);
      return info;
    }

    const std::string& CompilationFlags()
    {
      static const string info(COCOA_CXXFLAGS);
      return info;
    }

    const std::string& CompilationPlatform()
    {
      static const string info(COCOA_PLATFORM);
      return info;
    }

    std::string CoCoA_32BIT_LONG_Flag()
    {
#ifdef CoCoA_32BIT_LONG
      return "defined";
#else
      return "undefined (==> 64-bit longs)";
#endif
    }


    namespace // anonymous
    {
      std::string GetDefinesFromHeaderFile()
      {
        // See which options were chosen in PREPROCESSOR_DEFNS.H
        string ans;
#ifdef CoCoA_DEBUG
        ans += "CoCoA_DEBUG  ";
#endif
#ifdef CoCoA_MEMPOOL_DEBUG
        ans += "CoCoA_MEMPOOL_DEBUG  ";
#endif
#ifdef CoCoA_MEMPOOL_DISABLE
        ans += "CoCoA_MEMPOOL_DISABLE  ";
#endif
        return ans;
      }
    }

    const std::string& CompilationPreprocessorDefines()
    {
      static const string info = GetDefinesFromHeaderFile();
      return info;
    }


    void PrintAll(std::ostream& out)
    {
      if (!out) return;  // short-cut for bad ostreams
      // Below we use string literal juxtaposition, so that the build
      // info can readily be extracted from libcocoa.a by the command
      // strings libcocoa.a | grep -F "CoCoA::BuildInfo"
      out << endl
          << "CoCoA::BuildInfo Summary of build information for CoCoALib:" << endl
          << "CoCoA::BuildInfo CoCoALib Version: " << version() << endl
          << "CoCoA::BuildInfo Compiler: " << compiler() << endl
          << "CoCoA::BuildInfo Compilation Flags: " << CompilationFlags() << endl
          << "CoCoA::BuildInfo Compilation Platform: " << CompilationPlatform() << endl
          << "CoCoA::BuildInfo CoCoA_32BIT_LONG: " << CoCoA_32BIT_LONG_Flag() << endl
          << "CoCoA::BuildInfo Compilation Preprocessor Defines: " << CompilationPreprocessorDefines() << endl
          << "CoCoA::BuildInfo NumBits int = " << std::numeric_limits<unsigned int>::digits << endl
          << "CoCoA::BuildInfo NumBits long = " << std::numeric_limits<unsigned long>::digits << endl
          << "CoCoA::BuildInfo NumBits SmallExponent_t = " << std::numeric_limits<SmallExponent_t>::digits << endl
          << endl;
    }

  } // end of namespace CoCoA::BuildInfo

} // end of namespace CoCoA
