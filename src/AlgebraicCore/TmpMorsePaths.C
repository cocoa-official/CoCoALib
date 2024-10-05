//   Copyright (c)  2015  John Abbott,  Anna M. Bigatti
//   Original author: 2015  Mario Albert

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


#include "CoCoA/TmpMorsePaths.H"

using std::make_pair;

namespace CoCoA
{
  namespace Involutive
  {

    void MorsePaths::myAddPath(const ConstResIter& m, const RingElem& elem)
    {
      CoCoA_ASSERT(!IsZero(elem));
      const PathMap::iterator iter(myPaths.find(m));
      if (iter != myPaths.end())
      {
        // iter->second += elem;
        owner(elem)->myAdd(raw(iter->second), raw(iter->second), raw(elem));
        if (IsZero(iter->second))
        {
          myPaths.erase(iter);
        }
      } else {
        myPaths.insert(make_pair(m, elem));
      }
    }


    bool MorseIterCompare::operator()(const ConstResIter& lhs, const ConstResIter& rhs) const
    {
      return (lhs->first) < (rhs->first);
    }


    std::string MorsePaths::toStr() const
    {
      std::ostringstream os;
      for (PathMap::const_iterator i = myPaths.begin(); i != myPaths.end(); ++i)
      {
        os << (i->first)->first << "--- " << i->second << " --->" << std::endl;
      }
      return os.str();
    }

  } // end of namespace Involutive
} // end of namespace CoCoA
