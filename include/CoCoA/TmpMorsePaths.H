#ifndef CoCoA_TmpMorsePaths_H
#define CoCoA_TmpMorsePaths_H

//   Copyright (c) 2015  John Abbott,  Anna M. Bigatti
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

#include "CoCoA/TmpMorseElement.H"
#include <sstream>

namespace CoCoA
{
  namespace Involutive
  {

    /* This file defines the class MorsePaths and all necessary stuff for it.
     * A MorsePath is an edge in the corresponding MorseGraph. It consists
     * of a starting point, end point (both MorseElements) and a value (a RingElem)
     * We using an usual data structure:
     * A Path is accessible by its end point that is, we store for an MorseElement m
     * a list of paths where m is the end point. The class MorsePaths shall
     * represent this list. We stress out, that an instance of MorsePaths does not know
     * anything of m. These information is placed in the MorseGraph itself.
     */
    class MorsePaths;

    // This is an element in the MorseGraph
    typedef std::map<MorseElement, MorsePaths>::const_iterator ConstResIter;

    // the paths are stored in an ordered list. This struct defines this ordering
    struct MorseIterCompare
    {
      // The ordering is simply the ordering of the corresponding MorseElements
      bool operator()(const ConstResIter& lhs, const ConstResIter& rhs) const;
    };

    // This is one 'path' it maps the Element ConstResIter with RingElem to the
    // current MorseElement
    // We do not know in MorsePaths which is the current element this Element is the key
    // Element of a key-value pair, where the vlaue is MorsePaths
    typedef std::map<ConstResIter, RingElem, MorseIterCompare> PathMap;

    /**
     * This file defines the class MorsePaths and all necessary stuff for it.
     * A MorsePath is an edge in the corresponding MorseGraph. It consists
     * of a starting point, end point (both MorseElements) and a value (a RingElem)
     * We using an usual data structure:
     * A Path is accessible by its end point that is, we store for an MorseElement m
     * a list of paths where m is the end point. The class MorsePaths shall
     * represent this list. We stress out, that an instance of MorsePaths does not know
     * anything of m. These information is placed in the MorseGraph itself.
     */
    class MorsePaths
    {
    public:
      /**
       * Returns a constant reference to myPaths
       *
       * @return a constant reference to myPaths
       */
      inline const PathMap& myGetPaths() const { return myPaths; }

      /**
       * Add to the path with origin m the RingElem elem.
       * If there is no path from m, we add a new path.
       *
       * @param m
       * @param elem
       */
      void myAddPath(const ConstResIter& m, const RingElem& elem);

      /**
       * Test if there are any paths.
       *
       * @return returns true if no paths are present
       */
      inline bool IamEmpty() const { return myPaths.empty(); }

      // Return the path from MorseElement m.
      /**
       * Return the path from m. If there is no path it returns zero which is
       * in the polynomial ring of m.
       *
       * @param m
       * @return a RingElem
       */
      inline ConstRefRingElem myGetPath(const ConstResIter& m)
      {
        // find the corresponding path
        const PathMap::iterator iter(myPaths.find(m));
        // It does not exist a corresponding path -> return zero
        if (iter == myPaths.end())
        {
          return zero(owner((m->first).myGetBasisElement()));
        }
        return iter->second;
      }

      /**
       * Returns the MorsePaths as string
       *
       * @return a string
       */
      std::string toStr() const;
    private:
      //basic datastructure
      PathMap myPaths;

    };

    /**
     * MorsePaths to stream
     *
     * @param os
     * @param obj
     *
     * @return stream object
     */
    inline std::ostream& operator<<(std::ostream& os, const MorsePaths& obj)
    {
      os << obj.toStr();
      return os;
    }


  } // end of namespace Involutive
} // end of namespace CoCoA

#endif
