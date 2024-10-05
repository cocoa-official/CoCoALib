//   Copyright (c) 2009  John Abbott,  Anna Bigatti
//   Orig author: 2009  Giovanni Lagorio (lagorio@disi.unige.it)
//
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

#include "CompilationDate.H"

using std::string;

namespace CoCoA
{

  namespace // anon
  {
    string ReorderDate(const char* const USDate)
    {
      string date(11, ' ');
      date[0] = USDate[4];
      date[1] = USDate[5];
      date[3] = USDate[0];
      date[4] = USDate[1];
      date[5] = USDate[2];
      date[7] = USDate[7];
      date[8] = USDate[8];
      date[9] = USDate[9];
      date[10] = USDate[10];
      return date;
    }
  }

  const string& CompilationDate()
  {
    // Unfortunately the format of __DATE__ is broken (american) so we must reorder it
    const static string DateTime = ReorderDate(__DATE__) + " at " + string(__TIME__);
    return DateTime;
  }

} // end of namespace CoCoA

