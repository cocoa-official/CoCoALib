#ifndef CoCoA_time_H
#define CoCoA_time_H

//   Copyright (c)  2005,2020  John Abbott, and Anna M. Bigatti

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


namespace CoCoA
{

  double CpuTime() noexcept;  ///< Returns user CPU usage in seconds since start of program.
  double ElapsedTime() noexcept; ///< Number of seconds since process start

  double RealTime() noexcept; ///< Returns number of seconds elapsed since some fixed point in the past.
  void DateTime(long& date, long& time) noexcept; ///<  yyyy*10000+mm*100+dd,  hh*10000+mm*100+ss

} // end of namespace CoCoA

#endif
