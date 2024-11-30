#ifndef CoCoA_LogStream_H
#define CoCoA_LogStream_H

//   Copyright (c)  2017  John Abbott,  Anna M. Bigatti

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

#include <iosfwd>
//using std::ostream;

namespace CoCoA
{

  class LogStreamForThisBlock
  {
  public:
    LogStreamForThisBlock(std::ostream& out): myPrevLogStreamPtr(ourActiveLogStreamPtr) { mySetLogStreamPtr(&out); }
    ~LogStreamForThisBlock() { mySetLogStreamPtr(myPrevLogStreamPtr); }
  private:
    void mySetLogStreamPtr(std::ostream* out) { ourActiveLogStreamPtr = out; }
  private: // data member, static data member
    std::ostream* myPrevLogStreamPtr;
    static std::ostream* ourActiveLogStreamPtr;
    friend std::ostream& LogStream();
  };


  std::ostream& operator<<(std::ostream& out, const LogStreamForThisBlock& /*dummy*/);
  
  inline std::ostream& LogStream()
  {
    return *LogStreamForThisBlock::ourActiveLogStreamPtr;
  }

} // end of namespace CoCoA

#endif
