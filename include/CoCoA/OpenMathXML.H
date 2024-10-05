#ifndef CoCoA_OpenMathXML_H
#define CoCoA_OpenMathXML_H

//   Copyright (c)  2005  John Abbott,  Anna M. Bigatti

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


#include "CoCoA/OpenMath.H"
#include "CoCoA/BigInt.H"

#include <iostream>
// using std::ostream;
#include <stack>
// using std::stack;
#include <string>
// using std::string


namespace CoCoA
{

  class OpenMathOutputXML: public OpenMathOutputBase
  {
  public:
    OpenMathOutputXML(std::ostream& out);
    ~OpenMathOutputXML();
    void mySend(const MachineInt& n) override;
    void mySend(const BigInt& N) override;
    void mySend(const OpenMathSymbol& s) override;
    void mySendApplyStart() override;
    void mySendApplyEnd() override;
    void mySendObjectStart() override;
    void mySendObjectEnd() override;
    //void OutputAttribute();//??????
  private:
    std::ostream& myOut;  // UGLY, RISKY!!???!!
    long myLevel;  // nesting level
    std::stack<OpenMathTag> myTagStack;  // only really needed for debugging
    const char* myIndent();
  };


  class OpenMathInputXML: public OpenMathInputBase
  {
  public:
    OpenMathInputXML(std::istream& in); // just sets ref count to zero.
    virtual ~OpenMathInputXML();

    void advance() override;
    OpenMathTag myCurrTag() override;
    long NumDescendants() const override;
    bool myRecv(long & n) override; // SHORT INTEGER
//????    virtual bool myRecv(unsigned long & n); // UNSIGNED SHORT INTEGER
    bool myRecv(BigInt& N) override;    // LONG INTEGER
    bool myRecv(OpenMathSymbol& s) override;
  private: // data members
    enum class ReadStatus { NotYetRead, AlreadyRead };
    static constexpr ReadStatus NotYetRead = ReadStatus::NotYetRead;
    static constexpr ReadStatus AlreadyRead = ReadStatus::AlreadyRead;
    ReadStatus myReadStatus;
    OpenMathTag myCurrTagType;
    BigInt myIntValue;
    OpenMathSymbol mySymbol;
    std::istream& myIn;  // UGLY, RISKY!!???!!
    long myLevel;  // nesting level
    std::stack<OpenMathTag> myTagStack;  // only really needed for debugging
  private:
    char ReadChar();
    char SkipWSReadChar();
    void ReadNextNode();
    bool SkipMatch(const std::string& expected);
    bool ReadDecimalString(std::string& DecimalDigits);
    bool ReadStringInQuotes(std::string& QuotedString);
  };


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/OpenMathXML.H,v 1.12 2022/02/18 14:11:47 abbott Exp $
// $Log: OpenMathXML.H,v $
// Revision 1.12  2022/02/18 14:11:47  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.11  2021/10/30 17:09:00  abbott
// Summary: Used keyword override (redmine 1625)
//
// Revision 1.10  2021/03/03 22:09:32  abbott
// Summary: New enum class (redmine 894)
//
// Revision 1.9  2021/01/07 14:19:14  abbott
// Summary: Corrected copyright
//
// Revision 1.8  2011/11/09 13:47:56  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.7  2011/08/14 15:52:18  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.6  2011/08/12 15:49:28  abbott
// Added mem fns for BigInt; also new data member.
//
// Revision 1.5  2011/03/11 14:49:08  abbott
// Changed size_t into long.
//
// Revision 1.4  2010/03/18 13:53:48  abbott
// Minor rationalization to OpenMath implementation: moved op<<(OMOut,ZZ) to ZZ files.
//
// Revision 1.3  2008/12/16 21:10:32  abbott
// Replaced the various output fns for different sort of machine integers by a
// single one for MachineInt.
//
// Revision 1.2  2007/10/30 17:14:12  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.7  2006/12/06 17:14:12  cocoa
// -- removed #include "config.H"
//
// Revision 1.6  2006/11/27 14:22:04  cocoa
// -- ooops, missing semicolon
//
// Revision 1.5  2006/11/27 13:53:10  cocoa
// changed mySend(std::size_t n)  into  mySend(unsigned long n)
//
// Revision 1.4  2006/11/27 13:04:35  cocoa
// Added explicit OpenMath output operators for all integer types.
//
// Revision 1.3  2006/11/02 13:25:44  cocoa
// Simplification of header files: the OpenMath classes have been renamed.
// Many minor consequential changes.
//
// Revision 1.2  2006/10/06 14:04:15  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//

#endif
