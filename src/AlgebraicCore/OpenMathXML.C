//   Copyright (c)  2005,2011  John Abbott and Anna M. Bigatti

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


#include "CoCoA/OpenMathXML.H"
#include "CoCoA/error.H"
#include "CoCoA/assert.H"
#include "CoCoA/convert.H"
#include "CoCoA/utils.H"

#include <iostream>
using std::ostream;
using std::istream;
#include <vector>
using std::vector;
//#include <string>
using std::string;
#include "gmp.h"


namespace CoCoA
{

    OpenMathOutputXML::OpenMathOutputXML(ostream& out):
        myOut(out),
        myLevel(0),
        myTagStack()
    {}


    OpenMathOutputXML::~OpenMathOutputXML()
    {
      CoCoA_ASSERT(myLevel == 0); //??? how to avoid throwing inside dtor???
      // ??? Should send token to say "over and out"???
    }


    inline const char* OpenMathOutputXML::myIndent()
    {
      static const int SpacesPerLevel = 2;
      static const int MaxIndent = 20;
      static const string spaces(SpacesPerLevel*MaxIndent, ' ');
      static const char* const ptr = spaces.c_str();
      if (myLevel >= MaxIndent) return ptr;
      return ptr+(MaxIndent-myLevel)*SpacesPerLevel;
    }


    void OpenMathOutputXML::mySend(const MachineInt& n)
    {
      myOut << myIndent() << "<OMI> " << n << " </OMI>\n";
    }


    void OpenMathOutputXML::mySend(const BigInt& N)
    {
      myOut << myIndent() << "<OMI> " << N << " </OMI>\n";
    }

    void OpenMathOutputXML::mySend(const OpenMathSymbol& s)
    {
      myOut << myIndent() << "<OMS cd=\"" << CD(s) << "\" name=\"" << name(s) << "\" />\n";
    }


    void OpenMathOutputXML::mySendApplyStart()
    {
//       std::clog<<"ApplyStart at level " << myLevel << endl;
//       CoCoA_ASSERT(myLevel > 0);
      myOut << myIndent() << "<OMA>\n";
      myTagStack.push(OpenMathTag::Apply);
      ++myLevel;
    }

    void OpenMathOutputXML::mySendApplyEnd()
    {
//       CoCoA_ASSERT(myLevel > 1);
      --myLevel;
      CoCoA_ASSERT(!myTagStack.empty());
      CoCoA_ASSERT(myTagStack.top() == OpenMathTag::Apply);
      myTagStack.pop();
      myOut << myIndent() << "</OMA>\n";
    }


    void OpenMathOutputXML::mySendObjectStart()
    {
      CoCoA_ASSERT(myLevel == 0);
      myOut << "<OMOBJ>\n";
      myTagStack.push(OpenMathTag::Obj);
      ++myLevel;
    }

    void OpenMathOutputXML::mySendObjectEnd()
    {
      CoCoA_ASSERT(myLevel == 1);
      --myLevel;
      CoCoA_ASSERT(!myTagStack.empty());
      CoCoA_ASSERT(myTagStack.top() == OpenMathTag::Obj);
      myTagStack.pop();
      myOut << "</OMOBJ>\n";
    }



    OpenMathInputXML::OpenMathInputXML(istream& in):
        myReadStatus(AlreadyRead),
        myCurrTagType(OpenMathTag::ReachedEOF), // could be any value
        myIntValue(-1),     // could be any value
        mySymbol(),         // default "unset" OpenMathSymbol
        myIn(in),           // reference to input stream
        myLevel(0),         // nesting level initially 0
        myTagStack()        // stack initially empty
    {
      myIn.unsetf(std::ios::skipws);
    }


    OpenMathInputXML::~OpenMathInputXML()
    {
      CoCoA_ASSERT(myLevel == 0);
    }


    void OpenMathInputXML::advance()
    {
      if (myReadStatus == NotYetRead)
        ReadNextNode(); //??? discard next node???
      myReadStatus = NotYetRead;
    }


    OpenMathTag OpenMathInputXML::myCurrTag()
    {
      if (myReadStatus == NotYetRead)
        ReadNextNode();
      return myCurrTagType;
    }


    long OpenMathInputXML::NumDescendants() const
    { return 0; //????????
    }


    bool OpenMathInputXML::myRecv(long & n)
    {
      if (myCurrTag() != OpenMathTag::Int) return false;
//???        CoCoA_THROW_ERROR1("OpenMath node is not OMI");
      if (!IsConvertible(n, myIntValue))
        CoCoA_THROW_ERROR2(ERR::ArgTooBig, "OpenMathInputXML::InputInteger");
      return true;
    }


//???????
//     bool OpenMathInputXML::myRecv(unsigned long & n)
//     {
//       if (myCurrTag() != OpenMathInt) return false;
// //???        CoCoA_THROW_ERROR1("OpenMath node is not OMI");
//       if (!IsConvertible(n, myIntegerValue))
//         CoCoA_THROW_ERROR2(ERR::ArgTooBig, "OpenMathInputXML::InputInteger");
//       return true;
//     }


    bool OpenMathInputXML::myRecv(BigInt& N)
    {
      if (myCurrTag() != OpenMathTag::Int) return false;
//???        CoCoA_THROW_ERROR1("OpenMath node is not OMI");
      N = myIntValue;
      return true;
    }


    bool OpenMathInputXML::myRecv(OpenMathSymbol& s)
    {
      if (myCurrTag() != OpenMathTag::Sym) return false;
//???        CoCoA_THROW_ERROR1("OpenMath node is not OMS");
      s = mySymbol;
      return true;
 }

    char OpenMathInputXML::ReadChar()
    {
//???      if (myIn.eof()) return '\0';
      std::clog<<'['<<char(myIn.peek())<<']';
      return myIn.get(); //???? what about EOF????  BUG BUG BUG ???
    }

    char OpenMathInputXML::SkipWSReadChar()
    {
      char ch;
      myIn >> ch;
      if (!myIn.eof()) return ch;
      return '\0';
      // while (!myIn.eof())
      // {
      //   char ch = ReadChar();
      //   if (ch == ' ' || ch == '\t' || ch == '\n') continue;
      //   return ch;
      // }
      // return '\0';
    }


    // NB a space in `expected' is regarded as any amount of whitespace.
    bool OpenMathInputXML::SkipMatch(const string& expected)
    {
      if (myIn.eof()) return false;
      const long nchars = len(expected);
      for (long i = 0; i < nchars; ++i)
      {
        const char expected_ch = expected[i];
        if (expected_ch == ' ') { myIn >> std::ws; continue; }
        char ch = ReadChar();
//        std::clog<<"EXPECTING `"<<expected_ch<<"'   READ `"<<ch<<"'"<<std::endl;
        if (ch != expected_ch)
          return false;
      }
      return true;
    }

    bool OpenMathInputXML::ReadDecimalString(string& DecimalDigits)
    {
      myIn >> std::ws;
      if (myIn.eof()) return false; // hit EOF

      DecimalDigits.clear();
      // Check to see if number is negative; grab sign if so
      if (myIn.peek() == '-')
      {
        DecimalDigits += ReadChar();
        //??? myIn >> ws; // allow whitespace between sign first digit
      }
      bool ReadAtLeastOneDigit = false;
      while (!myIn.eof() && isdigit(myIn.peek()))
      {
        ReadAtLeastOneDigit = true;
        DecimalDigits += ReadChar();
      }
      return ReadAtLeastOneDigit;
    }


    bool OpenMathInputXML::ReadStringInQuotes(string& QuotedString)
    {
      if (!SkipMatch(" \"")) return false;
      while (!myIn.eof() && myIn.peek() != '"')
      {
        QuotedString += ReadChar();
      }
      return SkipMatch("\"");
    }


    void OpenMathInputXML::ReadNextNode()
    {
      myIn >> std::ws;
      std::clog<<"ReadNextNode: looking at `"<<char(myIn.peek())<<"'"<<std::endl;
      if (myIn.eof()) { myCurrTagType = OpenMathTag::ReachedEOF; return; }
      if (ReadChar() != '<'){std::clog<<"DIDN'T FIND `<' WHERE ONE SHOULD BE"<<std::endl;return;}
      myIn >> std::ws;
      if (myIn.peek() == '/')
      {
        if (!SkipMatch("/ OMA > ")){std::clog<<"DIDN'T FIND `</OMA>' AFTER READING `</'"<<std::endl;return;}
        if (myTagStack.empty() || myTagStack.top() != OpenMathTag::Apply){std::clog<<"MISPLACED OR TOO MANY `</OMA>'s"<<std::endl;return;}
        myCurrTagType = OpenMathTag::Apply;
        --myLevel;
      }
      if (!SkipMatch(" OM")){std::clog<<"DIDN'T FIND `<OM' WHERE ONE SHOULD BE"<<std::endl;return;}
      char ch = ReadChar();
      if (myIn.eof()){std::clog<<"UNEXPECTED EOF"<<std::endl;return;}
      switch (ch)
      {
      case 'A':
      {
        myCurrTagType = OpenMathTag::Apply;
        myTagStack.push(OpenMathTag::Apply);
        if (!SkipMatch(" > ")) {std::clog<<"`OMA' not followed by `>': instead found '"<<ReadChar()<<"'"<<std::endl;}
        ++myLevel;
        return;
      }
      case 'I':
      {
        if (!SkipMatch(" >")) {std::clog<<"`OMI' not followed by `>': instead found '"<<ReadChar()<<"'"<<std::endl;}
        std::clog<<"READING INTEGER...";
        myCurrTagType = OpenMathTag::Int;
        string decimal;
        if (!ReadDecimalString(decimal)) {std::clog<<"BAD DECIMAL NUMBER"<<std::endl; return;}

        std::clog<<"DECIMAL string is `"<<decimal<<"'"<<std::endl;
        mpz_set_str(mpzref(myIntValue), decimal.c_str(), 10);
        myIntValue = BigIntFromString(decimal);

        if (!SkipMatch(" < /OMI > ")) { std::clog<<"DIDN'T FIND '</OMI>' WHERE IT SHOULD BE"<<std::endl; return;}
        return;
      }
      case 'S':
      {
        myCurrTagType = OpenMathTag::Sym;
        std::clog<<"READING OMS...";
        SkipMatch(" cd = ");
        string cd;
        ReadStringInQuotes(cd);
        std::clog<<"CD part is `"<<cd<<"'   ";
        SkipMatch(" name = ");
        string name;
        ReadStringInQuotes(name);
        std::clog<<"NAME part is `"<<name<<"'"<<std::endl;
        mySymbol = OpenMathSymbol(cd, name);
        SkipMatch(" /> ");
        return;
      }
      default:
        std::clog<<"UNKNOWN NODE TYPE: OM" << ch << std::endl;
        return;
      }
    }

} // end of namespace CoCoA
