//   Copyright (c)  2010,2017  Anna Bigatti, John Abbott
//   Main author: Anna M Bigatti

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

#include "OnlineHelp.H"

#include "CoCoA/error.H"

#include <algorithm>
using std::transform;
#include <cctype>
//using std::tolower;
#include <fstream>
using std::ifstream;
#include <iostream>
using std::endl;
#include <sstream> // for ostringstream - nasty trick
#include <string>
using std::string;
#include <vector>
using std::vector;


extern string packageDir;

namespace CoCoA
{
namespace OnlineHelp
{

  namespace // anonymous
  {

    std::ifstream OpenXMLManual()
    {
      ifstream in(CoCoAManFileName());
      if (!in)
        CoCoA_THROW_ERROR("CoCoA-5 manual missing (or not readable); it should be in "+CoCoAManFileName(), "OnlineHelp::OpenXMLManual()");
      return in;
    }
    
  } // end of namespace anonymous

  class index
  {
  public:
    class entry
    {
    public:
      entry(const std::string& title, const std::string& inFileName):
        myTitleValue(title), myFileNameValue(inFileName) {}
      // ~entry(); // default is OK
      const std::string& myTitle() const {return myTitleValue;}
      const std::string& myFileName() const {return myFileNameValue;}
      const std::vector<std::string>& myKeys() const {return myKeysValue;}
      const std::vector<std::string>& myTypes() const {return myTypesValue;}
      const std::vector<std::string>& myRtnTypes() const {return myRtnTypesValue;}
      const std::vector<std::string>& mySeeAlso() const {return mySeeAlsoValue;}
      const std::string& mySyntax() const {return mySyntaxValue;}
      void myAddKey(const std::string& key);
      void myAddType(const std::string& key);
      void myAddRtnType(const std::string& key);
      void myAddSee(const std::string& key);
      void myAddKeys(std::ifstream& in);
      void myAddTypes(std::ifstream& in);
      void myAddSeeAlso(std::ifstream& in);
      void myAddSyntax(std::ifstream& in);

    private: // data members
      std::string myTitleValue;
      std::string myFileNameValue;
      std::vector<std::string> myKeysValue;
      std::vector<std::string> myTypesValue;
      std::vector<std::string> myRtnTypesValue;
      std::vector<std::string> mySeeAlsoValue;
      std::string mySyntaxValue;
    };

  public:
    index();
    index& UniqueIndex(void); // pseudoconstructor /// 2020-01-31
    // ~index(); // default is OK
    void myLoad(std::ostream &out);
    void myLoadExtraXML(std::ostream &out, const std::string& inFileName);
    const index::entry& myEntry(std::size_t i) const {return myVec[i];}
    const std::string& myTitle(std::size_t i) const {return myVec[i].myTitle();}
    const std::string& myFileName(std::size_t i) const {return myVec[i].myFileName();}
    const std::vector<std::string>& myKeys(std::size_t i) const {return myVec[i].myKeys();}
    const std::vector<std::string>& myTypes(std::size_t i) const {return myVec[i].myTypes();}
    const std::vector<std::string>& myRtnTypes(std::size_t i) const {return myVec[i].myRtnTypes();}
    const std::vector<std::string>& mySeeAlso(std::size_t i) const {return myVec[i].mySeeAlso();}
    const std::string& mySyntax(std::size_t i) const {return myVec[i].mySyntax();}
    std::size_t mySize() const {return myVec.size();}
    //  std::string myManCoCoAVersion() const {return myManCoCoAVersionValue;}
    std::string myManDate() const {return myManDateValue;}

  private:
    void myAddEntry(std::ifstream& in, const std::string& title, const std::string& ClosingTag, const std::string& inFileName);

  private: // data members
    std::vector<entry> myVec;
    std::string myManDateValue;
    //    std::string myManCoCoAVersionValue;
  };

  // ----------------------------------------------------------------------
  // fwd decl -- defined later in this file
  std::string LowerCase(const std::string& str);
  bool IsStr1BeforeStr2(std::ifstream& in, const std::string& s1, const std::string& s2);
  std::string StringInTag(const std::string& line, const std::string& XMLTag);
  std::string SkipToStringInTag(std::ifstream& in, const std::string& XMLTag);
  std::string CleanupLine(const std::string& line);
  std::string CleanupKeyword(const std::string& key);
  void PrintCommAndFunc(std::ostream &out, std::string& line);
  void PrintCommAndFuncRtn(std::ostream &out, std::string& line);
  void PrintKeys(std::ostream &out, const std::string& s);
  void PrintSeeAlso(std::ostream &out, const std::string& s);
  void PrintSyntax(std::ostream &out, const std::string& s);
  void PrintAllMatchesFor(std::ostream &out, const std::string& str);
  void PrintAllMatchesFor(std::ostream &out, const std::string& str,
                          const std::vector<long>& found);

  // ----------------------------------------------------------------------
  inline bool IsSubstring(const std::string& s, const std::string& line)
  { return line.find(s) != string::npos; }

  // ----------------------------------------------------------------------

  //-- index --

  //  const index& UniqueIndex(); // no, can be updated
  index& UniqueIndex();

  index::index()
  {
    std::ostringstream os; // does not print
    myLoad(os);
  }


  void index::myLoad(std::ostream &out)
  {
    out << "Loading CoCoAManual index ..." << std::flush;
    ifstream in = OpenXMLManual();
///    ifstream in(CoCoAManFileName());
///    if (!in)
///      CoCoA_THROW_ERROR("Cannot find XML file "+CoCoAManFileName(), "OnlineHelp::myLoad()");
//     string line; // used only inside while loop
//     string LastCommand;
//     string NewCommand = "a";
//     //    string ManCoCoAVersion;
//     string ManDate;
    myVec.clear();
    myLoadExtraXML(out, CoCoAManFileName());
    out << "... CoCoAManual index loaded" << endl;
  }


  void index::myLoadExtraXML(std::ostream &out, const std::string& inFileName)
  {
    out << "Loading manual file " << inFileName << "..." << std::flush;
    ifstream in(inFileName.c_str());
    if (!in)
      CoCoA_THROW_ERROR("Cannot find input file "+inFileName, "myLoadExtraXML");
    string line; // used only inside while loop
    string LastCommand;
    string NewCommand = "a";
    //    string ManCoCoAVersion;
    string ManDate;
    //myVec.clear();
    while (!in.eof())
    {
      getline(in, line);
      if (IsSubstring("<!--", line)) continue;
      if (IsSubstring("<command>", line))
      {
        LastCommand = NewCommand;
        NewCommand = SkipToStringInTag(in, "<title>");
        // check for proper sorting
        if (LowerCase(LastCommand).compare(LowerCase(NewCommand))>0)
          std::cout << "OnlineHelp: unsorted entries: "
                    << LastCommand << " -- " << NewCommand << std::endl;
        myAddEntry(in, CleanupLine(NewCommand), "</command>", inFileName);
      }
      if (IsSubstring("<section>", line))
        myAddEntry(in, CleanupLine(SkipToStringInTag(in, "<title>")), "</section>", inFileName);
    }
    out << "... loaded" << endl;
  }


  void index::myAddEntry(std::ifstream& in, const std::string& title, const std::string& ClosingTag, const std::string& inFileName)
  {
    myVec.push_back(entry(title, inFileName));
    entry& e(myVec.back());
    e.myAddKey(title);
    string line;
    getline(in, line);
    while (!in.eof() && !IsSubstring(ClosingTag, line))
    {
      if (IsSubstring("<keys>", line)) e.myAddKeys(in);
      if (IsSubstring("<types>", line)) e.myAddTypes(in);
      if (IsSubstring("<seealso>", line)) e.myAddSeeAlso(in);
      if (IsSubstring("<syntax>", line)) e.myAddSyntax(in);
      getline(in, line);
    }
  }


  //-- index::entry --
  void index::entry::myAddKey(const std::string& key)
  { myKeysValue.push_back(CleanupLine(LowerCase(key))); }

  void index::entry::myAddType(const std::string& t)
  { myTypesValue.push_back(t); }

  void index::entry::myAddRtnType(const std::string& t)
  { myRtnTypesValue.push_back(t); }

  void index::entry::myAddSee(const std::string& t)
  { mySeeAlsoValue.push_back(t); }

  void index::entry::myAddKeys(std::ifstream& in)
  {
    string line;
    getline(in, line);
    string s; // used only inside while loop
    while ((!in.eof()) && (!IsSubstring("</keys>",line)))
    {
      if ((s=StringInTag(line, "<key>")) != "") myAddKey(s);
      getline(in, line);
    }
  }

  void index::entry::myAddTypes(std::ifstream& in)
  {
    string line;
    getline(in, line);
    string s; // used only inside while loop
    while ((!in.eof()) && (!IsSubstring("</types>",line)))
    {
      if ((s=StringInTag(line, "<type>")) != "") myAddType(s);
      getline(in, line);
    }
  }


  void index::entry::myAddSeeAlso(std::ifstream& in)
  {
    string line;
    getline(in, line);
    string s; // used only inside while loop
    while ((!in.eof()) && (!IsSubstring("</seealso>",line)))
    {
      if ((s=StringInTag(line, "<see>")) != "") myAddSee(s);
      getline(in, line);
    }
  }


  void index::entry::myAddSyntax(std::ifstream& in)
  {
    string line;
    getline(in, line);
    while ((!in.eof()) && (!IsSubstring("</syntax>",line)))
    {
      if (IsSubstring("<type>", line))
        myAddType(StringInTag(line, "<type>"));
      if (IsSubstring("<rtn>", line))
        myAddRtnType(StringInTag(line, "<rtn>"));
      mySyntaxValue += "--> ";
      mySyntaxValue += CleanupLine(line);
      mySyntaxValue += "\n";
      getline(in, line);
    }
  }


  //----------------------------------------------------------------

  std::string LowerCase(const std::string& str)
  {
    string s(str); // slightly wasteful, but readable
    transform(s.begin(), s.end(), s.begin(), tolower);
    return s;
  }


  bool IsStr1BeforeStr2(std::ifstream& in, const std::string& s1, const std::string& s2)
  {
    string line;
    getline(in, line);
    while (!in.eof())
    {
      if ( IsSubstring(s2, line) ) return false;
      if ( IsSubstring(s1, line) ) return true;
      getline(in, line);
    }
    CoCoA_THROW_ERROR(CoCoA::ERR::ShouldNeverGetHere, "IsStr1BeforeStr2");
    return false;
  }


  // only for opening and closing tag in the same line
  std::string StringInTag(const std::string& line, const std::string& XMLTag)
  {
    size_t open;
    if ( (open=line.find(XMLTag)) == string::npos) return "";

    string ClosedXMLTag = XMLTag;
    ClosedXMLTag.replace(0, 1, "</");
    size_t close = line.find(ClosedXMLTag);
    if (close == string::npos)
      CoCoA_THROW_ERROR(XMLTag+" closing tag not found in this line", "StringInTag");
    size_t StartPos = open + XMLTag.length();
    return line.substr(StartPos, close-StartPos);
  }


  std::string SkipToStringInTag(std::ifstream& in, const std::string& XMLTag)
  {
    string line;
    getline(in, line);
    string s;
    while (!in.eof())
    {
      if ( (s=StringInTag(line, XMLTag)) != "") return s;
      getline(in, line);
    }
    return "";
  }


  void ReplaceWith(std::string& line, const std::string& SFrom, const std::string& STo)
  {
    size_t open;
    while ( (open=line.find(SFrom)) != string::npos)
      line.replace(open, SFrom.length(), STo);
  }


  void ReplaceWithQuotes(std::string& line, const std::string& s)
  { ReplaceWith(line, s, "\""); }


  std::string& ProcessLine(std::string& line)
  {
    ReplaceWithQuotes(line, "<quotes>"); ReplaceWithQuotes(line, "</quotes>");
    ReplaceWithQuotes(line, "<tt>");     ReplaceWithQuotes(line, "</tt>");
    ReplaceWithQuotes(line, "<code>");   ReplaceWithQuotes(line, "</code>");
    ReplaceWithQuotes(line, "<ttref>");  ReplaceWithQuotes(line, "</ttref>");
    ReplaceWithQuotes(line, "<ref>");    ReplaceWithQuotes(line, "</ref>");
    ReplaceWith(line, "<em>", "*");      ReplaceWith(line, "</em>", "*");
    ReplaceWith(line, "<b>", "*");       ReplaceWith(line, "</b>", "*");
    ReplaceWith(line, "<i>", "");        ReplaceWith(line, "</i>", "");
    ReplaceWith(line, "<verbatim>", ""); ReplaceWith(line, "</verbatim>", "");
    ReplaceWith(line, "<formula>", " "); ReplaceWith(line, "</formula>", " ");
    ReplaceWith(line, "&lt;", "<");
    ReplaceWith(line, "&gt;", ">");
    ReplaceWith(line, "&apos;", "'");
    ReplaceWith(line, "&amp;", "&");
    ReplaceWith(line, "<less_eq/>", "<=");
    ReplaceWith(line, "<par/>", "");    ReplaceWith(line, "<par />", "");
    ReplaceWith(line, "<ie/>", "i.e."); ReplaceWith(line, "<ie />", "i.e.");
    ReplaceWith(line, "<eg/>", "e.g."); ReplaceWith(line, "<eg />", "e.g.");
    ReplaceWith(line, "<BOOL/>", "BOOL");
    ReplaceWith(line, "<ERROR/>", "ERROR");
    ReplaceWith(line, "<FUNCTION/>", "FUNCTION");
    ReplaceWith(line, "<IDEAL/>", "IDEAL");
    ReplaceWith(line, "<INT/>", "INT");
    ReplaceWith(line, "<LIST/>", "LIST");
    ReplaceWith(line, "<MAT/>", "MAT");
    ReplaceWith(line, "<MODULE/>", "MODULE");
    ReplaceWith(line, "<MODULEELEM/>", "MODULEELEM");
    ReplaceWith(line, "<RAT/>", "RAT");
    ReplaceWith(line, "<RECORD/>", "RECORD");
    ReplaceWith(line, "<RING/>", "RING");
    ReplaceWith(line, "<RINGELEM/>", "RINGELEM");
    ReplaceWith(line, "<RINGHOM/>", "RINGHOM");
    ReplaceWith(line, "<STRING/>", "STRING");
    ReplaceWith(line, "<TAGGED/>", "TAGGED");
    ReplaceWith(line, "<TYPE/>", "TYPE");

    //ReplaceWith(line, "<cocoa_version/>", UniqueIndex().myManCoCoAVersion());
    //    ReplaceWith(line, "<cocoa_version/>", CoCoAVersion());
    //    ReplaceWith(line, "<cocoa_date/>", UniqueIndex().myManDate());
    return line;
  }


  std::string CleanupLine(const std::string& line)
  {
    // CleanupLine is called by the ctor UniqueIndex
    string s(line);
    ReplaceWithQuotes(s, "<quotes>"); ReplaceWithQuotes(s, "</quotes>");
    ReplaceWithQuotes(s, "<tt>");     ReplaceWithQuotes(s, "</tt>");
    ReplaceWith(s, "<type>", "");     ReplaceWith(s, "</type>", "");
    ReplaceWith(s, "<rtn>", "");      ReplaceWith(s, "</rtn>", "");
    ReplaceWith(s, "<em>", "*");      ReplaceWith(s, "</em>", "*");
    ReplaceWith(s, "&lt;", "<");
    ReplaceWith(s, "&gt;", ">");
    ReplaceWith(s, "&apos;", "'");
    ReplaceWith(s, "&amp;", "&");
    ReplaceWith(s, "<less_eq/>", "<=");
    return s;
  }


  void SkipComment(std::ifstream& in, std::string& line)
  {
    while (  (!in.eof()) && !IsSubstring("-->", line) )  getline(in, line);
    if ( in.eof() )
      CoCoA_THROW_ERROR("missing comment close tag: eof", "SkipComment");
    const size_t close = line.find("-->");
    if (close != line.size()-3)
      CoCoA_THROW_ERROR("text after closed comment: \""+line+"\"", "SkipComment");
  }


  void PrintExample(std::ostream &out, std::ifstream& in)
  {
    out << "----<<<<  example  >>>>----" << endl;
    string line;
    getline(in, line);
    while ( !IsSubstring("</example>", line) )
    {
      out << CleanupLine(line) << endl;
      getline(in, line);
    }
    out << "----==< end example >==----" << endl;
  }


  void PrintExampleWithoutOutput(std::ostream &out, std::ifstream& in)
  {
    out << "-----<<<  example  >>>-----" << endl;
    out << "use ZZ/(4)[SillyName];  R := \"clear R\";  P := \"clear P\";" << endl;
    string line;
    getline(in, line);
    while ( !IsSubstring("</example>", line) )
    {
      if (IsSubstring("/**/", line))  out << CleanupLine(line) << endl;
      getline(in, line);
    }
    out << "----==< end example >==----" << endl;
  }


  void PrintSyntax(std::ostream &out, const std::string& s)
  {
    const index& index(UniqueIndex());
    bool IsFirst=true;
    for (size_t i=0; i<index.mySize(); ++i)
      if ( LowerCase(index.myTitle(i)) == LowerCase(s)
           && !index.mySyntax(i).empty())
      {
        if (IsFirst)
        {
          IsFirst = false;
          //          out << "--====<  syntax  >====--" << endl;
        }
        out << index.mySyntax(i);  // << endl; // ends with newline
        break;
      }
  }


  void PrintSearchExact(std::ostream &out, const index::entry& ManEntry)
{
  const string title(ManEntry.myTitle());
  const string EntryFileName(ManEntry.myFileName());
  ifstream in(EntryFileName.c_str());
  if (!in)
  {
    std::cout << "Cannot find input file " << EntryFileName <<".  Aborting." << endl;
    abort();
  }
  string s = SkipToStringInTag(in, "<title>");
  while (!in.eof() && s != title)  s = SkipToStringInTag(in, "<title>");
  if (in.eof())
    CoCoA_THROW_ERROR("Not found exact title: " + title, "PrintSearchExact");
  out << "--<<<<<<<<<<<<<<( " << CleanupLine(title) << " )>>>>>>>>>>>>>>>--" << endl;
  PrintSyntax(out, title);
  string line;
  getline(in, line);
  while (!IsSubstring("<description>", line))  getline(in, line);
  //  out << "--====<  description  >====--" << endl;
  getline(in, line);
  while (!IsSubstring("</description>", line))
  {
    if ( IsSubstring("<obsolete_functions>", line))
      PrintAllMatchesFor(out, "[obsolete]");
    else
    {
    if ( IsSubstring("<obsolescent_functions>", line))
      PrintAllMatchesFor(out, "[obsolescent]");
    else
    {
    if ( IsSubstring("<commands_and_functions_for", line))
      PrintCommAndFunc(out, line);
    else
    {
      if ( IsSubstring("<commands_and_functions_rtn", line))
        PrintCommAndFuncRtn(out, line);
      else
      {
        if ( IsSubstring("<example>", line))  PrintExample(out, in);
        else
        {
          if ( IsSubstring("<!--", line))  SkipComment(in, line);
          else
            out << ProcessLine(line) << endl;
        }
      }
    }
    }
    }
    getline(in, line);
  }
  getline(in, line);
  while (!IsSubstring("<", line))   getline(in, line);
  PrintKeys(out, title);
  out << "--============( end " << CleanupLine(title) << " )=============--" << endl;
  PrintSeeAlso(out, title);
  //  out << "--============( end " << CleanupLine(title) << " )=============--" << endl;
}


#if 0
  // hopefully useless
  std::size_t SkipToLAngle(std::ifstream& in, std::string& line, std::size_t from)
  {
    size_t pos;
    while (!in.eof())
    {
      if ( (pos=line.find("<")) != string::npos) return pos;
      getline(in, line);
    }
    return string::npos;
  }
#endif


//----------------------------------------------------------------------

void SearchMatch(std::vector<long>& MatchingEntries, const std::string& s)
{
  const string ls(LowerCase(s));

  const index& index(UniqueIndex());
  for (size_t i=0; i<index.mySize(); ++i)
    for (size_t j=0; j<index.myKeys(i).size(); ++j)
      if (IsSubstring(ls, index.myKeys(i)[j]))
      {
        MatchingEntries.push_back(i);
        break;
      }
}


void SearchType(std::vector<long>& MatchingEntries, const std::string& s)
{
  const index& index(UniqueIndex());
  for (size_t i=0; i<index.mySize(); ++i)
    for (size_t j=0; j<index.myTypes(i).size(); ++j)
      if (s == index.myTypes(i)[j])
      {
        MatchingEntries.push_back(i);
        break;
      }
}


void SearchRtnType(std::vector<long>& MatchingEntries, const std::string& s)
{
  const index& index(UniqueIndex());
  for (size_t i=0; i<index.mySize(); ++i)
    for (size_t j=0; j<index.myRtnTypes(i).size(); ++j)
      if (s == index.myRtnTypes(i)[j])
      {
        MatchingEntries.push_back(i);
        break;
      }
}


  void PrintCommAndFunc(std::ostream &out, std::string& TypeLine)
  {
    ReplaceWith(TypeLine, "<commands_and_functions_for type=\"", "");
    ReplaceWith(TypeLine, "\"></commands_and_functions_for>", "");
    ifstream in = OpenXMLManual();
    std::vector<long> found;
    SearchType(found, TypeLine);
    const index& index(UniqueIndex());
    for (size_t i=0; i<found.size(); ++i)
    {
      in.seekg(0);///      ifstream in(CoCoAManFileName());  // if (!in) ...
      string s = SkipToStringInTag(in, "<title>");
      while (!in.eof() && s != index.myTitle(found[i]))
        s = SkipToStringInTag(in, "<title>");
      if (in.eof())
        CoCoA_THROW_ERROR("Not found exact title: " + index.myTitle(found[i]), "PrintCommAndFunc");
      out << "? " << s << " -- ";
      string line;
      getline(in, line);
      while (!IsSubstring("<short_description>", line))  getline(in, line);
      out << StringInTag(line, "<short_description>") << endl;
    }
    out << "--==================<>===================--" << endl;
  }


  void PrintCommAndFuncRtn(std::ostream &out, std::string& TypeLine)
  {
    ifstream in = OpenXMLManual();
    ReplaceWith(TypeLine, "<commands_and_functions_rtn type=\"", "");
    ReplaceWith(TypeLine, "\"></commands_and_functions_rtn>", "");
    std::vector<long> found;
    SearchRtnType(found, TypeLine);
    const index& index(UniqueIndex());
    for (size_t i=0; i<found.size(); ++i)
    {
      in.seekg(0); ///      ifstream in(CoCoAManFileName());  // if (!in) ...
      string s = SkipToStringInTag(in, "<title>");
      while (!in.eof() && s != index.myTitle(found[i]))
        s = SkipToStringInTag(in, "<title>");
      if (in.eof())
        CoCoA_THROW_ERROR("Not found exact title: " + index.myTitle(found[i]), "PrintCommAndFuncRtn");
      out << "? " << s << " -- ";
      string line;
      getline(in, line);
      while (!IsSubstring("<short_description>", line))  getline(in, line);
      out << StringInTag(line, "<short_description>") << endl;
    }
    out << "--==================<>===================--" << endl;
  }


namespace // anon for file local fn
{
  // >>> IMPORTANT <<<
  // Here double-quotes are regarded as being "whitespace"
  // See Redmine issue 1021.

  // This fn is used only in CleanupKeyword
  long SkipWhitespace(const string& str, long i)
  {
    const long n = str.size();
    while (i < n && (isspace(str[i]) || str[i] == '"'))
      ++i;
    return i;
  }
}

  string CleanupKeyword(const std::string& key)
  {
    using std::isspace;
    string keyword;
    const long n = key.size(); // BUG BUG BUG might overflow!!

    for (long i=SkipWhitespace(key,0); i < n; ++i)
    {
      bool InsertSpace = false;
      // Skip over any whitespace; set flag InsertSpace if there was whitespace.
      if (isspace(key[i]) || key[i] == '"')
      {
        i = SkipWhitespace(key,i);
        if (i == n) { break; }
        InsertSpace = true;
      }
      // Not at end of key, and curr char is not whitespace.
      const char ch = key[i];
      // Break out if looking at ";" or "//" or "--"
      if (ch == ';') break; // skip anything after ';'
      if (i != n-1) // skip comments
      {
        if (ch == '/' && key[i+1] == '/') break;
        if (ch == '-' && key[i+1] == '-') break;
      }
      if (InsertSpace) keyword += ' ';
      keyword += ch;
    }
    return keyword;
  }


void PrintMan(std::ostream &out, std::string keyword)
  {
    const index& index(UniqueIndex());
    //      std::cout << " --after UniqueIndex" << std::endl;
    keyword = CleanupKeyword(keyword);
    //      std::cout << " --after CleanupKeyword" << std::endl;
    enum { SingleQuery, DoubleQuery } HelpType = SingleQuery;

    // Check whether there is a 2nd '?'; if so, skip it & whitespace
    if (keyword[0]=='?')
    {
      HelpType = DoubleQuery;
      keyword = CleanupKeyword(keyword.substr(1)); // skip first char: '?'
    }

    // Empty manual query produces a short guide
    if (HelpType == SingleQuery && keyword.empty())
    {
      out << "--======( Quick guide to online help )======--" << endl;
      out << "  ? ring                     --> best match for `ring'" << endl;
      out << "  ?? ring                    --> all matches for `ring'" << endl;
      out << "  indent(starting(\"ring\"));  --> all functions *starting* with `ring' or `Ring'" << endl;
      out << endl;
      out << "--> Some helpful topics:" << endl;
      out << "  ? tutorial                 --> list of CoCoA-5 tutorials" << endl;
      out << "  ? syntax                   --> page about CoCoA-5 language commands" << endl;
      out << "  ? operators, shortcuts     --> page about CoCoA-5 language operators" << endl;
      out << "--======( End of quick guide )======--" << endl;
      return;
    }
    vector<long> AllMatches;
    SearchMatch(AllMatches, keyword);
    //    std::cout << " --after SearchMatch" << std::endl;
    if (AllMatches.empty())
    {
      out << "--====<  No matches for \"" << keyword << "\"  >====--" << endl;
      return; //-->>
    }
    const long NumAllMatches = AllMatches.size();
    if (HelpType == SingleQuery && NumAllMatches == 1)
    {
      PrintSearchExact(out, index.myEntry(AllMatches[0]));
      return; //-->>
    }
    long ExactMatch = -1;
    vector<long> ExactKeyMatches;
    if (HelpType == SingleQuery)
      for (long i=0; i<NumAllMatches; ++i)
      {
        long CurrentMatch = AllMatches[i];
        if (LowerCase(index.myTitle(CurrentMatch)) == LowerCase(keyword))
        {
          ExactMatch = CurrentMatch;
          PrintSearchExact(out, index.myEntry(CurrentMatch));
          break;
        }
        const std::vector<std::string>& keys_i = index.myKeys(CurrentMatch);
        for (long n=keys_i.size()-1; n>=0; --n)
          if (LowerCase(keys_i[n]) == LowerCase(keyword))
          {
            ExactKeyMatches.push_back(CurrentMatch);
            break;
          }
      }
    if (ExactMatch == -1)
    {
      if (ExactKeyMatches.size()==1)
      {
        ExactMatch = ExactKeyMatches[0];
        PrintSearchExact(out, index.myEntry(ExactMatch));
      }
      else
      {
        if (HelpType == SingleQuery)
          out << "--====<  No exact match for \""
              << keyword << "\"  >====--" << endl;
        PrintAllMatchesFor(out, keyword);
        return; //-->>
      }
    }
    if (NumAllMatches > 6)
    {
      out << "--> To see all " << NumAllMatches << " matches for \""
          << keyword << "\":" << endl;
      out << " ?? " << keyword << endl;
    }
    else // (NumAllMatches <= 6)
    {
      if (ExactMatch == -1)
      {
        out << "--> All further matches for \"" << keyword << "\":\n";
        for (int i=0; i<NumAllMatches; ++i)
          out << " ? " << index.myTitle(AllMatches[i]) << endl;
      }
      else // (ExactMatch != -1)
      {
        vector<long> FurtherMatches;
        const vector<string>& EMSeeAlso = index.mySeeAlso(ExactMatch);
        for (int i=0; i<NumAllMatches; ++i)
        {
          bool AlreadyCited = false;
          if (AllMatches[i] == ExactMatch) AlreadyCited = true;
          for (size_t s=0; s<EMSeeAlso.size(); ++s)
            if ( index.myTitle(AllMatches[i]) == EMSeeAlso[s] )
              AlreadyCited = true;
          if (!AlreadyCited) FurtherMatches.push_back(AllMatches[i]);
        }
        if (!FurtherMatches.empty())
        {
          if (FurtherMatches.size() == 1)
            out << "--> One further match for \"" << keyword << "\":\n";
          else
            out << "--> All " << FurtherMatches.size() << " further matches for \"" << keyword << "\":\n";
          for (size_t i=0; i<FurtherMatches.size(); ++i)
            out << " ? " << index.myTitle(FurtherMatches[i]) << endl;
        }
      }
    }
  }


void PrintAllMatchesFor(std::ostream &out,
                        const std::string& str,
                        const std::vector<long>& found)
  {
    if (found.empty())
    {
      out << "--> No matches for \"" << str << "\"" << endl;
      return;
    }
    const int NumFound = found.size();
    if (NumFound == 1)
      out << "--> The only match for \"" << str << "\":\n";
    else
      out << "--> All " << NumFound << " matches for \"" << str << "\":\n";
//    out << " (" << NumFound << ")" << endl;
    const index& index(UniqueIndex());
    for (int i=0; i<NumFound; ++i)
      out << " ? " << index.myTitle(found[i]) << endl;
  }


void PrintAllMatchesFor(std::ostream &out, const std::string& str)
  {
//     vector<string> found;
    vector<long> found;
    SearchMatch(found, str);
    PrintAllMatchesFor(out, str, found);
  }


  void ReloadMan(std::ostream &out)
  {
    UniqueIndex().myLoad(out);
  }


void ReloadMan(std::ostream &out, const std::vector<std::string>& inFileNames)
  {
    UniqueIndex().myLoad(out);
    const long NumFiles = inFileNames.size(); // silent cast from unsigned long
    for (long i=0; i < NumFiles; ++i)
      UniqueIndex().myLoadExtraXML(out, inFileNames[i]);
  }


void PrintKeys(std::ostream &out, const std::string& s)
  {
    const index& index(UniqueIndex());
    for (size_t i=0; i<index.mySize(); ++i)
      if ( LowerCase(index.myTitle(i)) == LowerCase(s) )
      {
        if (index.myKeys(i).size() != 1)
        {
          out << "-- SEARCH KEYS: " << index.myKeys(i)[1];
          for (size_t j=2; j<index.myKeys(i).size(); ++j)
            out << ", " << index.myKeys(i)[j];
          out << endl;
        }
        break;
      }
  }


void PrintSeeAlso(std::ostream &out, const std::string& s)
  {
    const index& index(UniqueIndex());
    for (size_t i=0; i<index.mySize(); ++i)
      if ( LowerCase(index.myTitle(i)) == LowerCase(s) )
      {
        if (index.mySeeAlso(i).size() != 0)
          out << "--> See also:" << endl;
        for (size_t j=0; j<index.mySeeAlso(i).size(); ++j)
          out << " ? " << index.mySeeAlso(i)[j] << endl;
        break;
      }
  }


void PrintAllExamples(std::ostream &out)
{
  ifstream in = OpenXMLManual();
  // ifstream in(CoCoAManFileName());
  // if (!in)
  // {
  //   std::cout << "Cannot find input file `CoCoAHelp.xml`.  Aborting." << endl;
  //   abort();
  // }
  string s; /// = "";
  while (!in.eof())
  {
    s = SkipToStringInTag(in, "<title>");
    if (in.eof()) break;
    out << "PrintLn \"-- " << s << " =============================\";" << endl;
    string line;
    getline(in, line);
    while (!IsSubstring("<description>", line))  getline(in, line);
    getline(in, line);
    while (!IsSubstring("</description>", line))
    {
      if ( IsSubstring("<example>", line) && !IsSubstring("<!--", line) )
        PrintExample(out, in);
      getline(in, line);
    }
  }
  out << "--- PrintAllExamples: done ---" << endl;
}


// ManExamples //////////////////////
void PrintAllExamplesWithoutOutput(std::ostream &out)
{
  ifstream in = OpenXMLManual();
  string s; /// = "";
  while (!in.eof())
  {
    s = SkipToStringInTag(in, "<title>");
    if (in.eof()) break;
    out << "PrintLn \"-- " << s << " =============================\";" << endl;
    out << "SetVerbosityLevel(0);" << endl;
    string line;
    getline(in, line);
    while (!IsSubstring("<description>", line))  getline(in, line);
    getline(in, line);
    while (!IsSubstring("</description>", line))
    {
      if ( IsSubstring("<example>", line) && !IsSubstring("<!--", line) )
        PrintExampleWithoutOutput(out, in);
      getline(in, line);
    }
  }
  out << "--- PrintAllExamples: done ---" << endl;
}


void PrintWordlist(std::ostream &out)
{
  ifstream in = OpenXMLManual();
  // ifstream in(CoCoAManFileName());
  // if (!in)
  // {
  //   std::cout << "Cannot find input file `CoCoAHelp.xml`.  Aborting." << endl;
  //   abort();
  // }
///???  string s = "";
  string line;
  getline(in, line);
  while (!in.eof())
  {
    while (!IsSubstring("<command>", line) && !in.eof())  getline(in, line);
    out << SkipToStringInTag(in, "<title>") << endl;
    getline(in, line);
  }
}


std::vector<std::string> WordList()
{
  ifstream in = OpenXMLManual(); 
  string line;
  std::vector<string> WL;
  getline(in, line);
  while (!in.eof())
  {
    while (!IsSubstring("<command>", line) && !in.eof())  getline(in, line);
    WL.push_back(SkipToStringInTag(in, "<title>"));
    getline(in, line);
  }
  return WL;
}


const std::string& CoCoAManFileName()
{
  static const string UniqueCopy(packageDir+"/../CoCoAManual/CoCoAHelp.xml");
  return UniqueCopy;
}


//const index& UniqueIndex()
index& UniqueIndex()
{
  static index UniqueCopy;
  return UniqueCopy;
}

} // namespace OnlineHelp
} // namespace CoCoA
