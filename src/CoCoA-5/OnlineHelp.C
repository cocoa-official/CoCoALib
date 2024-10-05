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



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/CoCoA-5/OnlineHelp.C,v 1.104 2024/03/08 20:00:38 abbott Exp $
// $Log: OnlineHelp.C,v $
// Revision 1.104  2024/03/08 20:00:38  abbott
// Summary: Removed unnecessary include
//
// Revision 1.103  2023/12/01 15:34:46  bigatti
// Summary: Added WordList() returning vector of string
//
// Revision 1.102  2023/06/25 15:20:32  bigatti
// Summary: resetting VerbosityLevel to 0 before each example
//
// Revision 1.101  2023/04/06 14:11:46  bigatti
// Summary: new look for header/closing
//
// Revision 1.100  2023/03/27 14:02:45  abbott
// Summary: Added better mesg when single query with no string is input
//
// Revision 1.99  2022/03/07 16:54:46  bigatti
// Summary: shortened line "SillyName" because readline truncates long echo-lines
//
// Revision 1.98  2022/03/03 14:25:28  bigatti
// Summary: added types
//
// Revision 1.97  2022/03/02 11:11:08  bigatti
// Summary: added: i b eg ie
//
// Revision 1.96  2022/02/22 20:39:27  abbott
// Summary: Updated copyright message (redmine 1555)
//
// Revision 1.95  2021/11/10 19:02:35  abbott
// Summary: Improved error mesg (redmine 1397)
//
// Revision 1.94  2020/10/23 08:28:11  bigatti
// Summary: renamed XMLFileName into CoCoAManFileName
// and other "FileName"s in OnlineHelp.C  for readability
//
// Revision 1.93  2020/10/14 13:10:46  abbott
// Summary: Changed XMLFileName so that result is a std::string
//
// Revision 1.92  2020/08/10 10:14:49  bigatti
// Summary: print also search keys
//
// Revision 1.91  2020/06/17 15:49:30  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.90  2019/10/10 16:38:07  bigatti
// -- text for clearing variables in only one line
//
// Revision 1.89  2019/10/09 12:57:15  abbott
// Summary: Changed default ring to have coeffs in ZZ/(4)
//
// Revision 1.88  2018/05/17 15:58:02  bigatti
// -- renamed VectorOperations --> VectorOps
//
// Revision 1.87  2018/03/13 14:57:18  abbott
// Summary: Minor cleaning
//
// Revision 1.86  2018/03/13 14:53:32  abbott
// Summary: Added c_str to args of ctor for ifstream
//
// Revision 1.85  2018/03/08 17:00:07  bigatti
// -- ReloadMan: now many manual files can be loaded
//
// Revision 1.84  2017/11/15 16:33:49  abbott
// Summary: Manual search now treats double-quotes as spaces
//
// Revision 1.83  2017/10/09 22:27:40  abbott
// Summary: Now resets vars R and P between tests; also sets sily curr ring
//
// Revision 1.82  2017/09/22 13:29:48  abbott
// Summary: Added two spaces for better formatted output
//
// Revision 1.81  2017/09/06 11:57:03  abbott
// Summary: Changed ERR::SERIOUS into ERR:ShouldNeverGetHere
//
// Revision 1.80  2017/06/29 10:36:06  abbott
// Summary: Major revision to CleanupKeyword (see redmine 1021)
//
// Revision 1.79  2017/04/10 12:40:52  bigatti
// -- using () instead of [] for titles
// -- added end line
//
// Revision 1.78  2016/11/21 17:58:27  abbott
// Summary: Minor change: print "the only" when there is just 1 hit
//
// Revision 1.77  2016/10/10 14:39:35  bigatti
// -- commented out references to version
//
// Revision 1.76  2016/06/22 14:22:49  abbott
// Summary: Corrected condition for checking if ifstream was created OK
//
// Revision 1.75  2016/02/18 08:10:16  bigatti
// -- removed printout on clog -- nasty trick: printing on ostringstream
//
// Revision 1.74  2016/01/29 13:35:15  bigatti
// -- removed some newlines, added some "-->" to highlight text
//    (overall the output looks more compact and readable: too many times
//    I thought I had no result from the manual because it scrolled too much)
//
// Revision 1.73  2015/07/27 14:40:55  abbott
// Summary: Special case message if there is exactly 1 "further match"
//
// Revision 1.72  2014/08/01 11:23:45  bigatti
// -- added skipping comments for ManExamples
//
// Revision 1.71  2014/07/18 12:48:07  bigatti
// -- improvement for make ManExamples
//
// Revision 1.70  2014/07/18 11:21:12  bigatti
// -- minor cleanup
//
// Revision 1.69  2014/07/18 11:08:35  bigatti
// -- removed unwanted debugging print
// -- fixed trailing spaces before comment
//
// Revision 1.68  2014/07/18 08:27:54  bigatti
// -- now "all further matches" does not print what was in seealso
//
// Revision 1.67  2014/07/16 13:24:57  bigatti
// -- added obsolete_functions and obsolescent_functions
//
// Revision 1.66  2014/05/14 07:29:37  bigatti
// -- improved (fixed) PrintExampleWithoutOutput
//
// Revision 1.65  2014/05/09 10:57:33  bigatti
// -- added -- before printing title so that it gets coloured as comments
//
// Revision 1.64  2014/04/24 11:37:34  abbott
// Summary: Added removal of <rtn> and </rtn> to CleanupLine
// Author: JAA
//
// Revision 1.63  2014/04/24 11:03:28  bigatti
// -- rtntype cahnged into rtn
//
// Revision 1.62  2014/04/23 16:48:44  bigatti
// -- fixed RtnTypes
//
// Revision 1.61  2014/04/23 14:15:38  bigatti
// ++ added myAddRetType
//
// Revision 1.60  2014/04/22 17:18:21  bigatti
// ++ removed backslash (no longer needed)
//
// Revision 1.59  2014/04/22 12:46:07  abbott
// Summary: Text inside <em>...</em> now printed between *...*
// Author: JAA
//
// Revision 1.58  2014/04/09 13:28:38  bigatti
// -- removed printing title for syntax and description
//
// Revision 1.57  2014/04/07 11:10:34  bigatti
// -- fixed matches for case ExactMatchKey==1
//
// Revision 1.56  2014/04/07 10:25:29  abbott
// Summary: Added {} as suggested by compiler; added condition to avoid duplicate in "other matches"
// Author: JAA
//
// Revision 1.55  2014/04/04 13:42:00  bigatti
// -- another fix
//
// Revision 1.54  2014/04/04 12:06:24  bigatti
// -- another fix
//
// Revision 1.53  2014/04/04 11:13:28  bigatti
// -- another fix
//
// Revision 1.52  2014/04/04 11:07:44  bigatti
// -- fix
//
// Revision 1.51  2014/04/04 11:01:04  bigatti
// -- now also searching if there is a single exact match in keys
//
// Revision 1.50  2014/03/28 15:23:41  bigatti
// -- removed "syntax" for "commands and functions for .."
//
// Revision 1.49  2014/03/27 17:33:58  bigatti
// -- added "No entry for "keyword"" message
//
// Revision 1.48  2014/03/26 16:28:04  abbott
// Summary: Improved messages printed out by online help
// Author: JAA
//
// Revision 1.47  2014/03/26 15:20:26  bigatti
// -- not printing "all matches" automatically if there are too many
//
// Revision 1.46  2014/03/26 12:02:01  abbott
// Summary: Major change to PrintMan: if search string starts with '?' then just list titles of matching manual entries
// Author: JAA
//
// Revision 1.45  2014/03/19 15:51:54  abbott
// Summary: Added blank line after syntax section; changed "See Also" --> "see also"
// Author: JAA
//
// Revision 1.44  2014/03/06 15:52:52  abbott
// Summary: Cleaned impl of PrintMan
// Author: JAA
//
// Revision 1.43  2014/03/04 14:27:54  bigatti
// -- added CleanupLine for myAddKeys
//
// Revision 1.42  2013/06/12 08:55:21  bigatti
// -- fixed "syntax" for error in entry "CartesianProduct"
//
// Revision 1.41  2013/02/27 10:44:24  bigatti
// -- added <tt> in CleanupLine (actually useless)
//
// Revision 1.40  2013/02/26 14:18:01  bigatti
// -- added field for syntax in "index"
//
// Revision 1.39  2013/02/26 13:40:31  bigatti
// -- added CleanupLine to "title"  (for the variable "it")
//
// Revision 1.38  2013/02/22 18:09:25  bigatti
// -- removed unused arg out from SkipComment
//
// Revision 1.37  2013/02/22 11:03:41  bigatti
// -- improved SkipComment: more robust
//
// Revision 1.36  2013/02/22 10:34:28  bigatti
// ++ added SkipComment (for multiline comments ONLY!)
//
// Revision 1.35  2013/02/19 18:52:50  abbott
// Added line to ignore HTML italic indications <i>...</i>.
//
// Revision 1.34  2012/06/19 15:28:47  bigatti
// -- aesthetics: changed delimiters for <example>
//
// Revision 1.33  2012/06/18 10:02:44  bigatti
// -- added mechanism to get version number and date from the xml file
//
// Revision 1.32  2012/06/04 09:34:05  bigatti
// -- added PrintWordlist
//
// Revision 1.31  2012/04/04 13:56:35  bigatti
// -- added PrintAllExamplesWithoutOutput
//
// Revision 1.30  2012/04/02 15:15:37  bigatti
// -- added check for trailing tab
//
// Revision 1.29  2012/03/20 14:34:09  bigatti
// -- added recognition of some xml tags
// -- minor cleaning
// -- Printing "matches for .." only when meaningful
//
// Revision 1.28  2012/03/13 15:59:12  bigatti
// -- xml source moved into CoCoAManual/
//
// Revision 1.27  2012/03/09 17:01:27  bigatti
// -- removed warning about old manual
// -- added "end example"
//
// Revision 1.26  2012/03/07 18:10:40  bigatti
// -- fixed interpretation of <backslash/>
//
// Revision 1.25  2012/02/24 13:10:12  bigatti
// -- added ReloadMan
//
// Revision 1.24  2012/01/25 14:33:15  bigatti
// -- added code for checking proper sorting (inactive)
//
// Revision 1.23  2011/11/02 14:41:16  bigatti
// -- added syntax
// -- removing trailing comments and spaces to searching string
// -- minor aesthetics change
//
// Revision 1.22  2011/07/06 15:52:13  bigatti
// -- improved IsSubstring
//
// Revision 1.21  2011/05/25 12:19:15  bigatti
// -- added case <par />
//
// Revision 1.20  2011/05/24 08:17:16  bigatti
// -- added warning also at the end of description
//
// Revision 1.19  2011/05/23 13:30:37  bigatti
// -- chencged include CoCoA/library --> CoCoA/error
//
// Revision 1.18  2011/05/13 10:26:37  bigatti
// -- added warning "thi is CoCoA-4 manual"
//
// Revision 1.17  2011/05/04 12:03:50  lagorio
// *** empty log message ***
//
// Revision 1.16  2011/03/23 17:31:22  bigatti
// -- added <code>..</code>
//
// Revision 1.15  2011/02/16 17:22:40  bigatti
// -- added "See Also"
//
// Revision 1.14  2011/02/16 16:14:09  bigatti
// -- class deefinition moved into .C file
// -- added storing of types, and function for <commands_and_functions_for>
// -- cleaning up
//
// Revision 1.13  2011/02/14 10:10:08  bigatti
// -- fixed &amp; &apos; <backslash/>  ...
// -- added function PrintAllExamples
//
// Revision 1.12  2011/01/25 15:45:23  bigatti
// -- added tags verbatim and formula
//
// Revision 1.11  2010/12/26 13:14:44  abbott
// Changed GlobalLogput to clog.
//
// Revision 1.10  2010/11/22 17:45:02  bigatti
// -- fixed bug in storing keys
//
// Revision 1.9  2010/11/22 17:37:44  bigatti
// -- improved message for entries with no matching
//
// Revision 1.8  2010/09/02 13:04:08  bigatti
// -- added conversion for <em>
//
// Revision 1.7  2010/09/01 13:24:48  bigatti
// -- moved all manual functions into CoCoA::OnlineHelp namespace
//
// Revision 1.6  2010/09/01 12:27:52  lagorio
// *** empty log message ***
//
// Revision 1.5  2010/09/01 09:23:19  bigatti
// -- added some "const"
// -- using CoCoA::GlobalLogput and CoCoA_ERROR instead of clog
//
// Revision 1.4  2010/09/01 08:28:35  lagorio
// *** empty log message ***
//
// Revision 1.3  2010/09/01 07:46:19  lagorio
// *** empty log message ***
//
// Revision 1.2  2010/08/31 14:55:02  bigatti
// -- pretty basic fixes
//

