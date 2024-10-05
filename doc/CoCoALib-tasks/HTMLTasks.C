//   Copyright (c)  2011 Anna Bigatti

//   This file is part of the source of CoCoALib, the CoCoA Library.

//   CoCoALib is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.

//   CoCoALib is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.

//   You should have received a copy of the GNU General Public License
//   along with CoCoALib.  If not, see <http://www.gnu.org/licenses/>.

//#include "CoCoA/time.H"
//using CoCoA::DateTime;  -- hard copied here
#include <algorithm>
using std::transform;
#include <cctype>
//using std::tolower;
#include <ctime>
using std::time;
#include <fstream>
using std::ifstream;
using std::ofstream;
#include <iostream>
using std::endl;
#include <string>
using std::string;
#include <vector>
using std::vector;
#include <string.h>

extern string packageDir;

namespace HTMLTasks
{

  //--------------------------------------------------
  //--------------------------------------------------  

  class index
  {
  public: 
    class entry
    {
    public:
      entry(const std::string& title);
      // ~entry(); // default is OK
      const std::string& myTitle() const {return myTitleValue;}
      //      const std::vector<std::string>& myAuthors() const {return myAuthorsValue;}
      const std::string& myAuthors() const {return myAuthorsValue;}
      const std::string& myDescription() const {return myDescrValue;}
      const std::string& myStatus() const {return myStatusValue;}
      const std::string& myStage() const {return myStageValue;}
      bool IamHighPriority() const {return IamHighPriorityValue;}
      const std::string& myComplDate() const {return myComplDateValue;}
      const std::string& myLastNewsDate() const {return myLastNewsDateValue;}
      void myAddDescription(std::ifstream& in);
      //      void myAddAuthors(std::ifstream& in);
      //      void myAddAuthor(const std::string& t);
      void myAddAuthors(const std::string& t);
      void myAddStatus(const std::string& s);
      void myAddStage(const std::string& s);
      void myAddHighPriority(const std::string& s);
      void myAddComplDate(const std::string& s);
      void myAddLastNewsDate(const std::string& s);
      void myReadEntry(std::ifstream& in, const std::string& ClosingTag);

    private: // data members
      std::string myTitleValue;
      std::string myDescrValue;
      //      std::vector<std::string> myAuthorsValue;
      std::string myAuthorsValue;
      std::string myStatusValue;
      std::string myStageValue;
      bool IamHighPriorityValue;
      std::string myComplDateValue;
      std::string myLastNewsDateValue;
    };
 
  public:
    index();
    // ~index(); // default is OK
    //    void myLoad();
    const index::entry& myEntry(std::size_t i) const {return myVec[i];}
    const std::string& myTitle(std::size_t i) const {return myVec[i].myTitle();}
    const std::string& myStatus(std::size_t i) const {return myVec[i].myStatus();}
    const std::string& myStage(std::size_t i) const {return myVec[i].myStage();}
    bool IamHighPriority(std::size_t i) const {return myVec[i].IamHighPriority();}
    const std::string& myComplDate(std::size_t i) const {return myVec[i].myComplDate();}
    const std::string& myLastNewsDate(std::size_t i) const {return myVec[i].myLastNewsDate();}
    const std::string& myDescription(std::size_t i) const {return myVec[i].myDescription();}
    //    const std::vector<std::string> myAuthors(std::size_t i) const {return myVec[i].myAuthors();}
    const std::string& myAuthors(std::size_t i) const {return myVec[i].myAuthors();}
    std::size_t mySize() const {return myVec.size();}
    void mySortAlphabetically();

  private:
    void myAddEntry(std::ifstream& in, const std::string& title, const std::string& ClosingTag);
    void myAddCompletedEntry(std::ifstream& in, const std::string& title, const std::string& ClosingTag);
  
  private: // data members
    std::vector<entry> myVec;
  };

  // ----------------------------------------------------------------------
  // fwd decl -- defined later in this file
  bool CmpTitles(const index::entry& t1, const index::entry& t2);
  std::string WithoutSpaces(const std::string& title);
  std::string LinkName(const std::string& title);
  std::string ClosedTag(const std::string& tag);
  void ReplaceWith(std::string& line, const std::string& SFrom, const std::string& STo);
  bool IsDateInFuture(const std::string& s);
  bool IsMoreThan6MonthsAgo(const std::string& s);
  bool IsMoreThan12MonthsAgo(const std::string& s);
  long YearInDate(const std::string& s);
  long YearToday();
  const char* XMLFileName();
  std::string LowerCase(const std::string& str);
  bool IsStr1BeforeStr2(std::ifstream& in, const std::string& s1, const std::string& s2);
  std::string StringInTag(const std::string& line, const std::string& XMLTag);
  //  std::string SkipToStringInTag(std::ifstream& in, const std::string& XMLTag);
  std::string SkipToStringInTag(std::ifstream& in, std::string& line, const std::string& XMLTag);
  void SkipToTagLine(std::ifstream& in, std::string& line, const std::string& XMLTag);
  //  void PrintCommAndFunc(std::string& line, std::ostream &out);
  void PrintEntry(const std::string& s, std::ostream &out);
  void PrintWarning(const std::string& where, const std::string& what);
  void PrintDate(std::ostream &out, const std::string& s, const std::string& tag);
  bool IsCommentLine(const std::string& line);
  bool IsEmptyLine(const std::string& line);

  // ----------------------------------------------------------------------
  inline bool IsSubstring(const std::string& s, const std::string& line) 
  { return line.find(s) != string::npos; }

#if defined(__MINGW32__) || defined(_WINDOWS)
  void DateTime(long& date, long& time)
  {
    std::time_t SecondsSinceEpoch;
    std::time(&SecondsSinceEpoch); // gets current time - Async-signal-safe
//   const string date = std::ctime_r(&SecondsSinceEpoch);
//   cout << "Time & date: " << date << endl;

    std::tm* CurrentTime;
    CurrentTime = localtime(&SecondsSinceEpoch); // thread-safe in windows
    const int year = CurrentTime->tm_year + 1900;
    const int month = CurrentTime->tm_mon + 1;
    const int day = CurrentTime->tm_mday;
    const int hour = CurrentTime->tm_hour;
    const int minute = CurrentTime->tm_min;
    const int second = CurrentTime->tm_sec;

    date = 10000*year + 100*month + day;
    time = 10000*hour + 100*minute + second;
  }
  
#else

  void DateTime(long& date, long& time)
  {
    std::time_t SecondsSinceEpoch;
    std::time(&SecondsSinceEpoch); // gets current time - Async-signal-safe
//   const string date = std::ctime_r(&SecondsSinceEpoch);
//   cout << "Time & date: " << date << endl;

    std::tm CurrentTime;
    localtime_r(&SecondsSinceEpoch, &CurrentTime); // thread-safe
    const int year = CurrentTime.tm_year + 1900;
    const int month = CurrentTime.tm_mon + 1;
    const int day = CurrentTime.tm_mday;
    const int hour = CurrentTime.tm_hour;
    const int minute = CurrentTime.tm_min;
    const int second = CurrentTime.tm_sec;

    date = 10000*year + 100*month + day;
    time = 10000*hour + 100*minute + second;
  }

#endif
  // ----------------------------------------------------------------------

  //-- index --

  const index& UniqueIndex();

  index::index()
  {
    ifstream in(XMLFileName());
    if (!in)
      throw "Cannot find input file `tasks.xml`";
    index::entry LastEntry("AAAA");
    string line="";
    //    std::clog << "Loading task table" << std::endl;
    SkipToTagLine(in, line, "<tasks>");
    std::cout << "Loading open tasks\n";
    while (!in.eof())
    {
      getline(in, line);
      //      if (IsSubstring("<!--", line)) continue;
      if (IsCommentLine(line)) continue;
      if (IsEmptyLine(line)) continue;
      if (IsSubstring("</tasks>", line)) break;
      std::cout << ".";
      if (!myVec.empty()) LastEntry = myVec.back();
      if (IsSubstring("<task>", line))
        myAddEntry(in, SkipToStringInTag(in, line, "<title>"), "</task>");
      else
        PrintWarning(line, "unrecognized line");
      if (CmpTitles(myVec.back(), LastEntry))
        PrintWarning(myVec.back().myTitle(), "wrong sorting");
    }
    SkipToTagLine(in, line, "<completed-tasks>");
    std::cout << "\nLoading completed tasks\n";
    while (!in.eof())
    {
      getline(in, line);
      //      if (IsSubstring("<!--", line)) continue;
      if (IsCommentLine(line)) continue;
      if (IsEmptyLine(line)) continue;
      if (IsSubstring("</completed-tasks>", line)) break;
      std::cout << ".";      
      if (IsSubstring("<task>", line))
        myAddCompletedEntry(in, SkipToStringInTag(in, line, "<title>"), "</task>");
      else
        PrintWarning(line, "unrecognized line");
    }
    std::cout << endl;
    //    std::clog << "done" << endl;
  }


  void index::myAddEntry(std::ifstream& in, const std::string& title, const std::string& ClosingTag)
  {
    myVec.push_back(entry(title));
    entry& e(myVec.back());
    e.myReadEntry(in, ClosingTag);
    // checks
    if (e.myStatus() == "completed")
      PrintWarning(e.myTitle(), "move task to <completed-tasks>");
    if (!IsDateInFuture(e.myComplDate()))
      PrintWarning(e.myTitle(), "UPDATE COMPLETION DATE");
  }


  void index::myAddCompletedEntry(std::ifstream& in, const std::string& title, const std::string& ClosingTag)
  {
    myVec.push_back(entry(title));
    entry& e(myVec.back());
    e.myAddStatus("<status>completed</status>");
    e.myReadEntry(in, ClosingTag);
    // checks
    if (e.myStatus() != "completed")
      PrintWarning(e.myTitle(), "non-completed task in completed table");
  }


  bool CmpTitles(const index::entry& t1, const index::entry& t2)
  {
    return t1.myTitle().compare(t2.myTitle()) < 0;
  }
  

  void index::mySortAlphabetically()
  {
    sort(myVec.begin(), myVec.end(), CmpTitles);
  }
  

  //-- index::entry --
  index::entry::entry(const std::string& title):
    myTitleValue(title)
  {
    IamHighPriorityValue = false;
  }

  void index::entry::myAddComplDate(const std::string& line)
  {
    myComplDateValue = WithoutSpaces(StringInTag(line, "<compl-date>"));
  }

  void index::entry::myAddLastNewsDate(const std::string& line)
  {
    myLastNewsDateValue = WithoutSpaces(StringInTag(line, "<lastnews-date>"));
  }

  void index::entry::myAddStatus(const std::string& line)
  { myStatusValue = LowerCase(StringInTag(line, "<status>")); }

  void index::entry::myAddStage(const std::string& line)
  { myStageValue = LowerCase(StringInTag(line, "<stage>")); }

  void index::entry::myAddHighPriority(const std::string& line)
  { IamHighPriorityValue = IsSubstring("true", LowerCase(StringInTag(line, "<high-priority>"))); }

//   void index::entry::myAddAuthor(const std::string& t)
//   { myAuthorsValue.push_back(t); }

//   void index::entry::myAddAuthors(std::ifstream& in)
//   {
//     string line;
//     getline(in, line);
//     string s; // used only inside while loop
//     while ((!in.eof()) && (!IsSubstring("</authors>",line)))
//     {
//       if ((s=StringInTag(line, "<author>")) != "") myAddAuthor(s);
//       getline(in, line);
//     }
//   }

  void index::entry::myAddAuthors(const std::string& line)
  { myAuthorsValue = StringInTag(line, "<authors>"); }


  void index::entry::myAddDescription(std::ifstream& in)
  {
    string line;
    getline(in, line);
    while ((!in.eof()) && (!IsSubstring("</description>",line)))
    {
      //      ProcessLine(line);
      myDescrValue += line;
      myDescrValue += "\n";
      getline(in, line);
    }
  }

  void index::entry::myReadEntry(std::ifstream& in, const std::string& ClosingTag)
  {
    string line;
    getline(in, line);
    while (!in.eof() && !IsSubstring(ClosingTag, line))
    {
      if (IsSubstring("<status>", line)) myAddStatus(line);
      if (IsSubstring("<stage>", line)) myAddStage(line);
      if (IsSubstring("<high-priority>", line)) myAddHighPriority(line);
      if (IsSubstring("<compl-date>", line)) myAddComplDate(line);
      if (IsSubstring("<lastnews-date>", line)) myAddLastNewsDate(line);
      //      if (IsSubstring("<authors>", line)) myAddAuthors(in);
      if (IsSubstring("<authors>", line)) myAddAuthors(line);
      if (IsSubstring("<description>", line)) myAddDescription(in);
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
    std::cerr << "IsStr1BeforeStr2";
    throw "IsStr1BeforeStr2";
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
    {
      std::cerr << XMLTag+" closing tag not found in this line";
      throw " closing tag not found in this line";
    }
    size_t StartPos = open + XMLTag.length();
    return line.substr(StartPos, close-StartPos);
  }


  void SkipToTagLine(std::ifstream& in, std::string& line, const std::string& XMLTag)
  {
    while (!in.eof())
    {
      if (IsSubstring(XMLTag, line)) return;
      getline(in, line);
    }
  }
  

  std::string SkipToStringInTag(std::ifstream& in, std::string& line, const std::string& XMLTag)
  {
    //    string line;
    //    getline(in, line);
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
  

  std::string WithoutSpaces(const std::string& title)
  {
    std::string s = title;
    ReplaceWith(s, " ", "");
    return s;
  }

  std::string LinkName(const std::string& title)
  {
    std::string s = WithoutSpaces(title);
    ReplaceWith(s, "(", ""); ReplaceWith(s, ")", "");
    ReplaceWith(s, ",", "");
    ReplaceWith(s, ".", "");
    ReplaceWith(s, "<code>", "");   ReplaceWith(s, "</code>", "");
    ReplaceWith(s, "<ttref>", "");  ReplaceWith(s, "</ttref>", "");
    ReplaceWith(s, "<ref>", "");    ReplaceWith(s, "</ref>", "");
    ReplaceWith(s, "<em>", "");     ReplaceWith(s, "</em>", "");
    return s;
//     ReplaceWith(s, "<verbatim>", ""); ReplaceWith(s, "</verbatim>", "");
//     ReplaceWith(s, "<formula>", " "); ReplaceWith(s, "</formula>", " ");
//     ReplaceWith(s, "&lt;", "<");
//     ReplaceWith(s, "&gt;", ">");
//     ReplaceWith(s, "&apos;", "'");
//     ReplaceWith(s, "&amp;", "&");
//     ReplaceWith(s, "<backslash/>;", "\\");
//     if (IsSubstring("<par/>", s)) s = "";
  }


  bool IsCommentLine(const std::string& line)
  {
    if (IsSubstring("<!--", line))
    {
      if ( line.substr(0,4) == "<!--" && line.substr(line.size()-3,3) == "-->" )
        return true;
      PrintWarning(line, "only full line comments are allowed");
    }
    return false;
  }
  

  bool IsEmptyLine(const std::string& line)
  {
    if ( line.size() == 0 ) return true;
    return false;
  }
  

  bool IsDateInFuture(const std::string& s)  
  {
    long date, time_notused;
    DateTime(date, time_notused);
    const int tyear = date/10000;
    const int tmonth = (date/100)%100;
    const int tday = date%100;
    
    const int year  = atoi((s.substr(0, 4)).c_str());
    if (year > tyear || year == 0) return true;
    if (year < tyear) return false;
    const int month = atoi((s.substr(5, 2)).c_str());
    if (month > tmonth || month == 0) return true;
    if (month < tmonth) return false;
    if (s.size()<10) return true;
    const int day  = atoi((s.substr(8, 2)).c_str());
    if (day > tday) return true;
    //    if (tday > day) 
    return false;
  }


  bool IsMoreThan6MonthsAgo(const std::string& s)  
  {
    long date, time_notused;
    DateTime(date, time_notused);
    int tyear = date/10000;
    int tmonth = (date/100)%100;
    const int year  = atoi((s.substr(0, 4)).c_str());
    const int month = atoi((s.substr(5, 2)).c_str());
    if (tmonth > 6)  tmonth -= 6;  else  { tmonth +=6; --tyear; }
    if (year < tyear) return true;
    if (year > tyear) return false;
    if (month < tmonth) return true;
    return false;
  }


  bool IsMoreThan12MonthsAgo(const std::string& s)  
  {
    long date, time_notused;
    DateTime(date, time_notused);
    const int lastyear = date/10000 - 1;
    const int tmonth = (date/100)%100;
    const int year  = atoi((s.substr(0, 4)).c_str());
    const int month = atoi((s.substr(5, 2)).c_str());
    if (year < lastyear) return true;
    if (year > lastyear) return false;
    if (month < tmonth) return true;
    return false;
  }


  long YearInDate(const std::string& s)  
  {
    const string y = s.substr(0, 4);
    return atoi(y.c_str());
  }


  long YearToday()  
  {
    long date, time_notused;
    DateTime(date, time_notused);
    return date/10000;
  }


  void PrintWarning(const std::string& where, const std::string& what)
  {
    std::cerr << "\n**** " << where << " ******* " << what << " *******\n";
  }
  

void PrintYearLinks(std::ostream &out)
{
  long y = YearToday();
  out <<
    "<br><br>\n"
    "<!-- ============================================================ -->\n"
    "<center><a name=\"top\"><table cellpadding=\"5\"><tbody><tr>\n";
  for (long a=2006; a<=y ; ++a)
    out << "<th><a href=\"#" << a << "\">" << a << "</a></th>\n";
  out << "</tr></tbody></table></a></center>\n"
    "<!-- ============================================================ -->\n\n";
}


  void PrintDate(std::ostream &out, const std::string& s, const std::string& tag)
  {
    out << tag;
    if (s=="") { out << ClosedTag(tag); return; }
    if (s.size()<10)
      if (s.size()!=7)
      {
        std::cerr << "************************* "+s+" incomplete date\n";
      throw "incomplete date";
      }
    string y = s.substr(0, 4);
    if (y=="0000") { out << "????" << ClosedTag(tag); return; }
    out << y;
    string m = s.substr(5, 2);
    if (m=="00") { out << "-??" << ClosedTag(tag); return; }
    out << "-" << m;
    string d;
    if (s.size()!=7)
    {
      string d = s.substr(8, 2);
      if (d!="00") out << "-" << d;
    }
    out << ClosedTag(tag);
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




  void PrintEntry(const std::string& s, std::ostream &out)
  {
    const index& index(UniqueIndex());
    for (size_t i=0; i<index.mySize(); ++i)
      if ( LowerCase(index.myTitle(i)) == LowerCase(s) )
      {
        //        out << "----------------------------" << endl;
        out << "Title: " << index.myTitle(i) << endl;
//         out << "Status: " << index.myStatus(i) << endl;
//         out << "Stage: " << index.myStage(i) << endl;
//         out << "high priority: " << index.IamHighPriority(i) << endl;
//         out << "Date last news: ";
//         PrintDate(out, index.myLastNewsDate(i));
//         out << endl << "Date completion: ";
//         PrintDate(out, index.myComplDate(i));
//         out << endl;
//         out << "Authors: ";
//         for (size_t j=0; j<index.myAuthors(i).size(); ++j)
//           out << " - " << index.myAuthors(i)[j];
//         out << endl;
//         out << "Description: " << endl << index.myDescription(i) << endl;
        // checks
        if (index.myStatus(i)!="completed" &&
            !IsDateInFuture(index.myComplDate(i)))
          out << "******* ERROR ******* should be completed *******" << endl;
        break;
      }
  }
  



//----------------------------------------------------------------------
//  HTML blocks
//----------------------------------------------------------------------

void Print_CSS(std::ostream &out)
{
  out <<
"<style type=\"text/css\">\n"
"     body { background: #3cc; }\n"
"     b { background: #ff6; }\n"
"     table, a:hover, a.button:hover { background: #188; }\n"
"     table, div {margin-bottom: 6pt; }\n"
"     th { text-align:left; background: #dff; }\n"
"     tr { background: #fff; }\n"
"     a:hover { color: #fff; }\n"
"     p { border: ridge #cdd 1pt;}\n"
"     p:hover { background: #dfd;}\n"
"     h1, h2 { color: red; }\n"
"     h1 { text-align: center; }\n"
"     div, h1, h2\n"
"     {\n"
"       background: white;\n"
"       border: ridge #3cc  2pt;\n"
"       padding: 4pt;\n"
"     }\n"
"\n"
"  .completed { color: #060; font-style: italic; background: #afa; }\n"
"  .active { color: #600; font-weight: bold; background: #faa; }\n"
"  .class { color: #070; font-weight: bold }\n"
"  .file { color: red; }\n"
"  .func { color: blue; }\n"
"  .version { color: red; font-weight: bold }\n"
"</style>\n";
}
  

void Print_HTMLHEAD(std::ostream &out)
{
  out <<
"<html>\n"
"<head>\n"
"<title>CoCoALib Task Table</title>\n"
"<link rel=\"CoCoA Home Page\" href=\"http://cocoa.dima.unige.it\"/>\n";
  Print_CSS(out);
  out << "</head>" << endl << endl;
}


void Print_AllCompletedTasksLink(std::ostream &out)
{
    out <<
      "\n<p>\n"
      "<a href=\"https://cocoa.dima.unige.it/redmine/projects/cocoalib/issues?query_id=7\"><i>Completed tasks (on Redmine)</i></a>\n"
      "<a href=\"CoCoALib-CompletedTasks.html\"><i>Completed tasks (old table)</i></a>\n"
      "</p>\n";
}

void Print_HTMLTAIL(std::ostream &out)
{
  out <<
"\n</div>\n\n"
"<script type=\"text/javascript\">\n"
"var gaJsHost = ((\"https:\" == document.location.protocol) ? \"https://ssl.\" : \"http://\\\n"
"www.\");\n"
"document.write(unescape(\"%3Cscript src='\" + gaJsHost + \"google-analytics.com/ga.js' \\" 
"\n"
"type='text/javascript'%3E%3C/script%3E\"));\n"
"</script>\n"
"<script type=\"text/javascript\">\n"
"try {\n"
"var pageTracker = _gat._getTracker(\"UA-11245099-8\");\n"
"pageTracker._trackPageview();\n"
"} catch(err) {}</script>\n"
"</body>\n"
"</html>\n";
}

void Print_InternalLinks(std::ostream &out, const std::string& link)
{
  out <<
    "<!-- ============================================================ -->\n"
    "<center><a name=\""<<link<<"\"><table cellpadding=\"5\"><tbody><tr>\n"
    "<th><a href=\"#task-table\">Task table</a></th>\n"
    "<th><a href=\"#task-descriptions\">Task descriptions</a></th>\n"
    "<th><a href=\"#RecentlyCompletedTasks\">Recently completed tasks</a></th>\n"
    "<th><a href=\"https://cocoa.dima.unige.it/redmine/projects/cocoalib/issues?query_id=7\"><i>Completed tasks (redmine)</i></a></th>\n"
    "<th><a href=\"CoCoALib-CompletedTasks.html\"><i>Completed tasks (old table)</i></a></th>\n"
    "</tr></tbody></table></a></center>\n"
    "\n";
}

void Print_LinkTable_Preamble(std::ostream &out)
{
  out <<
"<!-- ============================================================ -->\n"
"<!-- ============================================================ -->\n"
"     <h2>The Task Table</h2>\n"
"<!-- ============================================================ -->\n"
"<!-- ============================================================ -->\n"
"\n"
"<div>\n"
"This table summarises the state of current and imminent tasks in\n"
"CoCoALib.\n"
"The first tasks are considered to be the most pressing.\n"
"Full descriptions of the tasks are given <a href=\"#task-descriptions\">below</a>.\n"
"\n"
"<br>\n"
"<font color=\"blue\">All time estimates and completion dates are to be considered approximate.</font> \n"
"</div>\n"
"\n"
"<center>\n"
"<table cellpadding=\"5\">\n"
"<thead>\n"
"<tr><!-- ========================= --></tr>\n"
"<tr>\n"
"<th>Task</th>\n"
"<th>Who</th>\n"
"<th><a name=\"TblEntryStatus\" href=\"#TaskDescrStatus\">Status</a></th>\n"
"<th><a name=\"TblEntryProgress\" href=\"#TaskDescrProgress\">Progress</a></th>\n"
"<th>Last news</th>\n"
"<th>Completion</th>\n"
"</tr>\n"
"<tr><!-- ========================= --></tr>\n"
"</thead>\n"
"\n"
"<!-- ============================================================ -->\n"
"<tbody>\n";
}
  
void Print_LinkTable_Postamble(std::ostream &out)
{
  out <<
"</tbody></table>\n"
"</center>\n\n";
}
  
void Print_TasksDescription_Preamble(std::ostream &out)
{
  out <<
"<!-- ============================================================ -->\n"
"<h2>Task descriptions (Alphabetical order)</h2>\n"
"<!-- ============================================================ -->\n"
"\n"
"<div>\n"
"\n"
"This section contains a short description of each task.\n"
"<!-- Each task has a numerical ID distinct from those of other tasks; the -->\n"
"<!-- order of the numerical IDs has no particular meaning. -->\n"
"A time estimate for a task assumes \"full-time\" working; as ever for\n"
"software, time estimates can be wildly inaccurate (usually hopelessly\n"
"optimistic).\n\n";
 }

void Print_TasksDescription_Postamble(std::ostream &out)
{
  out <<
"</div>\n"
"\n"
"\n"
"\n"
"<a name=\"TaskDescrProgress\"></a>\n"
"<div>\n"
"<h3>Development stages</h3>\n"
"\n"
"To describe the progress of each subtask (contributing C++ code to the\n"
"library), I suggest using the following scheme for identifying what\n"
"has been done, and what is still to be done.  The intention is that\n"
"normally all lower numbered stages be completed before progressing to\n"
"the higher numbered stages.  For other types of task, the task description\n"
"should give an indication of progress.\n"
"\n"
"<ul>\n"
"<li>(1) subproject start (incl specification of subproject)</li>\n"
"<br/>\n"
"<li>(2a) abstract algm design, choice of data-structures, etc.</li>\n"
"<li>(2b) initial test set</li>\n"
"<li>(2c) early prototype, mostly works, maybe some bugs</li>\n"
"<br/>\n"
"<li>(3a) late prototype, fully working, almost no bugs</li>\n"
"<li>(3b) comprehensive test cases</li>\n"
"<br/>\n"
"<li>(4a) final C++ version (incl. relevant build commands)</li>\n"
"<li>(4b) maybe extra tests</li>\n"
"<li>(4c) portability and quality control issues</li>\n"
"<li>(4d) documentation</li>\n"
"<br/>\n"
"<li>(5a)-(5d) independent verification of stages (4a)-(4d)</li>\n"
"</ul>\n"
"</div>\n"
"\n"
"\n"
"<a name=\"TaskDescrStatus\"></a>\n"
"<div>\n"
"<h3>Status</h3>\n"
"\n"
"For each task, the table indicates the current \"status\".  Here are the main choices for the status:\n"
"<ul>\n"
"<li><b>active</b> being worked on currently</li>\n"
"<li><b>inactive</b> not being worked on currently (but work should begin fairly soon)</li>\n"
"<li><b>waiting</b> progress suspended until some other task has progressed enough</li>\n"
"<li><b>resting</b> temporarily suspended (e.g. holiday)</li>\n"
"<li><b>other</b> explanation of actual status is with task description</li>\n"
"</ul>\n"
"</div>\n\n";
}

void Print_CompletedTasks_Preamble(std::ostream &out, const std::string& title)
{
  out <<
    "<!-- ============================================================ -->\n";
  if (title=="all")
    out <<
      "<h2>All Completed Tasks \n"
      "<a href=\"https://cocoa.dima.unige.it/redmine/projects/cocoalib/issues?query_id=7\">new on redmine</a>\n"
      "(<a href=\"CoCoALib-Tasks.html\">Task Table</a>)</h2>\n";
  else if (title=="recent")
    out <<
      "<h2>Completed Tasks \n"
      "(<a href=\"CoCoALib-CompletedTasks.html\">all completed tasks</a>)</h2>\n";
  else
    PrintWarning(title, "unrecognized arg for Print_CompletedTasks_Preamble");
  out <<
    "<!-- ============================================================ -->\n"
    "\n"
    "<div>\n\n"
    "Tasks we believe to have been completed are moved from the table to\n"
    "here.\n"
    "<br/>\n"
    "The dates are only approximate and this table is updated less often\n"
    "than it should be (too much work on CoCoALib ;-).\n"
    "\n\n";
}

void Print_CompletedTasks_Postamble(std::ostream &out)
{
  out <<
"</div>\n"
"<!-- ============================================================ -->\n"
"\n";
 }

//----------------------------------------------------------------------


std::string ClosedTag(const std::string& tag)
{
  string ClosedXMLTag = tag;
  ClosedXMLTag.replace(0, 1, "</");
  size_t SpacePos = ClosedXMLTag.find(" ");
  if (SpacePos != string::npos)  return ClosedXMLTag.substr(0, SpacePos)+">";
  return ClosedXMLTag;
}


void PrintAuthors(std::ostream &out, const index::entry& task, const std::string& tag)
{
  //  const vector<string>& authors = task.myAuthors();
  string authors = task.myAuthors();
  if (authors.size()==0)
  {
    if (task.myComplDate().compare("2007-04-04") > 0)
      PrintWarning(task.myTitle(), "empty authors list");
    return;
  }
  out << tag;
//   out << authors[0];
//   for (size_t j=1; j < authors.size(); ++j)  out << "+" << authors[j];
  out << authors;  
  out << ClosedTag(tag) << endl;
}


void PrintTasksTable(std::ostream &out, long /*YearBound*/)
{
  const index& index(UniqueIndex());
  string link;
  out << "<h1>CoCoALib Task Table (Genova)</h1>" << endl << endl;
  Print_InternalLinks(out, "task-table");
  Print_LinkTable_Preamble(out);
  for (size_t i=0; i<index.mySize(); ++i)
  {
    if ( !index.IamHighPriority(i) )  continue;
    if ( index.myStatus(i) == "active" )
      out << "<tr class=active>";
    else 
      out << "<tr>";
    link = LinkName(index.myTitle(i));
    out << "<td><a name=\"Tbl:" <<link<< "\" href=\"#Descr:"
        <<link<<"\"\n       >"
        << index.myTitle(i) << "</a></td>" << endl;
    PrintAuthors(out, index.myEntry(i), "<td>");
    out << "<td>" << index.myStatus(i) << "</td>  ";
    out << "<td>Stage " << index.myStage(i) << "</td>  ";
    PrintDate(out, index.myLastNewsDate(i), "<td>");
    out << "  ";
    PrintDate(out, index.myComplDate(i), "<td>");    
    out << "\n</tr>" << endl << endl;
  }
  out << "<!-- === COMPLETED ===================================== -->\n\n";
  for (size_t i=0; i<index.mySize(); ++i)
  {
    if ( index.myStatus(i) != "completed" ) continue;
    if ( IsMoreThan12MonthsAgo(index.myComplDate(i) ) ) continue;
    out << "<tr class=completed>";
    link = LinkName(index.myTitle(i));
//     out << "<td><a name=\"Tbl:" <<link<< "\" href=\"#Descr:"<<link<<"\">"
//         << index.myTitle(i) << "</a></td>" << endl;
          out << "<td><a href=\"#RecentlyCompletedTasks\"\n       >"
              << index.myTitle(i) << "</a></td>" << endl;
    PrintAuthors(out, index.myEntry(i), "<td>");
    out << "<td>COMPLETED</td><td></td><td></td>";
    PrintDate(out, index.myComplDate(i), "<td>");    
    out << endl << "</tr>" << endl << endl;
  }
  out << "<tr><!-- ========================= --></tr>\n";
  out << "<tr><th colspan=6></th></tr>\n";
  out << "<tr><!-- ========================= --></tr>\n";
  out << "<!-- === LESS URGENT ===================================== -->\n";
  for (size_t i=0; i<index.mySize(); ++i)
  {
    if ( index.myStatus(i) == "completed" )  continue;
    if ( index.IamHighPriority(i) )  continue;
    if ( index.myStatus(i) == "active" )
      PrintWarning(index.myTitle(i), "active and low priority????");
    out << "<tr>\n";
    link = LinkName(index.myTitle(i));
    out << "<td><a name=\"Tbl:" <<link<< "\" href=\"#Descr:"
        <<link<<"\"\n       >"
        << index.myTitle(i) << "</a></td>" << endl;
    PrintAuthors(out, index.myEntry(i), "<td>");
    out << "<td>" << index.myStatus(i) << "</td>  ";
    out << "<td>Stage " << index.myStage(i) << "</td>  ";
    PrintDate(out, index.myLastNewsDate(i), "<td>");
    out << "  ";
    PrintDate(out, index.myComplDate(i), "<td>");
    out << "\n</tr>" << endl << endl;
  }
  Print_LinkTable_Postamble(out);
}


void PrintTasksDescription(std::ostream &out)
{
  index index(UniqueIndex());
  index.mySortAlphabetically();
  string link;
  Print_InternalLinks(out, "task-descriptions");
  Print_TasksDescription_Preamble(out);
  for (size_t i=0; i<index.mySize(); ++i)
  {
    if ( index.myStatus(i) == "completed" )  continue;
    out << "<p>";
    link = LinkName(index.myTitle(i));
    out << "<b><a name=\"Descr:" <<link<< "\" href=\"#Tbl:"<<link<<"\">"
        << index.myTitle(i) << "</a></b>" << endl;
    out << index.myDescription(i);
    PrintAuthors(out, index.myEntry(i), "<font color=\"blue\">");
    out << "</p>" << endl << endl;
  }
  Print_TasksDescription_Postamble(out);
}


void PrintCompletedTasks(std::ostream &out, long YearBound)
{
  const index& index(UniqueIndex());
  long y = 0;
  string link;
  if (YearBound!=0)
  {
    Print_InternalLinks(out, "RecentlyCompletedTasks");
    Print_CompletedTasks_Preamble(out, "recent");
  }
  else
  {
    Print_CompletedTasks_Preamble(out, "all");
    PrintYearLinks(out);
  }
  for (size_t i=0; i<index.mySize(); ++i)
  {
    if ( index.myStatus(i) != "completed" )  continue;
    if (YearInDate(index.myComplDate(i)) != y)
    {
      out << "</div>\n";
      y = YearInDate(index.myComplDate(i));
      if (y < YearBound)  return;
      out << 
        "<!-- =========================================================== -->\n"
        "<h2><a name=" << y << ">" << y << "</a></h2>\n"
        "<div>\n\n";
    }
    out << "<p>" << endl;
    PrintDate(out, index.myComplDate(i), "<i>");
    out << "\n<b>" << index.myTitle(i) << "</b>" << endl;
//     link = LinkName(index.myTitle(i));
//     out << "\n<b><a name=\"Descr:" <<link<< "\" href=\"#Tbl:"<<link<<"\">"
//         << index.myTitle(i) << "</a></b>" << endl;
    out << index.myDescription(i);
    PrintAuthors(out, index.myEntry(i), "<font color=\"blue\">");
    out << "</p>" << endl << endl;
    // checks
    if (IsDateInFuture(index.myComplDate(i)))
      PrintWarning(index.myTitle(i), "completion date is in the future");
  }
  Print_CompletedTasks_Postamble(out);
}


const char* XMLFileName()
{
  static const string UniqueCopy("tasks.xml");
  return UniqueCopy.c_str();
}


const index& UniqueIndex()
{
  static index UniqueCopy;
  return UniqueCopy;
}


} // namespace HTMLTasks


int main()
{
  // completed tasks
  {
    ofstream out("CoCoALib-CompletedTasks.html");
    if (!out)
    {
      std::cout << "Cannot use file `CoCoALib-CompletedTasks.html`.  Aborting." << endl;
      abort();
    }
    HTMLTasks::Print_HTMLHEAD(out);
    out << "<body>" << endl << endl;
    HTMLTasks::PrintCompletedTasks(out, 0);
    HTMLTasks::Print_HTMLTAIL(out);
  }
  // HTML
  {
    ofstream out("CoCoALib-tasks.html");
    if (!out)
    {
      std::cout << "Cannot use file `CoCoALib-tasks.html`.  Aborting." << endl;
      abort();
    }
    HTMLTasks::Print_HTMLHEAD(out);
    out << "<body>" << endl << endl;
    HTMLTasks::PrintTasksTable(out, HTMLTasks::YearToday()-1);
    HTMLTasks::PrintTasksDescription(out);
    HTMLTasks::PrintCompletedTasks(out, HTMLTasks::YearToday()-1);
    HTMLTasks::Print_AllCompletedTasksLink(out);
    HTMLTasks::Print_HTMLTAIL(out);
  }
}

// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/doc/CoCoALib-tasks/HTMLTasks.C,v 1.17 2019/09/16 14:24:48 abbott Exp $
// $Log: HTMLTasks.C,v $
// Revision 1.17  2019/09/16 14:24:48  abbott
// Summary: Commented out unused parameter
//
// Revision 1.16  2016/12/08 10:32:17  abbott
// Summary: Minor corrections (about testing whether a stream is bad)
//
// Revision 1.15  2014/09/26 16:11:27  bigatti
// *** empty log message ***
//
// Revision 1.14  2014/09/25 14:29:51  bigatti
// -- localtime instead of localtime_r for windows
//
// Revision 1.13  2014/09/25 13:41:40  bigatti
// -- disabled localtime_r for MINGW
//
// Revision 1.12  2012/05/28 13:57:18  bigatti
// -- added link to redmine
//
// Revision 1.11  2011/07/20 10:09:36  bigatti
// -- added include ctime
//
// Revision 1.10  2011/07/19 17:12:45  bigatti
// -- removed dependency for DateTime in time.H (problems with make -j2)
//
// Revision 1.9  2011/07/19 16:10:14  bigatti
// -- using cocoa DateTime
//
// Revision 1.8  2011/07/06 13:52:33  bigatti
// -- fixed IsDateInFuture
// -- improved IsSubstring
//
// Revision 1.7  2011/05/27 15:09:18  bigatti
// -- added include string.h
//
// Revision 1.6  2011/05/27 14:28:40  bigatti
// -- improved logging
//
// Revision 1.5  2011/05/26 06:31:48  bigatti
// -- commented out "Loading task table"
//
// Revision 1.4  2011/05/17 09:50:51  bigatti
// -- final: fixed output files
// -- warning for not sorted entries
//
// Revision 1.3  2011/05/16 16:53:00  bigatti
// -- simplified authors
//
// Revision 1.2  2011/05/16 15:19:33  bigatti
// -- evolution for tasks descriptions
//
// Revision 1.1  2011/05/16 13:30:35  bigatti
// -- first import
//

