//   Copyright (c)  2007  John Abbott and Anna Bigatti

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

#include "ServerOp.H"
#include "CoCoA/error.H"

#include <algorithm>
using std::find;
#include <iostream>
using std::endl;
using std::clog;
#include <map>
using std::map;
// #include <string>
using std::string;
#include <vector>  // only for PrintLibraries
using std::vector;

namespace CoCoA
{
  const std::string ServerOpBase::ourVarName4 = "MEMORY.PKG.CoCoA5.Result";


  // ---- LibraryInfo ----

  ServerOpBase::LibraryInfo::LibraryInfo(const std::string& name, const std::string& version, const std::string& group):
    myNameValue(name),
    myVersionValue(version),
    myGroupValue(group)
  {}


  void ServerOpBase::LibraryInfo::myOutputSelf(std::ostream& out) const
  {
    out << myNameValue << "-" << myVersionValue << " (" << myGroupValue << ")";
  }


  bool ServerOpBase::LibraryInfo::operator==(const ServerOpBase::LibraryInfo& li) const
  {
    return (myNameValue ==li.myNameValue) &&
      (myVersionValue == li.myVersionValue) &&
      (myGroupValue == li.myGroupValue);
  }
  

  std::ostream& operator<<(std::ostream& out, const ServerOpBase::LibraryInfo& l)
  {
    l.myOutputSelf(out);
    return out;
  }

  // ---- ServerOp ----

  std::ostream& operator<<(std::ostream& out, const ServerOp& o)
  {
    out << "[" << o->myLibrary() << "] \t";
    o->myOutputSelf(out);
    return out;
  }

  // ---- default value for ServerOp ----
  class VoidOperation: public ServerOpBase
  {
  public:
    VoidOperation(): ServerOpBase(ServerOpBase::LibraryInfo("VoidLibrary", "", "")) {};
    ~VoidOperation() {};
    void myOutputSelf(std::ostream& out) const {out << "VoidOperation";}
    void myReadArgs(std::istream& /*in*/, int /*NumArgs*/) { CoCoA_THROW_ERROR("non-existent operator","VoidOperation::myReadArgs"); }
    void myCompute() { CoCoA_THROW_ERROR("non-existent operator","VoidOperation::myCompute"); }
    void myWriteResult(std::ostream& /*out*/) const {}
    void myClear() {}
  };


  ServerOp::ServerOp():
    mySmartPtr(new VoidOperation())
  {}


  bool IsVoidOperation(const ServerOp& o)
  {
    return dynamic_cast<const VoidOperation*>(o.myRawPtr()) != 0;
  }


  //---- registry ------------------------------------------------------
  class registry
  {
  public:
    explicit registry() {}
    ~registry() {}

    friend void RegisterOp(const std::string& s, ServerOp o);
    friend ServerOp& GetOperation(const std::string& s);
    friend void PrintOperations(std::ostream& out);
    friend void PrintLibraries(std::ostream& out);
    friend void CheckOperationRegistry();
  private:
    std::map<std::string, ServerOp> myOperationMap;
    std::string myMultiplyDefinedString;
  };


  registry& GetReferenceToGlobalRegistry()
  {
    static registry GlobalRegistry;
    return GlobalRegistry;
  }


  ServerOp& GetOperation(const std::string& s)
  {
    return GetReferenceToGlobalRegistry().myOperationMap[s];
  }


  void PrintOperations(std::ostream& out)
  {
    out << "Provides the following operations:" << endl;
    map<string, ServerOp>& r = GetReferenceToGlobalRegistry().myOperationMap;
    for (map<string, ServerOp>::const_iterator it=r.begin(); it!=r.end(); ++it)
      out << "  " << it->second << endl;
  }
  

  void PrintLibraries(std::ostream& out)
  {
    out << "Provides operations defined in the following libraries:" << endl;
    // this is not very efficient, but it is called only once and the list is small
    vector<ServerOpBase::LibraryInfo> libs;
    for (map<string, ServerOp>::const_iterator it=GetReferenceToGlobalRegistry().myOperationMap.begin(); 
         it!=GetReferenceToGlobalRegistry().myOperationMap.end();
         ++it)
    {
      const ServerOpBase::LibraryInfo& itLib = (it->second)->myLibrary();
      if (find(libs.begin(), libs.end(), itLib) == libs.end())
      {
        libs.push_back(itLib);
        out << "  " << itLib << endl;
      }
    }
  }
  

  void RegisterOp(const std::string& s, ServerOp o)
  {
//     if (s == "")
//       GetReferenceToGlobalRegistry().myMultiplyDefinedString = "cannot set empty string";
//     else if ( !IsVoidOperation(GetOperation(s)) )
//       GetReferenceToGlobalRegistry().myMultiplyDefinedString = s;
//     else
      GetReferenceToGlobalRegistry().myOperationMap.insert(std::make_pair(s, o));
  }


  void CheckOperationRegistry()
  {
    if (GetReferenceToGlobalRegistry().myMultiplyDefinedString != "")
      CoCoA_THROW_ERROR("Multiply defined string: " +
                  GetReferenceToGlobalRegistry().myMultiplyDefinedString,
                  "CheckOperationRegistry");
  }
}

// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/server/ServerOp.C,v 1.3 2022/02/22 20:40:04 abbott Exp $
// $Log: ServerOp.C,v $
// Revision 1.3  2022/02/22 20:40:04  abbott
// Summary: Updated copyright message (redmine 1555)
//
// Revision 1.2  2020/06/17 15:49:31  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.1  2013/05/27 12:57:39  abbott
// Moved all server-related code into src/server/
//
// Revision 1.7  2009/10/29 18:42:13  abbott
// Added necessary include directive for <algorithm>
//
// Revision 1.6  2008/11/13 12:13:17  bigatti
// -- NumArgs is now used in every myReadArgs, but usually just a dumb check
//
// Revision 1.5  2008/09/22 16:07:02  bigatti
// -- tested (and fixed) communication with cocoa-4 (number of arguments
//    passed to myReadArgs   )
//
// Revision 1.4  2007/10/30 17:14:07  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.3  2007/09/27 10:31:42  bigatti
// -- added function PrintLibraries
//
// Revision 1.2  2007/09/25 16:28:31  bigatti
// -- ServerOp includes infos about the library it is defined in
// -- CoCoAServer no longer prints its own version
// -- CoCoAServer prints all offered operations
// -- CoCoAServer can change stat_level (Max's verbosity) at runtime
//
// Revision 1.1  2007/04/27 14:54:22  bigatti
// -- content of CoCoAServer.C split into dedicated files
// -- new registration mechanism (through include "RegisterServerOps.H")
//
