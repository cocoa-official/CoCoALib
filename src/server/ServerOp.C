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
