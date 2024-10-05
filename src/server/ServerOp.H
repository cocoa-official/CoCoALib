#ifndef CoCoA_ServerOp_H
#define CoCoA_ServerOp_H

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

#include "CoCoA/SmartPtrIRC.H"
#include <iosfwd>
// using std::istream;
// using std::ostream;
#include <string>
// using std::string;

namespace CoCoA
{
  //----------------------------------------------------------------------
  class ServerOpBase; // forward declaration for SmartPtrIRC

  //----------------------------------------------------------------------
  class ServerOp
  {
  public:
    explicit ServerOp(); // defined later: returns VoidOperation
    explicit ServerOp(ServerOpBase* ptr): mySmartPtr(ptr) {}
    // assignment disabled because SmartPtrIRC has no assignment
    ServerOpBase* operator->() { return mySmartPtr.operator->(); }  ///< Allow member fns to be called.
    const ServerOpBase* operator->() const { return mySmartPtr.operator->(); }  ///< Allow const member fns to be called.
    const ServerOpBase* myRawPtr() const { return mySmartPtr.myRawPtr(); } ///< Used by "downcasting" functions IsVoidOperation
    bool operator==(const ServerOp& f) const { return mySmartPtr==f.mySmartPtr; }

  private: // data members
    SmartPtrIRC<ServerOpBase> mySmartPtr;
  };


  void RegisterOp(const std::string& s, ServerOp o);
  ServerOp& GetOperation(const std::string& s);
  void PrintOperations(std::ostream& out);
  void PrintLibraries(std::ostream& out);
  void CheckOperationRegistry();

  //----------------------------------------------------------------------
  class ServerOpBase: protected IntrusiveReferenceCount
  {
    friend class SmartPtrIRC<ServerOpBase>; // morally "friend ServerOp", so it can alter reference count
    //---------------
  public:
    class LibraryInfo
    {
    public:
      LibraryInfo(const std::string& name, const std::string& version, const std::string& group);
      void myOutputSelf(std::ostream&) const;
      bool operator==(const LibraryInfo&) const;
    private:
      std::string myNameValue; ///< name of the library
      //      unsigned int myVersionMajorValue; ///< library major version
      //      unsigned int myVersionMinorValue; ///< library minor version
      std::string myVersionValue; ///< for now the version is a string (see BuildInfo.H)
      std::string myGroupValue; ///< family of functions with the same flavour
    };
    //---------------
    
  protected:
    ServerOpBase(const LibraryInfo& lib): myLibraryValue(lib) {};
  public:
    virtual ~ServerOpBase() {};
    const LibraryInfo& myLibrary() const  { return myLibraryValue; }
    virtual void myOutputSelf(std::ostream&) const = 0;
    virtual void myReadArgs(std::istream& in, int NumArgs) = 0;
    virtual void myCompute() = 0;
    virtual void myWriteResult(std::ostream&) const = 0;
    virtual void myClear() = 0;

    static const std::string ourVarName4; ///< name of CoCoA4 variable which will contain the result of the computation
  private:
    const LibraryInfo myLibraryValue; ///< the library (or sub-library) implementing this operation
};

  bool IsVoidOperation(const ServerOp& o);
  std::ostream& operator<<(std::ostream& out, const ServerOp& o);
  std::ostream& operator<<(std::ostream& out, const ServerOpBase::LibraryInfo& l);
}

// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/server/ServerOp.H,v 1.2 2022/02/22 20:40:04 abbott Exp $
// $Log: ServerOp.H,v $
// Revision 1.2  2022/02/22 20:40:04  abbott
// Summary: Updated copyright message (redmine 1555)
//
// Revision 1.1  2013/05/27 12:57:39  abbott
// Moved all server-related code into src/server/
//
// Revision 1.5  2008/09/22 16:07:02  bigatti
// -- tested (and fixed) communication with cocoa-4 (number of arguments
//    passed to myReadArgs   )
//
// Revision 1.4  2007/10/30 17:14:11  abbott
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

#endif
