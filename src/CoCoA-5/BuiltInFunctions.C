//   Copyright (c) 2010-2017 Giovanni Lagorio, John Abbott, Anna M. Bigatti
//   Authors: 2010-2011 Giovanni Lagorio
//   Authors: 2011-2017 John Abbott, Anna M. Bigatti
//
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

#include "BuiltInFunctions.H"

#include "CompilationDate.H"
#include "OnlineHelp.H"
#include "VersionInfo.H"
#include "globals.H"


#include <boost/iostreams/filter/newline.hpp>
// #include <boost/version.hpp> in ExternalLibs.C
#include <cstdlib>
// using C++ system command
#include <thread>
// using std::this_thread::sleep_for

using namespace std;
using namespace boost;
using namespace boost::iostreams;
using namespace CoCoA::AST;
using namespace CoCoA::LexerNS;
using namespace CoCoA::ParserNS;

namespace CoCoA {
namespace InterpreterNS {

//----------------------------------------------------------------------
// builtIns is essentially a global variable
//   but with better controlled initialization.
std::vector<NameFunPair>& builtIns()
{
  static std::vector<NameFunPair> GlobalVar;
  return GlobalVar;
}
//---- initialization ----
void RuntimeEnvironment::initBuiltInFunctions()
{
  for (NameFunPair &p: builtIns())
    this->setTopLevelVar(p.first, p.second, VariableSlot::VSF_SystemProtected);
}
//----------------------------------------------------------------------

//---- BOOL

DECLARE_STD_BUILTIN_FUNCTION(not, 1) {
  intrusive_ptr<BOOL> arg = runtimeEnv->evalArgAs<BOOL>(ARG(0));
  return Value::from(!(arg->theBool));
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(assert, 1) {
  intrusive_ptr<BOOL> arg = runtimeEnv->evalArgAs<BOOL>(ARG(0));
  if (!(arg->theBool)) 
    throw RuntimeException("assertion failed", invocationExpression);
  return VoidValue::theInstance;
}
END_STD_BUILTIN_FUNCTION


//---- SYSTEM
DECLARE_STD_BUILTIN_FUNCTION(SystemCommand, 1)
{
  using std::system;
  if (std::system(nullptr) == 0)
    throw RuntimeException("low-level \"std::system\" command unavailable", invocationExpression);
  if (!SystemCommandPermit::IsEnabled())
    throw RuntimeException("SystemCommand not enabled (see CoCoA interpreter launch options)", invocationExpression);
  intrusive_ptr<STRING> arg = runtimeEnv->evalArgAs<STRING>(ARG(0));
  std::cout << flush; // required according to man page at cppreference.com
  return Value::from(static_cast<long>(std::system(arg->theString.c_str())));
}
END_STD_BUILTIN_FUNCTION


//---- ZIP
//---- ex:   ZipRead("prova.zip", "tests/try_upon.cocoa5");
#if 0  // 2 functions for accessing zip files, currently disabled
DECLARE_STD_BUILTIN_FUNCTION(ZipFileList, 1) {
	const Argument &arg = ARG(0);
	const string filename = runtimeEnv->evalArgAs<STRING>(arg)->theString;
	int err;
	struct zip *zf = ::zip_open(filename.c_str(), ZIP_CHECKCONS, &err);
	if (!zf)
		throw RuntimeException("Cannot read Zip file \""+filename+"\"", arg.exp);
	intrusive_ptr<LIST> result(new LIST);
	const int numEntries = zip_get_num_files(zf);
	for(int a=0; a<numEntries; ++a)
		result->addValue(new STRING(zip_get_name(zf, a, ZIP_FL_UNCHANGED)));
	::zip_close(zf);
	return result;
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(ZipRead, 2) {
	const Argument &argFilename = ARG(0);
	const Argument &argEntryName = ARG(1);
	const string filename = runtimeEnv->evalArgAs<STRING>(argFilename)->theString;
	const string entryName = runtimeEnv->evalArgAs<STRING>(argEntryName)->theString;
	int err;
	struct zip *zf = ::zip_open(filename.c_str(), ZIP_CHECKCONS, &err);
	if (!zf)
		throw RuntimeException("Cannot read Zip file \""+filename+"\"", argFilename.exp);
	const int index = zip_name_locate(zf, entryName.c_str(), /* flags */ 0);
	if (index<0) {
		::zip_close(zf);
		throw RuntimeException("Cannot find an entry named \""+entryName+" inside the Zip file \""+filename+"\"", argEntryName.exp);
	}
	struct zip_file *entry = ::zip_fopen_index(zf, index, /* flags */ 0);
	if (!entry) {
		::zip_close(zf);
		throw RuntimeException("Error reading the entry named \""+entryName+" inside the Zip file \""+filename+"\"", argEntryName.exp);
	}
	string result;
	for(;;) {
		char buf[1024];
		int nRead = ::zip_fread(entry, buf, sizeof(buf));
		if (nRead<=0) {
			::zip_fclose(entry);
			::zip_close(zf);
			return new STRING(result);
		}
		result.append(buf, buf+nRead);
	}
}
END_STD_BUILTIN_FUNCTION
#endif


//---- LIST

DECLARE_ARITYCHECK_FUNCTION(NewList) { return (1<=nArg) && (nArg<=2); }
DECLARE_BUILTIN_FUNCTION(NewList) {
	invocationExpression->checkNumberOfArgs(1,2);
  //	intrusive_ptr<INT> N = runtimeEnv->evalArgAs<INT>(ARG(0));
  const long n = runtimeEnv->evalArgAsLong(ARG(0));
  if (n<0) throw RuntimeException("invalid length", ARG(0).exp);
	intrusive_ptr<LIST> returnValue = new LIST();
  intrusive_ptr<RightValue> val = INT::zero;
  if (invocationExpression->args.size()==2)
    val = runtimeEnv->evalArgAs<RightValue>(ARG(1));
  for(long i=0; i<n; ++i)  returnValue->addValue(val);
	return returnValue;
}

DECLARE_STD_BUILTIN_FUNCTION(IsEmpty, 1)
{
  intrusive_ptr<LIST> L = runtimeEnv->evalArgAs<LIST>(ARG(0));
  if (L->size() == 0)
    return BOOL::trueValue;
  return BOOL::falseValue;
}
END_STD_BUILTIN_FUNCTION


DECLARE_ARITYCHECK_FUNCTION(first) { return (1<=nArg) && (nArg<=2); }
DECLARE_BUILTIN_FUNCTION(first) {
	invocationExpression->checkNumberOfArgs(1,2);
	intrusive_ptr<LIST> list = runtimeEnv->evalArgAs<LIST>(ARG(0));
  if (invocationExpression->args.size()==1)
  {
    if (list->size() == 0) throw RuntimeException("list is empty", ARG(0).exp);
    return list->getValue(0);
  }
  //	intrusive_ptr<INT> N = runtimeEnv->evalArgAs<INT>(ARG(1));
	long n = runtimeEnv->evalArgAsLong(ARG(1));
  if (n<0) throw RuntimeException("invalid length", ARG(1).exp);
  if ((unsigned long)n > list->size())
    throw RuntimeException("value greater than list length", ARG(1).exp);
	intrusive_ptr<LIST> returnValue = new LIST();
	for(long i=0; i<n; ++i)  returnValue->addValue(list->getValue(i));
	return returnValue;
}

DECLARE_ARITYCHECK_FUNCTION(last) { return (1<=nArg) && (nArg<=2); }
DECLARE_BUILTIN_FUNCTION(last) {
	invocationExpression->checkNumberOfArgs(1,2);
	intrusive_ptr<LIST> list = runtimeEnv->evalArgAs<LIST>(ARG(0));
  ///////// here
  LIST::ContainerType::size_type s = list->size();
  if (invocationExpression->args.size()==1) 
  if (invocationExpression->args.size()==1)
  {
    if (list->size() == 0) throw RuntimeException("list is empty", ARG(0).exp);
    return list->getValue(s-1);
  }
  //	intrusive_ptr<INT> N = runtimeEnv->evalArgAs<INT>(ARG(1));
  long n = runtimeEnv->evalArgAsLong(ARG(1));
  if (n<0) throw RuntimeException("invalid length", ARG(1).exp);
  if ((unsigned long)n>s)
    throw RuntimeException("value greater than list length", ARG(1).exp);
	intrusive_ptr<LIST> returnValue = new LIST();
	for(unsigned long i=s-n; i<s; ++i)  returnValue->addValue(list->getValue(i));
	return returnValue;
}

DECLARE_ARITYCHECK_FUNCTION(sum) { return (1<=nArg) && (nArg<=2); }
DECLARE_BUILTIN_FUNCTION(sum) {
	invocationExpression->checkNumberOfArgs(1,2);
	const Argument &a0 = ARG(0);
	intrusive_ptr<LIST> list = runtimeEnv->evalArgAs<LIST>(a0);
	LIST::ContainerType::size_type size = list->size();
  intrusive_ptr<RightValue> result = INT::zero;
  if (invocationExpression->args.size()==1)
  {
    if (size==0) throw RuntimeException("empty sum", ARG(0).exp);
    result = list->getValue(--size);
  }
  else
    result = runtimeEnv->evalArgAs<RightValue>(ARG(1));
  for(LIST::ContainerType::size_type a=size; a!=0; /**/)
    try	{
			result = runtimeEnv->binaryOperatorDispatch(list->getValue(--a), result, RuntimeEnvironment::opPlusMap, a0.exp->getBegin(), a0.exp->getEnd());
	} catch (const InterruptException &) {
		throw;
	} catch (const RuntimeException &) {
		throw RuntimeException("Some elements have incompatible types for the sum", a0.exp);
	}
	return result;
}

DECLARE_ARITYCHECK_FUNCTION(product) { return (1<=nArg) && (nArg<=2); }
DECLARE_BUILTIN_FUNCTION(product) {
	invocationExpression->checkNumberOfArgs(1,2);
	const Argument &a0 = ARG(0);
	intrusive_ptr<LIST> list = runtimeEnv->evalArgAs<LIST>(a0);
	LIST::ContainerType::size_type size = list->size();
	intrusive_ptr<RightValue> result = INT::one;
  if (invocationExpression->args.size()==1)
  {
    if (size==0)  throw RuntimeException("empty product", ARG(0).exp);
    result = list->getValue(--size);
  }
  else
    result = runtimeEnv->evalArgAs<RightValue>(ARG(1));
  for(LIST::ContainerType::size_type a=size; a!=0; /**/)
		try	{
			result = runtimeEnv->binaryOperatorDispatch(list->getValue(--a), result, RuntimeEnvironment::opStarMap, a0.exp->getBegin(), a0.exp->getEnd());
	} catch (const InterruptException &) {
		throw;
	} catch (const RuntimeException &) {
		throw RuntimeException("Some elements have incompatible types for the product", a0.exp);
	}
	return result;
}

DECLARE_ARITYCHECK_FUNCTION(concat) { (void)nArg; return true; }
DECLARE_BUILTIN_FUNCTION(concat) {
	intrusive_ptr<LIST> returnValue = new LIST();
	for (const Argument &arg: invocationExpression->args)
        {
          intrusive_ptr<LIST> list = runtimeEnv->evalArgAs<LIST>(arg);
          int size = list->size();
          for(int a=0; a<size; ++a)
            returnValue->addValue(list->getValue(a));
	}
	return returnValue;
}

DECLARE_STD_BUILTIN_FUNCTION(ConcatLists, 1) {
	intrusive_ptr<LIST> listOfList = runtimeEnv->evalArgAs<LIST>(ARG(0));
	intrusive_ptr<LIST> result(new LIST);
	const LIST::ContainerType::size_type nLists = listOfList->size();
	for (LIST::ContainerType::size_type a=0; a<nLists; ++a) {
		intrusive_ptr<LIST> l = dynamic_pointer_cast<LIST>(listOfList->getValue(a));
		if (!l)
			throw RuntimeException("The argument is not a list of lists", ARG(0).exp);
		const LIST::ContainerType::size_type len = l->size();
		for(LIST::ContainerType::size_type b=0; b<len; ++b)
			result->addValue(l->getValue(b));
	}
	return result;
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(append, 2) {
	intrusive_ptr<RightValue> elem = runtimeEnv->evalArgAs<RightValue>(ARG(1));
	intrusive_ptr<LIST> list = runtimeEnv->obtainUnshared<LIST>(ARG(0));
	list->addValue(elem);
	return VoidValue::theInstance;
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(reversed, 1) {
	intrusive_ptr<LIST> list = runtimeEnv->evalArgAs<LIST>(ARG(0));
	intrusive_ptr<LIST> returnValue = new LIST();
	for(LIST::ContainerType::size_type a=list->size(); a!=0; /**/)
		returnValue->addValue(list->getValue(--a));
	return returnValue;
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(reverse, 1) {
	intrusive_ptr<LIST> list = runtimeEnv->obtainUnshared<LIST>(ARG(0));
	const LIST::ContainerType::size_type size = list->size();
	if (size) {
		LIST::ContainerType::size_type left = 0, right = size-1;
		while(left < right) {
			const intrusive_ptr<RightValue> tmp(list->getValue(left));
			list->setValue(left++, list->getValue(right));
			list->setValue(right--, tmp);
		}
	}
	return VoidValue::theInstance;
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(remove, 2) {
	const Argument indexArg = ARG(1);
	intrusive_ptr<RightValue> shouldBeAnIndex = runtimeEnv->evalArgAs<RightValue>(indexArg);
	intrusive_ptr<LIST> list = runtimeEnv->obtainUnshared<LIST>(ARG(0));
	LIST::ContainerType::size_type index = list->checkIndex(shouldBeAnIndex, indexArg.exp);
	list->removeValue(index);
	return VoidValue::theInstance;
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(insert, 3) {
	const Argument &indexArg = ARG(1);
	intrusive_ptr<RightValue> shouldBeAnIndex = runtimeEnv->evalArgAs<RightValue>(indexArg);
	intrusive_ptr<RightValue> elem = runtimeEnv->evalArgAs<RightValue>(ARG(2));
	intrusive_ptr<LIST> list = runtimeEnv->obtainUnshared<LIST>(ARG(0));
  //	LIST::ContainerType::size_type index = list->checkIndex(shouldBeAnIndex, indexArg.exp);
	LIST::ContainerType::size_type index = list->checkIndexOrEnd(shouldBeAnIndex, indexArg.exp);
	list->insertValue(index, elem);
	return VoidValue::theInstance;
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(ConcatStrings, 1) {
	intrusive_ptr<LIST> listOfString = runtimeEnv->evalArgAs<LIST>(ARG(0));
        string str;
	const LIST::ContainerType::size_type n = listOfString->size();
	for (LIST::ContainerType::size_type a=0; a<n; ++a) {
		intrusive_ptr<STRING> s = dynamic_pointer_cast<STRING>(listOfString->getValue(a));
		if (!s)
			throw RuntimeException("The argument is not a list of strings", ARG(0).exp);
                str += s->theString;
	}
	return Value::from(str);
}
END_STD_BUILTIN_FUNCTION




//---- RECORD
DECLARE_STD_BUILTIN_FUNCTION(fields, 1) {
	intrusive_ptr<RECORD> r = runtimeEnv->evalArgAs<RECORD>(ARG(0));
	return r->fieldNames();
}
END_STD_BUILTIN_FUNCTION


//---- I/O
DECLARE_STD_BUILTIN_FUNCTION(sprint, 1) {
	const Argument &arg = ARG(0);
	intrusive_ptr<OutputStringStreamValue> outs(new OutputStringStreamValue);
	outs->print(runtimeEnv, runtimeEnv->evalArgAs<RightValue>(arg), arg.exp);
	return outs->close(runtimeEnv, invocationExpression);
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(GetErrMesg, 1) {
	const Argument &arg = ARG(0);
	intrusive_ptr<ERROR_MESG> err = runtimeEnv->evalArgAs<ERROR_MESG>(arg);
	return new STRING(err->message);
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(close, 1) {
  int which;
  intrusive_ptr<RightValue> v = runtimeEnv->evalArgAsT1orT2<OSTREAM, ISTREAM>(ARG(0), which);
  switch (which) {
  case 1: 
    return intrusive_ptr_cast<OSTREAM>(v)->close(runtimeEnv, invocationExpression);

  case 2:
    return intrusive_ptr_cast<ISTREAM>(v)->close(runtimeEnv, invocationExpression);

  default: /*NEVER GET HERE*/return Value::from(false); // just to keep the compiler quiet
  }
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(OpenOString, 0) {
	(void)runtimeEnv; // to avoid a compilation warning
	return new OutputStringStreamValue();
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(GetEnv, 1) {
	const Argument &arg = ARG(0);
	intrusive_ptr<STRING> varName = runtimeEnv->evalArgAs<STRING>(arg);
	char *value = ::getenv(varName->theString.c_str());
	if (!value)
		return STRING::empty;
	return new STRING(string(value));
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(OpenOFile, 1) {
	const Argument &arg = ARG(0);
	intrusive_ptr<STRING> filename = runtimeEnv->evalArgAs<STRING>(arg);
	return new OutputFileStreamValue(filename, arg.exp);
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(StandardOutput, 0) {
  (void)(runtimeEnv);  // just to keep compiler quiet
  return new CppOSTREAM(std::cout);
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(OpenIFile, 1) {
	const Argument &arg = ARG(0);
	intrusive_ptr<STRING> filename = runtimeEnv->evalArgAs<STRING>(arg);
	return new InputFileStreamValue(filename, arg.exp);
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(OpenIString, 1) {
	const Argument &arg = ARG(0);
	intrusive_ptr<STRING> str = runtimeEnv->evalArgAs<STRING>(arg);
	return new InputStringStreamValue(str);
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(StandardInput, 0) {
  (void)(runtimeEnv);  // just to keep compiler quiet
  return new CppISTREAM(std::cin);
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(OpenSocket, 2)
{
  const string HostName = runtimeEnv->evalArgAs<STRING>(ARG(0))->theString;
  const long PortNum = runtimeEnv->evalArgAsLong(ARG(1));
  if (PortNum < 0 || PortNum > 65535)  throw RuntimeException("invalid port number (must be between 0 and 65535)", ARG(0).exp);
  std::ostringstream PortNumStr;
  PortNumStr << PortNum;
///  using boost::asio::ip::tcp;
  std::shared_ptr<boost::asio::ip::tcp::iostream> ClientSocketPtr(new boost::asio::ip::tcp::iostream);
  ClientSocketPtr->connect(HostName, PortNumStr.str()); // PortNo has to be a string :-/
  if (!*ClientSocketPtr)
    throw RuntimeException("failed to create socket connection", invocationExpression);
ClientSocketPtr->rdbuf()->pubsetbuf(0, 0); // forces a flush() after every output!!

  intrusive_ptr<RECORD> rec(new RECORD);
  rec->setFieldNoCheck("HostName", Value::from(HostName));
  rec->setFieldNoCheck("PortNum", Value::from(PortNum));

  rec->setFieldNoCheck("send", (new OutputBoostSocketStreamValue(ClientSocketPtr)));
  rec->setFieldNoCheck("recv", (new InputBoostSocketStreamValue(ClientSocketPtr)));
  return rec;
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(GetLine, 1) {
	const Argument &arg = ARG(0);
	intrusive_ptr<ISTREAM> istream = runtimeEnv->evalArgAs<ISTREAM>(arg);
        if (istream->eof()) throw RuntimeException("EOF: end of input reached", invocationExpression); // see redmine 1487
        const std::string line = istream->getline();
        return new STRING(line);
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(IsAtEOF, 1) {
	const Argument &arg = ARG(0);
	intrusive_ptr<ISTREAM> istream = runtimeEnv->evalArgAs<ISTREAM>(arg);
        return Value::from(istream->eof());
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(ContentsOfFile, 1) {
	const Argument &arg = ARG(0);
	intrusive_ptr<STRING> filename = runtimeEnv->evalArgAs<STRING>(arg);
	filtering_istream input;
	input.push(newline_filter(newline::posix));
	const file_source fs(filename->theString.c_str());
	if (!fs.is_open())
		throw RuntimeException("Cannot open file \""+filename->theString+"\" for reading.", invocationExpression);
	input.push(fs);
	string wholeFile;
	for(;;) {
		string line;
		getline(input, line);
		if (input.bad())
			throw RuntimeException("Cannot read from \""+filename->theString+"\".", invocationExpression);
		wholeFile += line;
		if (input.eof())
			break;
		wholeFile += '\n';
	}
	return new STRING(wholeFile);
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(substring, 3) {
	const Argument &arg0 = ARG(0);
	const Argument &arg1 = ARG(1);
	const Argument &arg2 = ARG(2);
	const string str = runtimeEnv->evalArgAs<STRING>(arg0)->theString;
        const long first = runtimeEnv->evalArgAsLong(arg1);
        const long length = runtimeEnv->evalArgAsLong(arg2);
        if (first < 1 || first > len(str)) throw RuntimeException("substring start posn not inside string", ARG(1).exp);
        if (length < 0) throw RuntimeException("substring length must be non-neg", ARG(2).exp);
	return new STRING(str.substr(first-1, length));
}
END_STD_BUILTIN_FUNCTION


//---- TAGGED
const string invalidTagError("The tag must be a valid identifier");

DECLARE_STD_BUILTIN_FUNCTION(tagged, 2) {
	const Argument &argValue = ARG(0);
	const Argument &argTag = ARG(1);
	intrusive_ptr<RightValue> v = runtimeEnv->evalArgAs<RightValue>(argValue);
	const string tag = runtimeEnv->evalArgAs<STRING>(argTag)->theString;
	// TODO: must be fixed
	//if (!LexerNS::Lexer::isValidIdentifier(tag))
	//	throw RuntimeException(invalidTagError, argTag.exp);
	return new TaggedValue(v->untagged(), invocationExpression->packageName+"."+tag);
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(untagged, 1) {
	return runtimeEnv->evalArgAs<RightValue>(ARG(0))->untagged();
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(tag, 1) {
	const Argument &argValue = ARG(0);
	intrusive_ptr<RightValue> v = runtimeEnv->evalArgAs<RightValue>(argValue);
	if (intrusive_ptr<TaggedValue> taggedValue = dynamic_pointer_cast<TaggedValue>(v)) {
		const string::size_type l = invocationExpression->packageName.length()+1;
		if (taggedValue->tag.substr(0, l)==invocationExpression->packageName+".")
			return new STRING(taggedValue->tag.substr(l, string::npos));
		return new STRING(taggedValue->tag);
	}
	return STRING::empty;
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(TAGGED, 1) {
	const Argument &arg = ARG(0);
	const string tag = runtimeEnv->evalArgAs<STRING>(arg)->theString;
	if (!LexerNS::Lexer::isValidIdentifier(tag))
		throw RuntimeException(invalidTagError, arg.exp);
	return TYPE::tagType(tag);
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(type, 1) {
	return runtimeEnv->evalArgAs<RightValue>(ARG(0))->getType();
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(types, 0) {
	(void)runtimeEnv; // to avoid a compilation warning
	throw RuntimeException("Cannot return a list of all types because they're infinite (CoCoA 4 lies about this, I don't). You might want to use CurrentTypes() to get a list of currently instantiated types", invocationExpression);
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(CurrentTypes, 0) {
	return runtimeEnv->currentTypes();
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(error, 1) {
	const Argument &arg = ARG(0);
	const string message = runtimeEnv->evalArgAs<STRING>(arg)->theString;
	throw RuntimeException(message, invocationExpression);
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(Error, 1) {
	const Argument &arg = ARG(0);
	const string message = runtimeEnv->evalArgAs<STRING>(arg)->theString;
	throw RuntimeException(message, invocationExpression);
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(ascii, 1) {
  int which;
  intrusive_ptr<RightValue> v = runtimeEnv->evalArgAsT1orT2orT3<STRING, INT, LIST>(ARG(0), which);
  switch (which) {
  case 1: {
    intrusive_ptr<LIST> l(new LIST());
		const string s(RefTo<string>(v));
		for (const unsigned char c: s)
                  l->addValue(new INT(static_cast<int>(c)));
		return l; }
	case 2: {
		long n;
		if (!IsConvertible(n, RefTo<BigInt>(v)) || n<0 || n>255)
			throw RuntimeException("When the argument is an INT, it must be in the range [0..255]", ARG(0).exp);
		return Value::from(string(1, static_cast<unsigned char>(n))); }
	case 3: {
		string retValue;
    intrusive_ptr<LIST> list = dynamic_pointer_cast<LIST>(v);
		LIST::ContainerType::size_type size = list->size();
		for(LIST::ContainerType::size_type a=0; a<size; ++a) {
			intrusive_ptr<INT> N = dynamic_pointer_cast<INT>(list->getValue(a));
			long n;
			if (!N || !IsConvertible(n, N->theBigInt) || n<0 || n>255)
				throw RuntimeException("When the argument is a LIST, it must be a list of INT only in the range [0..255]", ARG(0).exp);
			retValue += static_cast<unsigned char>(n);
		}
		return Value::from(retValue);
  }
  default: throw RuntimeException(ERRORMissingCode(v),invocationExpression); 
  }
}
END_STD_BUILTIN_FUNCTION

//---- general ----

DECLARE_STD_BUILTIN_FUNCTION(TopLevelFunctions, 0) { // GL
	return runtimeEnv->TopLevelFunctions();
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(PackageOf, 1) { // AMB
	const Argument &arg = ARG(0);
	string name = runtimeEnv->evalArgAs<STRING>(arg)->theString;
	return runtimeEnv->PackageOf(name);
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(CocoaPackagePath, 0) { // GL
	(void)runtimeEnv; // keeps the compiler happy (or, at least, silent ;-) )
	return Value::from(packageDir);
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(PkgName, 0) {
	(void)runtimeEnv; // to avoid a compilation warning
	return new STRING(invocationExpression->packageName);
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(aliases, 0) {
	runtimeEnv->topLevelAliases->dump(runtimeEnv->getOutputStream());
	return VoidValue::theInstance;
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(VersionInfo, 0) { // AMB
	(void)runtimeEnv; // keeps the compiler happy (or, at least, silent ;-) )
  intrusive_ptr<RECORD> rec(new RECORD);
  rec->setFieldNoCheck("CoCoALibVersion", Value::from(BuildInfo::version()));
  rec->setFieldNoCheck("CoCoAVersion", Value::from(CoCoAVersion()));
  rec->setFieldNoCheck("CoCoAGitRevision", Value::from(BuildInfo::GitRevision()));
  rec->setFieldNoCheck("CompilationDate", Value::from(CompilationDate()));
  rec->setFieldNoCheck("CompilationFlags", Value::from(BuildInfo::CompilationFlags()));
  rec->setFieldNoCheck("CompilationPlatform", Value::from(BuildInfo::CompilationPlatform()));
  rec->setFieldNoCheck("CoCoA_32BIT_LONG", Value::from(BuildInfo::CoCoA_32BIT_LONG_Flag()));
  rec->setFieldNoCheck("CompilationPreprocessorDefines", Value::from(BuildInfo::CompilationPreprocessorDefines()));
  rec->setFieldNoCheck("Compiler", Value::from(BuildInfo::compiler()));
  rec->setFieldNoCheck("NumBits_int", new INT(std::numeric_limits<unsigned int>::digits));
  rec->setFieldNoCheck("NumBits_long", new INT(std::numeric_limits<unsigned long>::digits));
  rec->setFieldNoCheck("NumBits_SmallExponent_t", new INT(std::numeric_limits<SmallExponent_t>::digits));
  const int NumLibs = len(ExternalLibs());
  intrusive_ptr<LIST> ExtLibList(new LIST);
  intrusive_ptr<RECORD> ExtLib(new RECORD);
  ExtLib->setFieldNoCheck("name", Value::from(std::string("Boost")));
  ExtLib->setFieldNoCheck("version", Value::from(ToString(BOOST_VERSION)));
  ExtLibList->addValue(ExtLib);
  for (int i=0; i < NumLibs; ++i)
  {
    intrusive_ptr<RECORD> ExtLib(new RECORD);
    ExtLib->setFieldNoCheck("name", Value::from(ExternalLibs()[i].myName));
    ExtLib->setFieldNoCheck("version", Value::from(ExternalLibs()[i].myVersion));
		ExtLibList->addValue(ExtLib);
  }
  rec->setFieldNoCheck("ExternalLibs", ExtLibList);
  return rec;
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(SetStackSize, 1) {
//   intrusive_ptr<INT> N = runtimeEnv->evalArgAs<INT>(ARG(0));
//   const BigInt NewSize = N->theBigInt;
//   long n;
//   if (NewSize < 2 || !IsConvertible(n, NewSize)) throw RuntimeException("Ridiculous stack size", ARG(0).exp);
//   return new INT(runtimeEnv->ResizeStack(n, ARG(0).exp));
  long NewSize = runtimeEnv->evalArgAsLong(ARG(0));
  if (NewSize < 2) throw RuntimeException("Ridiculous stack size", ARG(0).exp);
  return Value::from(runtimeEnv->ResizeStack(NewSize, ARG(0).exp));
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(pause, 0) {
#ifdef C5IDE
	runtimeEnv->interpreter->singleStepExecution = true;
#else // #ifdef C5IDE
	(void)runtimeEnv; // keeps the compiler happy (or, at least, silent ;-) )
#endif // #ifdef C5IDE
	return VoidValue::theInstance;
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(SleepFor, 1) {
  int which;
	intrusive_ptr<RightValue> v = runtimeEnv->evalArgAsT1orT2<INT, RAT>(ARG(0), which);
  switch (which) {
  case 1: 
	 {
		long n;
		if (!IsConvertible(n, RefTo<BigInt>(v)) || n<0 || n > 65535)
			throw RuntimeException("Sleep period must be in range 0..65535 (seconds)", ARG(0).exp);
                std::this_thread::sleep_for(std::chrono::seconds(n));
	return VoidValue::theInstance;
         }
	case 2: {
          const BigRat q = RefTo<BigRat>(v);
          if (q < 0 || q > 65535) throw RuntimeException("Sleep period must be in range 0..65535 (seconds)", ARG(0).exp);
          const long n = ConvertTo<long>(CoCoA::floor(1000*q)); // range checked above
          std::this_thread::sleep_for(std::chrono::milliseconds(n));
	return VoidValue::theInstance;
  }
  default: throw RuntimeException(ERRORMissingCode(v),invocationExpression); 
  }
}
END_STD_BUILTIN_FUNCTION

#ifdef C5IDE
DECLARE_STD_BUILTIN_FUNCTION(sleep, 1) {
	const Argument &arg = ARG(0);
	intrusive_ptr<INT> N = runtimeEnv->evalArgAs<INT>(arg);
	long l;
	if (!IsConvertible(l, N->theBigInt) || l<0 || l>1000)
		throw RuntimeException("The argument must be in the range [0..1000]", arg.exp);
        std::this_thread::sleep_for(std::chrono::seconds(static_cast<int>(l)));
	return VoidValue::theInstance;
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(debug, 1) {
	const Argument &arg = ARG(0);
	string filename = runtimeEnv->evalArgAs<STRING>(arg)->theString;
	try {
		runtimeEnv->interpreter->singleStepExecution = true;
		runtimeEnv->interpreter->readAndExecute(filename, false, false);
	} catch (const RuntimeException &) {
		throw;
	} catch (const BaseException &be) {
		throw RuntimeException(be.reason, invocationExpression->getBegin(), invocationExpression->getEnd());
	}
	return VoidValue::theInstance;
}
END_STD_BUILTIN_FUNCTION
#endif // #ifdef C5IDE


//---- manual CoCoAHelp ----

DECLARE_STD_BUILTIN_FUNCTION(WordList, 0) { // AMB
  (void)(runtimeEnv); // just to keep compiler quiet
	return Value::from(OnlineHelp::WordList());
}
END_STD_BUILTIN_FUNCTION


//DECLARE_STD_BUILTIN_FUNCTION(ReloadMan, 0) { // AMB
// variable number of args
DECLARE_ARITYCHECK_FUNCTION(ReloadMan) { return (0<=nArg) && (nArg<=1); }
DECLARE_BUILTIN_FUNCTION(ReloadMan) {  // AMB 2018-03
  invocationExpression->checkNumberOfArgs(0,1);// variable number of args
  std::vector<std::string> FileNames;
  if (invocationExpression->args.size()==1)
  {
    FileNames = runtimeEnv->evalArgAsListOf<STRING>(ARG(0));
  }
	try
  {
    ostringstream os;
    OnlineHelp::ReloadMan(os, FileNames);
    runtimeEnv->getOutputStream()->print(os.str())->flush();
	}
  catch (const std::exception& err)
  {
		throw RuntimeException(err.what(), invocationExpression->getBegin(), invocationExpression->getEnd());
	}
	return VoidValue::theInstance;
} 
// ---^^^--- variable number of args ---^^^---


} // namespace InterpreterNS
} // namespace CoCoA
