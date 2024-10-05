//   Copyright (c)  2009 John Abbott,  Anna M Bigatti
//   Orig author:  2009 Giovanni Lagorio (lagorio@disi.unige.it)
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

#include "CoCoA/error.H"
#include "CoCoA/GlobalManager.H"
//#include "CoCoA/library.H" // already included by AST->Parser->Interpreter


#include "globals.H"
#include "Interpreter.H"
#include "CompilationDate.H"
#include "Banner.H"
#include "OnlineHelp.H"
#include "VersionInfo.H"

#include <boost/filesystem.hpp>
#include <boost/scope_exit.hpp>
////#include <cstdlib>
////#include <fstream>
#include <iostream>
////#include <sys/types.h>
#include <dirent.h>

using namespace std;
using namespace boost;
using namespace CoCoA;
using namespace CoCoA::LexerNS;
using namespace CoCoA::ParserNS;
using namespace CoCoA::InterpreterNS;
using namespace CoCoA::AST;


namespace
{

// Banner not needed for GUI, so hide it C5IDE is defined.
#ifndef C5IDE
  void banner(const string &message)
  {
    cout <<
      "\n\n"
      "##########################################################################\n"
         << message << '\n' <<
      "##########################################################################\n"
         << endl;
  }
#endif

  void printUsage(const char * const progname)
  {
    cout << "Usage:  " << "CoCoAInterpreter" << " [<option> ... <option>] [<filename> ... <filename>]\n\n"
      "If no input filenames are given, then an interactive session is started.\n\n"
      "Options:\n"
      "  --do-nothing         Exit immediately (may not be used with other options)\n"
      "  --enable-system-cmd  Enable use of \"SystemCommand\" to run arbitrary shell commands\n"
      "  --fullCoCoALibError  Enables a verbose error reporting for problems involving CoCoALib\n"
      "  --help               Prints this help\n"
      "  --interactive        Starts an interactive session after the input files have been processed\n"
      "  --no-preamble        No initial banner and loading info\n"
      "  --no-prompt          Suppresses the prompt (silently prevails over --prompt)\n"
      "  --no-readline        Do not use readline (and its history mechanism)\n"
      "  --packageDir <dir>   Sets the package directory (default: \"packages\")\n"
      "  --prompt <prompt>    Sets the prompt suffix (default: \"#\")\n"
      "  --stacksize N        Sets the stacksize to N (must be >= 100)\n"
      "  --status-log FILE    Output status log messages to FILE\n"
      "  --version            Prints out version number then exits (other options are effectively ignored)\n"
      "  -w0, -w1, -w2, -w3   Sets the warning level: -w0 none, -w1 low, -w2 normal (the default), -w3 pedantic\n"
      "  -dwc5                Don't warn about features that may be removed in the final CoCoA 5\n\n"
         << "Executable is " << progname << "\n";
  }

  const string cocoaExtension(".cpkg5");
  const string::size_type cocoaExtLen = cocoaExtension.length();

  int findPackages(const string &pkgDir, const int RootLen, vector<string> &packageList, std::ostream& StrOut)
  {
    DIR *dir = opendir(pkgDir.c_str());
    if (!dir) {
      perror(pkgDir.c_str());
      return EXIT_FAILURE;
    }
    while (struct dirent *entry = readdir(dir))
    {
      if (*entry->d_name=='.')
        continue;
      const string filename(entry->d_name);
      const string::size_type fLen = filename.length();
      const string fullname(pkgDir+"/"+filename);
      struct stat buf;
      if (stat(fullname.c_str(), &buf))
      {
        perror(fullname.c_str());
        return EXIT_FAILURE;
      }
      if (S_ISDIR(buf.st_mode))
      {
        if (filename == "CVS") continue;
        if (findPackages(fullname, RootLen, packageList, StrOut))
          return EXIT_FAILURE;
        continue;
      }
      if (fLen<=cocoaExtLen || filename.substr(fLen-cocoaExtLen, cocoaExtLen)!=cocoaExtension)
      {
        if (//printPreamble &&
          fullname.substr(RootLen, fullname.length()-1) != "init.cocoa5")
          StrOut << "Ignoring file  : "
                 << fullname.substr(RootLen, fullname.length()-1) << endl;
        continue;
      }
      packageList.push_back(fullname);
    }
    closedir(dir);
    return EXIT_SUCCESS;
  }

  int findInitFile(const string &pkgDir, vector<string> &packageList)
  {
    DIR *dir = opendir(pkgDir.c_str());
    if (!dir) {
      perror(pkgDir.c_str());
      return EXIT_FAILURE;
    }
    const string fullname(pkgDir+"/init.cocoa5");
    struct stat buf;
    if (stat(fullname.c_str(), &buf)) {
      perror(fullname.c_str());
      return EXIT_FAILURE;
    }
    packageList.push_back(fullname);
    closedir(dir);
    return EXIT_SUCCESS;
  }


  intrusive_ptr<Interpreter> interpreter;

// #ifndef _WINDOWS
//  	extern "C"
// #endif
  // void signalHandler(int)
  // {
  // 	interpreter->controlC = true;
  // }

}


const string DefaultPackageDir1 = "./packages";
const string DefaultPackageDir2 = "../packages";
string packageDir = DefaultPackageDir1;

const string INIT_NOT_FOUND("!!!!! WARNING: CoCoA init file not found !!!!!");
const string PACKAGES_NOT_FOUND("!!!!! WARNING: Standard CoCoA packages not found !!!!!");
const string PACKAGE_AUTOLOAD_LOADING("Loading package: ");
const string PACKAGE_AUTOLOAD_FAILURE_MESSAGE("\n==========================================\n>>> Package auto-loading has failed!!! <<<\n>>>   See messages above the banner    <<<\n>>> The system may respond erratically <<<\n==========================================\n\n");
const string PACKAGE_AUTOLOAD_SKIPPING_PKG_DUE_TO_FAILURE("Due to failure, skipping auto-load of package: ");

int main(int argc, char **argv)
{
//  Call below to sync_with_stdio allows automatic prompt suppression if we compile with g++ [disabled because actual behaviour is platform dependent]
//  cin.sync_with_stdio(false); // see GetlineLineProvider::doReadNextLine in LineProviders.C
  GlobalManager CoCoAFoundations;
  CoCoA::SignalWatcher MonitorSIGINT(SIGINT);
  CoCoA::SignalWatcher MonitorSIGABRT(SIGABRT);
  CoCoA::SignalWatcher MonitorSIGTERM(SIGTERM);
  long MaxStackSize = 5000; // default max stack size
  bool warnAboutCocoa5 = true;
  bool interactive = false;
  bool UseReadline = true; // ignored if CoCoA_WITH_READLINE is not set
  bool fullCoCoALibError = false;
  bool printPreamble = true;
  WarningSeverity warningLevel = WS_NORMAL;
  int paramsIndex = 1;
  bool IsDefaultPackageDir = true;
  ostringstream PkgOut;  // loaded pkg at start: what should we do with this???
  while (paramsIndex<argc) {
    string arg(argv[paramsIndex]);
    if (arg[0]!='-')
      break;
    ++paramsIndex;
    if (arg=="--do-nothing")
    {
      if (argc==2 && paramsIndex==2) return EXIT_SUCCESS;
      cout << "ERROR: --do-nothing must be sole option\n";
      return EXIT_FAILURE;
    }
    if (arg=="--version") // just prints out version "number" and exits immediately
    {
      cout << CoCoAVersion() << endl;
      return EXIT_SUCCESS;
    }
    if (arg=="--packageDir") {
      IsDefaultPackageDir =false;
      if (paramsIndex==argc) {
        cout << "Missing --packageDir argument\n";
        return EXIT_FAILURE;
      }
      packageDir = argv[paramsIndex++];
      continue;
    }
    if (arg=="--enable-system-cmd")
    {
      SystemCommandPermit::EnableCommand();
      continue;
    }
    if (arg=="--stacksize") // "Quick fix" soln to the problem of limited stack size
    {
      // OBSOLESCENT -- because did not work!!
      if (paramsIndex==argc) {
        cout << "Missing --stacksize argument (size in kilobytes)\n";
        return EXIT_FAILURE;
      }
      MaxStackSize = atol(argv[paramsIndex++]); // better to use istringstream???
      if (MaxStackSize < 100) MaxStackSize = 100;
      continue;
    }
    if (arg=="--status-log") // waiting/running status for GUI
    {
      if (GlobalStatusLogStream.is_open()) { cout << "--status-log may be specified only once\n"; return EXIT_FAILURE; }
      if (paramsIndex==argc)
      { cout << "Missing --status-log argument (full path of status file)\n"; return EXIT_FAILURE; }
      const string StatusFile = argv[paramsIndex++];
      GlobalStatusLogStream.open(StatusFile);
      if (!GlobalStatusLogStream) { cout << "status log file not writable\n"; return EXIT_FAILURE; }
      continue;
    }
    if (arg=="--prompt") {
      if (paramsIndex==argc) {
        cout << "Missing --prompt argument\n";
        return EXIT_FAILURE;
      }
      InteractiveLineProvider::promptSuffix = argv[paramsIndex++];
      continue;
    }
    if (arg=="--no-prompt") {
      // Silently override --prompt option!!
      GlobalFlag_SuppressPrompt = true;
      continue;
    }
    if (arg=="--interactive") {
      interactive = true;
      continue;
    }
    if (arg=="--no-readline") {
      UseReadline = false;
      continue;
    }
    if (arg=="--fullCoCoALibError") {
      fullCoCoALibError = true;
      continue;
    }
    if (arg=="--no-preamble") {
      printPreamble = false;
      continue;
    }
    if (arg=="-w0") {
      warningLevel = WS_NONE;
      continue;
    }
    if (arg=="-w1") {
      warningLevel = WS_LOW;
      continue;
    }
    if (arg=="-w2") {
      warningLevel = WS_NORMAL;
      continue;
    }
    if (arg=="-w3") {
      warningLevel = WS_PEDANTIC;
      continue;
    }
    if (arg=="-dwc5") {
      warnAboutCocoa5 = false;
      continue;
    }
    printUsage(*argv);
    if (arg=="-help" || arg=="--help" || arg=="-?" || arg=="--?")
      return EXIT_SUCCESS;
    cout << "\nUnknown option: " << arg << endl;
    return EXIT_FAILURE;
  }

  // Load the (standard) CoCoA packages...
  vector<string> packageList;

  using boost::filesystem::path;
  using boost::filesystem::exists;
  if (IsDefaultPackageDir)
  {
    if (!is_directory(path(DefaultPackageDir1)))
      packageDir = DefaultPackageDir2;
    else
      packageDir = DefaultPackageDir1;
  }
  
  if (IsDefaultPackageDir && !(exists(path(packageDir)) && is_directory(path(packageDir))))
  {
    string LINE; LINE.resize(PACKAGES_NOT_FOUND.size(), '!');
    cout << "*** Cannot access default package directory: " << path(packageDir) << " ***" << endl
         << endl
         << LINE << endl
         << PACKAGES_NOT_FOUND << endl
         << LINE << endl
         << endl;
  }
  else
  {
    if (printPreamble) cout << "\nLoading CoCoA packages from directory:";
    if (printPreamble) cout << "\n  " << packageDir << endl;
    ///// provare anche ..
    if (findPackages(packageDir, packageDir.length()+1, packageList, PkgOut))
      return EXIT_FAILURE;
    findInitFile(packageDir, packageList); // no error if file is missing
    //    if (printPreamble) cout << "...done" << endl;
  }
#ifdef C5IDE
  (void)(interactive); // to avoid compiler warning about unused parameter
  (void)(UseReadline); // to suppress compiler warning about unused variable
  return IDE::launchTheIDE(warningLevel, warnAboutCocoa5, packageList, fullCoCoALibError);
#else // #ifdef C5IDE
  try {
    intrusive_ptr<CppOSTREAM> output(new CppOSTREAM(cout));
    intrusive_ptr<CppISTREAM> input(new CppISTREAM(cin));
    intrusive_ptr<ErrorReporter> errorReporter(new DefaultErrorReporter(warningLevel, output));
#ifdef CoCoA_WITH_READLINE
    intrusive_ptr<LineProvider> lineProvider;
    if (UseReadline)
      lineProvider = new ReadlineLineProvider();
    else
      lineProvider = new GetlineLineProvider();
#else
    (void)(UseReadline); // to suppress compiler warning about unused variable
    intrusive_ptr<LineProvider> lineProvider(new GetlineLineProvider());
#endif
    interpreter = new Interpreter(warnAboutCocoa5, lineProvider, errorReporter, output, /*input,*/ fullCoCoALibError, MaxStackSize);
    BOOST_SCOPE_EXIT( (&interpreter) ) // needs at least 1 arg!!
    {
      interpreter = nullptr;
    } BOOST_SCOPE_EXIT_END
        bool loadOk = true;
    for (const string &fullname: packageList)
    {
      if (loadOk)
      {
        //				cout << PACKAGE_AUTOLOAD_LOADING << fullname << endl;
        //        if (printPreamble)
        PkgOut << PACKAGE_AUTOLOAD_LOADING << fullname.substr(packageDir.length()+1, fullname.length()-1) << endl;
        interpreter->errorReporter->resetErrorCounts();
        interpreter->readAndExecute(fullname, true, true);
        if (interpreter->errorReporter->getErrorCount())
        {
          loadOk = false;
// JAA 2015-07-29          cout << PACKAGE_AUTOLOAD_FAILURE_MESSAGE << endl;
        }
      }
// JAA 2015-07-29     else
// JAA 2015-07-29       cout << PACKAGE_AUTOLOAD_SKIPPING_PKG_DUE_TO_FAILURE << fullname << endl;
    }
//     if (InitFileNotFound)
//       cout << INIT_NOT_FOUND << endl;
    if (paramsIndex==argc) {
      // No filenames given ==> starting interactive session
      PkgOut << "Starting interactive session\n";
      if (printPreamble) cout << CoCoA5Banner() << endl;
      if (!loadOk) cout << PACKAGE_AUTOLOAD_FAILURE_MESSAGE;
      return interpreter->run();
    }
    while (paramsIndex<argc) {
      const string filename(argv[paramsIndex++]);
      if (printPreamble)
        banner("Executing file " + filename);
      interpreter->readAndExecute(filename, true, true);
      if (errorReporter->getErrorCount())
        return EXIT_FAILURE;
    }
    if (interactive)
      return interpreter->run();
  } catch (const Ciao &) {
    /* nop */
  } catch (const InterruptException &) {
    return EXIT_FAILURE;
  } catch (const BaseException &e) {
    cout << e.reason << '\n';
    return EXIT_FAILURE;
  } catch (const ErrorInfo& err) {
    cout << "***ERROR***  UNCAUGHT CoCoA error";
    ANNOUNCE(cout, err);
    return EXIT_FAILURE;
  } catch (const std::exception &e) {
    cout << "Fatal error: " << e.what() << '\n';
    return EXIT_FAILURE;
  } catch(...) {
    cout << "***ERROR***  UNCAUGHT UNKNOWN EXCEPTION" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
#endif // #ifdef C5IDE
}

