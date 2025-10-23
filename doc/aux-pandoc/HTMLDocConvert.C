#include <iostream>
#include <fstream>
#include <string>


using namespace std;

int main(int argc, char**argv)
{
  char c;
  if (argc!=2) 
  {
    cout << argc << endl;
    cout << string(argv[1]) << endl;
    cout << "arg must be file name without extension" << endl;
    exit(1);
  }
  string FileBaseName = argv[1];
  string InName = FileBaseName + ".txt";
  string OutName = FileBaseName + ".conv.txt";

  ifstream fin;
  ofstream fout;
  
  fin.open(InName.c_str(), ios::in);
  fout.open(OutName.c_str(), ios::out);
  
  char CurrentLine[2][256];
  string CurrLineStr, AuthorLineStr, AuthorStr;

  getline(fin, CurrLineStr); // title
  fout << CurrLineStr << endl;
  getline(fin, AuthorLineStr); // author
  fout << AuthorLineStr << endl;
  getline(fin, CurrLineStr); // gnu
  fout << CurrLineStr << endl;
  getline(fin, CurrLineStr); // include
  fout << CurrLineStr << endl;
  getline(fin, CurrLineStr); // TeXTitle

  cout << AuthorLineStr << endl;
  AuthorStr = AuthorLineStr.substr(21);
  cout << AuthorStr << endl;
  fout << "      TeXTitle{" << FileBaseName << " (" << AuthorStr << ")}" << endl;

  getline(fin, CurrLineStr);
  for (int i=0 ; !fin.eof() ; ++i )
  {
    cout << i << " ";
    fout << CurrLineStr << endl;
    getline(fin, CurrLineStr);
  }

  fin.close();
  fout.close();
}
