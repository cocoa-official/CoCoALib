// Copyright (c) 2012  Anna Bigatti
//------------------------------------------------------------------ 
// This program is not a CoCoALib example:
// it is just some simple code to extract from the files ex-***.C
// the ShortDescription and the LongDescription strings
//------------------------------------------------------------------ 
// You are free to use any part of this example in your own programs.

#include<iostream>
#include<string>
#include<cstdio>  // for getchar
#include<cstdlib>  // for exit

using namespace std;


void CopyString(void)
{
  int NumQuotes = 0;
  char c;

  while ( (c=getchar()) != '=' && c != EOF ) {}
  if (c == EOF) exit(1);
  while ( (c=getchar()) != '\"' && c != ';' && c != EOF ) {}
  if (c != '\"') exit(1);
  while ( (c != ';') || (NumQuotes%2!=0) )
  {
    switch ( c )
    {
    case EOF: exit(1);
    case '\"': ++NumQuotes; break;
    case '\n':
      if (NumQuotes%2==0) { while ( (c=getchar()) != '\"' && c != EOF){} if (c == EOF) exit(1); ++NumQuotes;}
      else cout << c;
      break;
    case '\\':
      switch ( c=getchar() )
      {
      case 'n':  cout << '\n'; break;
      case '\"': cout << '\"'; break;
      default:   cout << '\\' << c;
      }
      break;
    default:   cout << c;
    }
    c = getchar();
  }
  cout << endl;
}


int main(int argc, char**argv)
{
  for ( int i=1 ; i<argc ; ++i )
    if (argv[i]!=string("ShortDescription") && argv[i]!=string("LongDescription"))
    {
      cerr << "***ERROR***  Unrecognised flag \"" << argv[i] << "\"" << endl;
      exit(1);
    }

  string str;
  while ( true )
  {
    cin >> str;
    if ( cin.eof() ) break;
    if ( str == argv[1])
    {
      CopyString();
      cout << endl;
    }
    if ( str == "program()")
      break;
  }
}

