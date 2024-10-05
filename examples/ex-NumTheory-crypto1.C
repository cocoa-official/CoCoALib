// Copyright (c) 2006  John Abbott,  Anna M. Bigatti
// Authors: Chiara Della Corte (2010), Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "Interactive program for investigating cryptosystems.  \n"
  "Implementation using BigInt.    \n";

const string LongDescription =
  " DiffieHellman   \n"
  " ElGamal         \n"
  " RSA             \n"
  " Rabin           \n"
  " ====                                                      \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  BigInt InputBigInt(const char* prompt);
  BigInt InputPrime(char ch, const char* str, long MaxNumBits=1000);
  BigInt InputBigIntInInterval(long x, BigInt y, const char* prompt);


  void DiffieHellman()
  {
    cout <<
      "====[  DIFFIE-HELLMAN algorithm for key generation  ]====\n"
      "  Alice and Bob pick together and publicly \n"
      "  a prime  P  and a generator  G  for (FF_P)*\n"
         << endl;
    BigInt P = InputPrime('P', "[PUBLIC]", 90);
    int G = PrimitiveRoot(P);
    cout << "  A generator for (FF_p)* is  G = "<< G << " [PUBLIC]" << endl;
    cout << endl;
    cout << "Alice picks  a [PRIVATE KEY] in 2,...," << P-1 << endl;
    BigInt a = InputBigIntInInterval(2,P-1, "a");
    BigInt A = PowerMod(G, a, P);
    cout << "  then Alice computes  G^a mod P, and makes it public:" << endl
         << "  A = PowerMod(G, a, P) = " << A << " [PUBLIC KEY]" << endl;
    cout << endl;
    cout << "Bob picks  b [PRIVATE KEY] in 2,...," << P-1 << endl;
    BigInt b = InputBigIntInInterval(2,P-1, "b");
    BigInt B = PowerMod(G, b, P);
    cout << "  then Bob computes  G^b mod P, and makes it public:" << endl
         << "  B = PowerMod(G, b, P) = " << B << " [PUBLIC KEY]" << endl;
    cout << endl;
    cout << "~~< Alice and Bob secret key >~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    cout << "  Now Alice and Bob compute G^(a*b) mod P = "
         << PowerMod(G,a*b,P) << endl
         << "  by doing respectively:" << endl
         << "    Alice: B^a mod P = PowerMod(B, a, P) = "
         << PowerMod(B, a, P) << endl 
         << "    Bob:   A^b mod P = PowerMod(A, b, P) = "
         << PowerMod(A, b, P) << endl
         << "  this is their SECRET key" << endl;
  }



  void ElGamal()
  {
    cout <<
      "====[  EL GAMAL algorithm  ]====\n"
      "  Alice and Bob pick together and publicly \n"
      "  a prime  P  and a generator  G  for (FF_P)*\n"
         << endl;
    BigInt P = InputPrime('P', "[PUBLIC]", 90);
    int G = PrimitiveRoot(P);
    cout << "  A generator for (FF_p)* is  G = "<< G << " [PUBLIC]" << endl;
    cout << endl;
    cout << "Alice picks  a [PRIVATE KEY] in 2,...," << P-1 << endl;
    BigInt a = InputBigIntInInterval(2,P-1, "a");
    BigInt A = PowerMod(G,a,P);
    cout << "  then computes  A = G^a mod P, and makes it public:" << endl
         << "    A = PowerMod(G, a, P) = " << A << " [PUBLIC KEY]" << endl;
    cout << endl;
    cout << "Bob picks  b [PRIVATE KEY] in 2,...," << P-1 << endl;
    BigInt b = InputBigIntInInterval(2,P-1, "b");
    BigInt B = PowerMod(G,b,P);
    cout << "  Bob computes  B = G^b mod P, and makes it public:" << endl
         << "    B = PowerMod(G, b, P) = " << B << " [PUBLIC KEY]" << endl;
    cout << endl;
    cout << "~~~~< Alice ciphers >~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    cout << "Alice wants to send Bob a message m (her phone number,..)" << endl
         << "in 0,...," << P-1 << endl;
    BigInt m = InputBigIntInInterval (0,P-1, "m");
    cout << endl;
    cout << "then picks a number k [PRIVATE] in 0,...," << P-1 << endl;
    BigInt k = InputBigIntInInterval (0,P-1, "k");
    BigInt C = LeastNNegRemainder(m*PowerMod(B, k, P), P);
    BigInt GK = PowerMod(G, k, P);
    cout << "  so she computes:" << endl
         << "    C = m * B^k mod P = " << C << endl
         << "    GK = G^k mod P = " << GK << endl
         << "  and sends (C, GK) [PUBLIC] =\n  ("
         << C << ", " << GK << ")" << endl;
    cout << endl;  
    cout << "~~~~< Bob deciphers >~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    cout << "  Now Bob, without knowing k, can compute  B^k as (GK)^b" << endl
         << "  and then deciphers  m = C * g^{-kb}" << endl
         << "  by doing C * InvMod(PowerMod(GK,b,P), P) mod P = "
         << LeastNNegRemainder( C * InvMod(PowerMod(GK,b,P), P), P) << endl;
  }



//----------------------------------------------------------------------
  void RSA()
  {
    cout <<
      "====[  RSA algorithm  ]====\n"
      "  Bob picks two distinct prime numbers p and q"
        << endl;
    BigInt p = InputPrime('p', "[PRIVATE]");
    cout <<" " << endl;
    BigInt q = InputPrime('q', "[PRIVATE]");
    const  BigInt N = p*q;
    const  BigInt f = p*q-p-q+1; 
    cout << "  and computes  N = p*q = " << N
         << " [PUBLIC]" << endl; 
    cout << "  and computes  f = phi(n) = (p-1)(q-1) = " << f << endl;
    cout << "~~~~[PUBLIC KEY]~~~~" << endl;
    cout << "Then he picks  E [PUBLIC] such that  gcd(E,f) = 1:" << endl;
    BigInt E = InputBigInt("E");
    while (gcd(E,f) != 1)
    {
      cout << "gcd(E,f) = " << gcd(E,f) << ".  Try again" << endl;
      E = InputBigInt("E");
    }
    cout << "gcd(E,f) = 1" << ".  Good!" << endl
      << "  PUBLIC KEY is  (N, E) = "
      <<  "\n  ( " << N
      << ",\n    "<< E << " )" << endl;
  
    cout << "~~~~[PRIVATE KEY]~~~~" << endl;
    cout
      << "  Finally he determines  d  such that  E*d is equivalent to 1 mod f"
      << endl;
    BigInt d = InvMod(E,f);
    cout
      << "  d = InvMod(E,f) =\n    = " << d << endl
      << "  PRIVATE KEY is (p, q, d) = "
      <<  "\n  ( " << p
      << ",\n    "<< q
      << ",\n    "<< d << " )" << endl;
    cout << endl;
    cout << "~~~~< Alice ciphers >~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    cout << "Alice wants to send Bob a message m (her phone number,..)" << endl
         << "in 0,...," << N-1 << endl;
    BigInt m = InputBigIntInInterval (0,N-1,"m");
    BigInt C = PowerMod(m, E, N);
    cout << "  so she ciphers using Bob's public key:  C = m^E mod N" << endl
         << "    C = PowerMod(m, E, N) = " << C << endl
         << "  and sends it [PUBLIC]" << endl;
    cout << endl;  
    cout << "~~~~< Bob deciphers >~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    cout << "  Now Bob deciphers using  d [PRIVATE]:  m = C^d mod N" << endl
         << "    PowerMod(C, d, N) = "
         << PowerMod(C, d, N) << endl;
  }


//----------------------------------------------------------------------
  void Rabin()
  {
    cout <<
      "====[  Rabin algorithm  ]====\n"
      "  Bob picks two distinct prime numbers p and q (same size)\n"
      "  such that p and q  mod 4 = 3" << endl;
    BigInt p = InputPrime('p', "");
    if (LeastNNegRemainder(p,4) != 3)
      cout << "<< p mod 4  should be  3:  picking next good prime >>" << endl;
    while (LeastNNegRemainder(p,4) != 3)      p = NextProbPrime(p);
    BigInt q = NextProbPrime(p + p/10); //InputPrime('q');
    while (LeastNNegRemainder(q,4) != 3)      q = NextProbPrime(q);
    cout << "  This is the first prime, congruent to 3 mod 4:" << endl
         << "    p = " << p << "[PRIVATE]" << endl
         << "  This is the another prime, similar size, congruent to 3 mod 4:" << endl
         << "    q = " << q << "[PRIVATE]" << endl;
    
    cout << "~~~~[PRIVATE KEY]~~~~" << endl;
    cout << "  PRIVATE KEY is (p,q) = "
         <<  "\n  ( " << p
         << ",\n    " << q << " )" << endl;

    cout << "~~~~[PUBLIC KEY]~~~~" << endl;
    const BigInt N = p*q ;
    cout << "  PUBLIC KEY is N = p*q = " << N << endl;

    //calcolo le chiavi

//     cout <<" " << endl;
//     cout << "Calcolo le chiavi ricordando che la chiave pubblica e' n = p*q e quella privata la coppia (p,q)" << endl<<endl;
//     cout << "La chiave privata e'(p,q) = " << "(" << p << "," << q << ")" << endl<<endl;
    //    cout << "La chiave pubblica e' n = " << n << endl<<endl;
    cout << endl;
    cout << "~~~~< Alice ciphers >~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    cout << "Alice wants to send Bob a message m (her phone number,..)" << endl
         << "in 0,...," << N-1 << endl;
//     cout << "inserisci il messaggio da cifrare, che deve essere un numero in {0...."<< n-1 << "}" <<  endl;
    BigInt m = InputBigIntInInterval(0,N-1,"m");
    BigInt C = PowerMod(m, 2, N);
    cout << "  so she ciphers using Bob's public key:  C = m^2 mod N" << endl
         << "    C = PowerMod(m, 2, N) = " << C << endl
         << "  and sends it [PUBLIC]" << endl;

//     cout << "La funzione di cifratura che un altro utente B deve usare per comunicare un messaggio m con A e' m^2 mod n" << endl<<endl;
//     BigInt c = PowerMod(m,2,n);
//     cout
//       << "Il messaggio cifrato e' c= " << c << endl<<endl;
  
    //decodifica
    cout << endl;
    cout << "~~~~< Bob deciphers >~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    BigInt mp = PowerMod(C, (p+1)/4, p);
    BigInt mq = PowerMod(C, (q+1)/4, q);
    cout <<
      "  The message may be deciphered by computing the square roots of C mod N\n"
      "  but Bob, knowing p,q, does it in polynomial-time as follows:\n"
    
//     cout << " A questo punto A deve risalire al messaggio m a partire dal testo cifrato c e ha proprio il problema di calcolare le radici di c mod n ma, dato che conosce p e q, riuscira' a calcolare le 4 radici di c mod n in tempo polinomiale. Per la decifrazione, l'utente A deve procedere come segue: " << endl<<endl;
      "  computes  mp = C^(p+1)/4   and  mq = C^(q+1)/4:" << endl;
    cout << "    mp = " << mp << endl;
    cout << "    mq = " << mq << endl;
    cout << "  +mp and -mp  are the square roots of  C  mod p" << endl;
    cout << "  +mq and -mq  are the square roots of  C  mod q" << endl;
//     cout << "calcola mp = c^(p+1)/4 e mq = c^(q+1)/4 e osserviamo che +/-mp e +/-mq sono le radici quadrate di c in Z_p e Z_q " << endl<<endl;
//     cout <<" " << endl;
    BigInt a = InvMod(p, q);
    BigInt b = InvMod(q, p);
    cout << "  Then computes a and b such that ap + bq = 1" << endl;
    cout << "    a = InvMod(p, q) = " << a << endl;
    cout << "    b = InvMod(q, p) = " << b << endl;
    cout << endl;
    cout << "  Finally he computes the four candidate solutions:" << endl;
    BigInt r1 = LeastNNegRemainder( a*p*mq + b*q*mp, N);
    BigInt r2 = LeastNNegRemainder( a*p*mq - b*q*mp, N);
    BigInt r3 = LeastNNegRemainder(-a*p*mq + b*q*mp, N);
    BigInt r4 = LeastNNegRemainder(-a*p*mq - b*q*mp, N);
    cout << "    r1 = ( a*p*mq + b*q*mp) mod N  = " << r1 << endl;
    cout << "    r2 = ( a*p*mq - b*q*mp) mod N  = " << r2 << endl;
    cout << "    r3 = (-a*p*mq + b*q*mp) mod N  = " << r3 << endl;
    cout << "    r4 = (-a*p*mq - b*q*mp) mod N  = " << r4 << endl;
    cout << "  One of which is  m";    
  }


//----------------------------------------------------------------------
//----------------------------------------------------------------------
//----------------------------------------------------------------------
  BigInt InputBigInt(const char* s)
  {
    BigInt n;
    cout << "number " << s << " > ";
    string line;
    while (true)
    {
      getline(cin, line);
      // This code converts from string to number safely.
      stringstream myStream(line);
      if (myStream >> n)  break;
      if (!line.empty())
      {
        cout << "Invalid number, please try again" << endl;
        cout << "number > ";
      }
    }
    return n;
  }
  

  BigInt InputPrime(char ch, const char* str, long MaxNumBits)
  {
    BigInt NumBits;
    cout << "How many bits for " << ch << "?" << endl;
    NumBits = InputBigInt("bits");
    while (NumBits < 2 || NumBits > MaxNumBits)  
    {
      cout << "Please, choose a number between 3 and " << MaxNumBits << endl; 
      NumBits = InputBigInt("bits");
    }
    if (NumBits > MaxNumBits)
      cout << "WARNING: might take a long time." << endl;
    BigInt min = power(2, NumBits-1);
    BigInt max = power(2, NumBits)-1;
  
    BigInt p = RandomBigInt(GlobalRandomSource(), min, max);
    while (!IsProbPrime(p))  p = RandomBigInt(GlobalRandomSource(),min, max);
  
    cout << "  A prime with " << NumBits << " bits is  "
         << ch << " = " << p << " " << str << endl;
    return p;
  }


  BigInt InputBigIntInInterval(long x, BigInt y, const char* prompt)
  {
    BigInt n;
    n = InputBigInt(prompt);
    while ((n<x) || (n>y))
    {
      cout << "This number is not in the given interval" << endl;
      cout << "Try again" << endl;
      n = InputBigInt(prompt);
    }
    return n;
  }
  

  void PrintMenu()
  {
    cout << "  'd' for  Diffie-Hellman " << endl;
    cout << "  'e' for  ElGamal " << endl;
    cout << "  'r' for  Rabin " << endl;
    cout << "  's' for  RSA " << endl;
    cout << "  'q' for quitting " << endl;
  }

  
  void program()
  {
    GlobalManager CoCoAFoundations;
    
    cout << ShortDescription << endl;
    cout << boolalpha; // so that bools print out as true/false
    SetVerbosityLevel(90);

    //\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/

    PrintMenu();
    char ch = 'a';
  
    while (ch != 'q')
    {
      cout << endl
           << "------------------------------------------------------------"
           << endl;
      cout << "Type 'm' for menu or another key" << endl;
      cout << "> ";
    
      cin >> ch ; 
      switch (ch)
      {
      case 'm': PrintMenu(); break;
      case 'd': DiffieHellman(); break;
      case 'e': ElGamal(); break;
      case 's': RSA(); break;
      case 'r': Rabin(); break;
      case 'q': break;
      
      default:
        cout << "La lettera inserita non appartiene a nessun caso, riprova" << endl;
        break;
      }
    }
  }

} // end of namespace CoCoA

//----------------------------------------------------------------------
// Use main() to handle any uncaught exceptions and warn the user about them.
int main()
{
  try
  {
    CoCoA::program();
    return 0;
  }
  // catch (const CoCoA::InterruptReceived& intr)
  // {
  //   cerr << endl
  //        << "------------------------------" << endl
  //        << ">>>  CoCoALib interrupted  <<<" << endl
  //        << "------------------------------" << endl
  //        << "-->>  " << intr << "  <<--" << endl;
  //   return 2;
  // }
  catch (const CoCoA::ErrorInfo& err)
  {
    cerr << "***ERROR***  UNCAUGHT CoCoA error";
    ANNOUNCE(cerr, err);
  }
  catch (const std::exception& exc)
  {
    cerr << "***ERROR***  UNCAUGHT std::exception: " << exc.what() << endl;
  }
  catch(...)
  {
    cerr << "***ERROR***  UNCAUGHT UNKNOWN EXCEPTION" << endl;
  }

  CoCoA::BuildInfo::PrintAll(cerr);
  return 1;
}


//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-NumTheory-crypto1.C,v 1.1 2018/05/16 14:05:42 bigatti Exp $
// $Log: ex-NumTheory-crypto1.C,v $
// Revision 1.1  2018/05/16 14:05:42  bigatti
// -- first import
//
