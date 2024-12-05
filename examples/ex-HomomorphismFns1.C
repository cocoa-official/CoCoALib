// Copyright (c) 2017  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This is a template file for example programs.  \n"
  "The program itself does nothing whatsoever.    \n";

const string LongDescription =
  "Make a copy of this file (called \"foo.C\", say) and put your code \n"
  "inside the procedure \"program\".                                  \n"
  "To compile your file in the examples directory just do this:       \n"
  "  make foo                                                         \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void PrintHomInfo(const RingHom& phi, ConstRefRingElem y)
  {
    cout << endl
         << "-------------------------------------------------------" << endl
         << "Homomorphism is phi = " << phi << endl
         << "IsInjective(phi) = " << IsInjective(phi) << endl
         << "IsSurjective(phi) = " << IsSurjective(phi) << endl;
    if (IsInImage(phi, y))
    {
      cout << y << " is in the image of phi" << endl
           << "preimage(phi,y) = " << preimage(phi,y) << endl
           << "ker(phi) = " << ker(phi) << endl;
    }
    else
    {
      cout << y << " is not in the image of phi" << endl;
    }
  }

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;
    cout << boolalpha; // so that bools print out as true/false

    //-------------------------------------------------------
    // Example of a RingHom from a poly ring to a poly ring
    ring P1 = NewPolyRing(RingQQ(), symbols("x,y,z"));
    RingHom phi1 = IdentityHom(P1);
    RingElem y1(P1, "x^2+y^2+z^2-1");
    PrintHomInfo(phi1, y1);

    //-------------------------------------------------------
    // Example of a RingHom from a poly ring to a quotient ring
    ring P2 = NewPolyRing(RingQQ(), symbols("a,b"));
    ideal I2 = ideal(RingElem(P2,"a^2-b^2"));
    QuotientRing Q2 = NewQuotientRing(P2, I2);
    vector<RingElem> images2;
    images2.push_back(RingElem(Q2, "a"));
    images2.push_back(RingElem(Q2, "b^2"));
    RingHom phi2 = PolyAlgebraHom(P2,Q2, images2);
    PrintHomInfo(phi2, RingElem(Q2,"a^3*b^2"));

    //-------------------------------------------------------
    // Example of RingHom from quot ring to poly ring (continued from above)
    vector<RingElem> images3;
    images3.push_back(RingElem(P2,"a+b"));
    images3.push_back(RingElem(P2,"a+b"));
    RingHom psi3 = PolyAlgebraHom(P2,P2, images3);
    RingHom phi3 = InducedHom(Q2, psi3);
    PrintHomInfo(phi3, RingElem(P2,"a-b"));
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
