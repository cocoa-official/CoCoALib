// Copyright (c) 2023  John Abbott,  Anna M. Bigatti
// [incl major contributions from Nico Mexis]
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.


#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
    "This example program shows how to check for cyclotomicity. \n";

const string LongDescription =
    "This example program shows how to check for cyclotomicity.          \n"
    "The input is always expected to be a univariate polynomial over ZZ. \n"
    "Here, two cyclotomic and one non-cyclotomic polynomial are being    \n"
    "examined.";
//----------------------------------------------------------------------

namespace CoCoA
{

    void program()
    {
        GlobalManager CoCoAFoundations;
        SignalWatcher MonitorInterrupt(SIGINT); // you must also call CheckForInterrupt every so often

        cout << ShortDescription << endl;
        cout << boolalpha; // so that bools print out as true/false

        SparsePolyRing P = NewPolyRing(RingZZ(), symbols("x")); // ZZ[x];

        // Calculate the 8-th cyclotomic polynomial Phi(8)
        RingElem f = cyclotomic(8, indet(P,0));
        unsigned long fc = CyclotomicIndex(f);
        if (fc)
            cout << "The polynomial " << f << " is the " << fc << "-th cyclotomic." << endl;
        else
            cout << "The polynomial " << f << " is not cyclotomic." << endl;

        // This polynomial is the 20-th cyclotomic polynomial Phi(20)
        RingElem g = RingElem(P, "x^8-x^6+x^4-x^2+1");
        unsigned long gc = CyclotomicIndex(g);
        if (gc)
            cout << "The polynomial " << g << " is the " << gc << "-th cyclotomic." << endl;
        else
            cout << "The polynomial " << g << " is not cyclotomic." << endl;

        // This polynomial is not cyclotomic (Lehmer's polynomial)
        RingElem h = RingElem(P, "x^10+x^9-x^7-x^6-x^5-x^4-x^3+x+1");
        unsigned long hc = CyclotomicIndex(h);
        if (hc)
            cout << "The polynomial " << h << " is the " << hc << "-th cyclotomic." << endl;
        else
            cout << "The polynomial " << h << " is not cyclotomic." << endl;
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

    catch (const CoCoA::InterruptReceived &intr)
    {
        cerr << endl
             << "------------------------------" << endl
             << ">>>  CoCoALib interrupted  <<<" << endl
             << "------------------------------" << endl
             << "-->>  " << intr << "  <<--" << endl;
        return 2;
    }
    catch (const CoCoA::ErrorInfo &err)
    {
        cerr << "***ERROR***  UNCAUGHT CoCoA error";
        ANNOUNCE(cerr, err);
    }
    catch (const std::exception &exc)
    {
        cerr << "***ERROR***  UNCAUGHT std::exception: " << exc.what() << endl;
    }
    catch (...)
    {
        cerr << "***ERROR***  UNCAUGHT UNKNOWN EXCEPTION" << endl;
    }

    CoCoA::BuildInfo::PrintAll(cerr);
    return 1;
}
