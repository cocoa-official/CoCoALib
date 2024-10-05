
#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
    "This example program shows how to check a LRS for degeneracy. \n";

const string LongDescription =
    "This example program shows how to check a LRS for degeneracy.       \n"
    "In order to achieve this, the characteristic polynomial can be used \n"
    "used to determine whether it has two distinct roots whose quotient  \n"
    "is a root of unity. As a simple example, two LRS-degenerate         \n"
    "polynomials and one not-degenerate are examined.                    \n";
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

        // This polynomial is LRS-degenerate with degree 3
        RingElem f = RingElem(P, "x^11+x^7-x^6+x^3+x^2+x");
		unsigned long fn = IsLRSDegenerate(f);
        if (fn)
            cout << "The polynomial " << f << " is " << fn << "-LRS-degenerate." << endl;
        else
            cout << "The polynomial " << f << " is not LRS-degenerate." << endl;

        // This polynomial is LRS-degenerate with degree 2
        RingElem g = RingElem(P, "x^4+2*x^3+4*x^2+4*x+4");
		unsigned long gn = IsLRSDegenerate(g);
        if (gn)
            cout << "The polynomial " << g << " is " << gn << "-LRS-degenerate." << endl;
        else
            cout << "The polynomial " << g << " is not LRS-degenerate." << endl;

        // This polynomial is not LRS-degenerate
        RingElem h = RingElem(P, "x^4+x^3-3*x^2+x+1");
		unsigned long hn = IsLRSDegenerate(h);
        if (hn)
            cout << "The polynomial " << h << " is " << hn << "-LRS-degenerate." << endl;
        else
            cout << "The polynomial " << h << " is not LRS-degenerate." << endl;
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
