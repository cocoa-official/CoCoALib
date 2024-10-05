
#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
    "This example program shows how to apply Graeffe's method. \n";

const string LongDescription =
    "This example program shows how to apply Graeffe's method.           \n"
    "When applying Graeffe's method to a polynomial in ZZ[x], the roots  \n"
    "of a polynomial get raised to some power. The graeffe function is   \n"
    "for squaring the roots, the graeffe3 function cubes them and the    \n"
    "GraeffeN function raises them to a custom power. The GraeffeSeq     \n"
    "class provides a way to iterate over all the powers up to a certain \n"
    "degree. In this example, the various functions are being applied    \n"
    "to a polynomial with roots that are easy to calculate.              \n";
//----------------------------------------------------------------------

namespace CoCoA
{

    void program()
    {
        GlobalManager CoCoAFoundations;
        SignalWatcher MonitorInterrupt(SIGINT); // you must also call CheckForInterrupt every so often

        cout << ShortDescription << endl;
        cout << boolalpha; // so that bools print out as true/false

        // Take a polynomial with nice roots
        SparsePolyRing P = NewPolyRing(RingZZ(), symbols("x")); // ZZ[x];
        RingElem f = RingElem(P, "(x-1)*(x-2)*(x+3)");
        cout << "f := " << f << endl;

        cout << "Using default Graeffe functions..." << endl;

        RingElem g2 = graeffe(f);
        cout << "Graeffe transformation: " << g2 << endl;

        RingElem g3 = graeffe3(f);
        cout << "Cubic Graeffe transformation: " << g3 << endl;

        RingElem g5 = GraeffeN(f, 5);
        cout << "5-th power Graeffe transformation: " << g5 << endl;

        cout << "Using GraeffeSeq iterator..." << endl;

        GraeffeSeq graeffeSeq(f, 5);
        for (int i = 1; i <= 5; i++)
        {
            ++graeffeSeq;
            cout << i << "-th power Graeffe transformation: " << *graeffeSeq << endl;
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
