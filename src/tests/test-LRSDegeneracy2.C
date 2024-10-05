
#include <fstream>
#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
    "Tests the current implementation of IsLRSDegenerate\n";
//----------------------------------------------------------------------

namespace CoCoA
{

    // Speed and correctness
    void CheckLRSDegeneracyOrder(const ring &Px, const string &input, long expected)
    {
        const RingElem f = RingElem(Px, input);

        const long calc = LRSDegeneracyOrder(f);
        CoCoA_ASSERT_ALWAYS(calc == expected);
    }


    void testDeg()
    {
        const SparsePolyRing Px = NewPolyRing(RingZZ(), symbols("x"));

        // SetVerbosityLevel(80);

        // File "LRSdegenerate.txt"

        // For k=3

        CheckLRSDegeneracyOrder(Px, "x^4+x^3+2*x^2-x+1", 3);     // Record[Exponents := [1], Factors := [x^4+x^3+2*x^2-x+1], RemainingFactor := 1]
        CheckLRSDegeneracyOrder(Px, "x^4+3*x^3+8*x^2+3*x+1", 3); // Record[Exponents := [1], Factors := [x^4+3*x^3+8*x^2+3*x+1], RemainingFactor := 1]

        // For k=4

        CheckLRSDegeneracyOrder(Px, "x^2+2*x+2", 4); // Record[Exponents := [1], Factors := [x^2+2*x+2], RemainingFactor := 1]

        // For k=5

        CheckLRSDegeneracyOrder(Px, "x^4+x^3+6*x^2-4*x+1", 5);   // Record[Exponents := [1], Factors := [x^4+x^3+6*x^2-4*x+1], RemainingFactor := 1]
        CheckLRSDegeneracyOrder(Px, "x^4+2*x^3+4*x^2+3*x+1", 5); // Record[Exponents := [1], Factors := [x^4+2*x^3+4*x^2+3*x+1], RemainingFactor := 1]

        // For k=6

        // this example is really just phi_3(x+1)--apparently the only example in deg=2
        CheckLRSDegeneracyOrder(Px, "x^2+3*x+3", 6); // Record[Exponents := [1], Factors := [x^2+3*x+3], RemainingFactor := 1]

        // Found none of deg =3
        CheckLRSDegeneracyOrder(Px, "x^4+3*x^3+2*x^2-3*x+1", 6); // Record[Exponents := [1], Factors := [x^4+3*x^3+2*x^2-3*x+1], RemainingFactor := 1]
        CheckLRSDegeneracyOrder(Px, "x^4+3*x^3+x^2-6*x+4", 6);   // Record[Exponents := [1], Factors := [x^4+3*x^3+x^2-6*x+4], RemainingFactor := 1]
        CheckLRSDegeneracyOrder(Px, "x^4+3*x^3+5*x^2+6*x+4", 6); // Record[Exponents := [1], Factors := [x^4+3*x^3+5*x^2+6*x+4], RemainingFactor := 1]
        CheckLRSDegeneracyOrder(Px, "x^4+3*x^3-9*x+9", 6);       // Record[Exponents := [1], Factors := [x^4+3*x^3-9*x+9], RemainingFactor := 1]

        // For k=7

        // Found none of deg < 6
        // phi(7) only example of height <= 3
        CheckLRSDegeneracyOrder(Px, "x^6+x^5+x^4+x^3+x^2+x+1", 7);

        // For k=8

        // phi_8(x+1)
        CheckLRSDegeneracyOrder(Px, "x^4+2*x^2+4*x+2", 8);       // Factors := [x^4+2*x^2+4*x+2]
        CheckLRSDegeneracyOrder(Px, "x^4+4*x^3+4*x^2+8", 8);     // derivable from poly above: (1/2)*rev(f(2*x))
        CheckLRSDegeneracyOrder(Px, "x^4+2*x^2+2", 2);           // of form f(x^2) so see also case k=4
        CheckLRSDegeneracyOrder(Px, "x^4+2*x^3-4*x-4", 2);       // Factors := [x^2+2*x+2, x^2-2]
        CheckLRSDegeneracyOrder(Px, "x^4+2*x^3+4*x^2+4*x+4", 2); // Factors := [x^2+2*x+2, x^2+2]
        CheckLRSDegeneracyOrder(Px, "x^4+4*x^3+6*x^2+4*x+2", 8); // irred, = phi_8(x+1)

        // For k=9

        // Found none of deg < 6
        // phi_9, also phi_18, also
        CheckLRSDegeneracyOrder(Px, "x^6+2*x^3+4", 3); // 4*phi_3((1/2)*x^3) Factors := [x^6+3*x^5+3*x^4-8*x^3+6*x^2-3*x+1]

        // For k=10

        CheckLRSDegeneracyOrder(Px, "x^4+5*x^3+5*x^2-5*x+5", 10); // Record[Exponents := [1], Factors := [x^4+5*x^3+5*x^2-5*x+5], RemainingFactor := 1]
        CheckLRSDegeneracyOrder(Px, "x^4+5*x+5", 10);             // Record[Exponents := [1], Factors := [x^4+5*x+5], RemainingFactor := 1]

        // For k=12

        CheckLRSDegeneracyOrder(Px, "x^4+3*x^2+3", 2);            // irred, see also k=6
        CheckLRSDegeneracyOrder(Px, "x^4+x^3+2*x^2+x+1", 2);      // Factors := [x^2+1, x^2+x+1]
        CheckLRSDegeneracyOrder(Px, "x^4+2*x^3+2*x^2+4*x+4", 3);  // irred
        CheckLRSDegeneracyOrder(Px, "x^4+2*x^3+5*x^2+4*x+1", 12); // irred, rev(phi_12(x+1))
        CheckLRSDegeneracyOrder(Px, "x^4+3*x^3-9*x-9", 2);        // Factors := [x^2+3*x+3, x^2-3]

        // For k=15

        // NOTE THIS EXAMPLE HAS DEGREE LOWER THAN deg(phi_15)
        CheckLRSDegeneracyOrder(Px, "x^6+2*x^5+3*x^4+3*x^3+3*x^2+2*x+1", 3); // Factors = phi_3*phi_5
        // (also phi_6*phi_10 works as it is just phi_3(-x)*phi5(-x))
        CheckLRSDegeneracyOrder(Px, "x^6-2*x^5+3*x^4-3*x^3+3*x^2-2*x+1", 3);
        // In deg 8, I found only phi_15 and multiples of (variants of) the poly above
        CheckLRSDegeneracyOrder(Px, "x^8-x^7+x^5-x^4+x^3-x+1", 3);

        // For k=18

        CheckLRSDegeneracyOrder(Px, "x^6+x^5+5*x^3+3*x^2-2*x+6", 3);           // Factors := [x^2-x+1, x^4+2*x^3+x^2+4*x+6]
        CheckLRSDegeneracyOrder(Px, "x^6+2*x^5-6*x^4+x^2-5*x-2", 3);           // Factors := [x^2-x+1, x^4+3*x^3-4*x^2-7*x-2]
        CheckLRSDegeneracyOrder(Px, "x^6+3*x^5+6*x^4+6*x^3+9*x^2+9*x+3", 18);  // irred, LRS 18
        CheckLRSDegeneracyOrder(Px, "x^6+7*x^5+5*x^4-2*x^3+8*x^2-9*x+1", 0);   // irred

        // File "CuriousPolynomials"

        CheckLRSDegeneracyOrder(Px, "x^4-2", 2);                    // And G is [[x+1, 4], [x^2+1, 4], [16, 1]]
        CheckLRSDegeneracyOrder(Px, "x^4+2*x^2+4", 2);              // And G is [[x^2-x+1, 2], [x^2+x+1, 2], [x+1, 4], [256, 1]]
        CheckLRSDegeneracyOrder(Px, "x^4+3*x^2+9", 2);              // And G is [[x^2+x+1, 2], [x^2-x+1, 2], [x+1, 4], [6561, 1]]
        CheckLRSDegeneracyOrder(Px, "x^4+x^3+6*x^2-4*x+1", 5);      // And G is [[x^4+11*x^3+46*x^2-4*x+1, 1], [x^4-4*x^3+46*x^2+11*x+1, 1], [x^4+x^3+x^2+x+1, 1]]
        CheckLRSDegeneracyOrder(Px, "x^4+2*x^3+4*x^2+3*x+1", 5);    // And G is [[x^4+x^3+x^2+x+1, 1], [x^4+x^3+6*x^2-4*x+1, 1], [x^4-4*x^3+6*x^2+x+1, 1]]
        CheckLRSDegeneracyOrder(Px, "x^4+3*x^3+19*x^2+7*x+1", 5);   // And G is [[x^4+11*x^3+321*x^2-29*x+1, 1], [x^4-29*x^3+321*x^2+11*x+1, 1], [x^4+x^3+x^2+x+1, 1]]
        CheckLRSDegeneracyOrder(Px, "x^4+4*x^3+19*x^2-12*x+9", 3);  // And G is [[9*x^4-66*x^3+475*x^2-66*x+9, 1], [x^2+x+1, 2], [3*x^2+22*x+3, 2], [81, 1]]
        CheckLRSDegeneracyOrder(Px, "x^4+6*x^3+16*x^2+16*x+16", 5); // And G is [[x^4+x^3+6*x^2-4*x+1, 1], [x^4-4*x^3+6*x^2+x+1, 1], [x^4+x^3+x^2+x+1, 1], [65536, 1]]

        // Testing up to height 4
        // Doing prefix [1, 0]
        CheckLRSDegeneracyOrder(Px, "x^4+2*x^2-4*x+2", 8); // and K=8
        CheckLRSDegeneracyOrder(Px, "x^4+2*x^2+4*x+2", 8); // and K=8
        // Doing prefix [1, 1]
        CheckLRSDegeneracyOrder(Px, "x^4+x^3-x^2+2*x+4", 3);   // and K=3
        CheckLRSDegeneracyOrder(Px, "x^4+x^3+x^2+x+1", 5);     // and K=5  (obviously!)
        CheckLRSDegeneracyOrder(Px, "x^4+x^3+2*x^2-x+1", 3);   // and K=3
        CheckLRSDegeneracyOrder(Px, "x^4+x^3+6*x^2-4*x+1", 5); // and K=5
        // Doing prefix [1, 2]
        CheckLRSDegeneracyOrder(Px, "x^4+2*x^3+2*x^2-4*x+4", 4); // and K=4
        CheckLRSDegeneracyOrder(Px, "x^4+2*x^3+2*x^2-2*x+1", 4); // and K=4
        CheckLRSDegeneracyOrder(Px, "x^4+2*x^3+2*x^2+4*x+4", 3); // and K=3,4  (no other 3,4 polys with height <= 99, except F(2*x))
        CheckLRSDegeneracyOrder(Px, "x^4+2*x^3+4*x^2+3*x+1", 5); // and K=5
        // Doing prefix [1, 3]
        CheckLRSDegeneracyOrder(Px, "x^4+3*x^3+2*x^2-3*x+1", 6); // and K=6
        // Doing prefix [1, 4]
        CheckLRSDegeneracyOrder(Px, "x^4+4*x^3+6*x^2+4*x+2", 8); // and K=8

        CheckLRSDegeneracyOrder(Px, "x^4+5*x+5", 10);             // and K=10
        CheckLRSDegeneracyOrder(Px, "x^4+5*x^3+5*x^2-5*x+5", 10); // and K=10

        // NONE in deg=5 (except of the form x^5-n)

        // DEG=6

        CheckLRSDegeneracyOrder(Px, "x^6-x^4-2*x^3+x^2+x+1", 3); // And G is [[x^12-3*x^11+2*x^10-3*x^9+15*x^8-24*x^7+25*x^6-24*x^5+15*x^4-3*x^3+2*x^2-3*x+1, 1], [x^6+3*x^5+7*x^4+9*x^3+7*x^2+3*x+1, 2], [x^2+x+1, 3]]
        CheckLRSDegeneracyOrder(Px, "x^6-x^4+2*x^3+x^2-x+1", 3); // And G is [[x^12-3*x^11+2*x^10-3*x^9+15*x^8-24*x^7+25*x^6-24*x^5+15*x^4-3*x^3+2*x^2-3*x+1, 1], [x^6+3*x^5+7*x^4+9*x^3+7*x^2+3*x+1, 2], [x^2+x+1, 3]]
        CheckLRSDegeneracyOrder(Px, "x^6+x^3-1", 3);             // And G is [[x^6+3*x^3+1, 3], [x^2+x+1, 6]]
        // >>>Lots of the form  x^6+A*x^3+B<<<

        CheckLRSDegeneracyOrder(Px, "x^6+x^4+2*x^3+x^2+x+1", 3);    // And G is [[x^12-3*x^11+4*x^10-5*x^9+5*x^8+2*x^7-7*x^6+2*x^5+5*x^4-5*x^3+4*x^2-3*x+1, 1], [x^6+3*x^5+5*x^4+5*x^3+5*x^2+3*x+1, 2], [x^2+x+1, 3]]
        CheckLRSDegeneracyOrder(Px, "x^6+2*x^4+2*x^3+2*x^2+1", 0); // And G is [[x^6+4*x^5+8*x^4+6*x^3+8*x^2+4*x+1, 1], [x^6+x^5-4*x^4-9*x^3-4*x^2+x+1, 2], [x^6+2*x^4+2*x^3+2*x^2+1, 2]]

        // These polys have non trivial GCD with Resultant(F(x*y), F(y), y)
        CheckLRSDegeneracyOrder(Px, "x^6+x^5+2*x^4+3*x^3+2*x^2+x+1", 0);   // And G is [[x^6+3*x^5+2*x^4-x^3+2*x^2+3*x+1, 1], [x^3-x^2-x-1, 2], [x^3+x^2+x-1, 2], [x^6+x^5+2*x^4+3*x^3+2*x^2+x+1, 2]]
        CheckLRSDegeneracyOrder(Px, "x^6+x^5+3*x^4-x^3+3*x^2+x+1", 0);     // And G is [[x^6+5*x^5+17*x^4+17*x^3+17*x^2+5*x+1, 1], [x^6-x^5+3*x^4+x^3+3*x^2-x+1, 2], [x^6+x^5-9*x^4-19*x^3-9*x^2+x+1, 2]]
        CheckLRSDegeneracyOrder(Px, "x^6+3*x^5+x^4-3*x^3+x^2+3*x+1", 0);   // And G is [[x^6-7*x^5+21*x^4-23*x^3+21*x^2-7*x+1, 1], [x^3-3*x^2-2*x-1, 2], [x^3+2*x^2+3*x-1, 2], [x^6+3*x^5+x^4-3*x^3+x^2+3*x+1, 2]]
        CheckLRSDegeneracyOrder(Px, "x^6+3*x^5+2*x^4-x^3+2*x^2+3*x+1", 0); // And G is [[x^6-5*x^5+14*x^4-9*x^3+14*x^2-5*x+1, 1], [x^3+x^2+3*x-1, 2], [x^3-3*x^2-x-1, 2], [x^6+3*x^5+2*x^4-x^3+2*x^2+3*x+1, 2]]
        CheckLRSDegeneracyOrder(Px, "x^6+3*x^5+3*x^4+x^3+3*x^2+3*x+1", 0); // And G is [[x^6-3*x^5+9*x^4+x^3+9*x^2-3*x+1, 1], [x^3-3*x^2-1, 2], [x^3+3*x-1, 2], [x^6+3*x^5+3*x^4+x^3+3*x^2+3*x+1, 2]]

        CheckLRSDegeneracyOrder(Px, "x^6+2*x^5+4*x^4+x^3+2*x^2-3*x+1 ", 7);   // with K=7 (only example with height<=4, except Phi_7)
        CheckLRSDegeneracyOrder(Px, "x^6+2*x^5+4*x^4+x^3+2*x^2-3*x+1   ", 7); // with K=7 (searched up to height 9)
        CheckLRSDegeneracyOrder(Px, "x^6+2*x^5+4*x^4+8*x^3+9*x^2+4*x+1 ", 7); // with K=7 (searched up to height 9)

        // Exhaustive search up to height 4
        CheckLRSDegeneracyOrder(Px, "x^6-2*x^4-4*x^3+4*x^2+4*x+4", 3);       // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6-2*x^4-2*x^3+4*x^2+2*x+1", 3);       // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6-2*x^4+2*x^3+4*x^2-2*x+1", 3);       // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6-2*x^4+4*x^3+4*x^2-4*x+4", 3);       // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6-x^4-2*x^3+x^2+x+1", 3);             // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6-x^4+x^2-3*x+3", 6);                 // and K=6
        CheckLRSDegeneracyOrder(Px, "x^6-x^4+x^2+3*x+3", 6);                 // and K=6
        CheckLRSDegeneracyOrder(Px, "x^6-x^4+2*x^3+x^2-x+1", 3);             // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6-4*x^3-4", 3);                       // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6-4*x^3-3", 3);                       // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6-4*x^3-2", 3);                       // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6-4*x^3+1", 3);                       // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6-4*x^3+2", 3);                       // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6-3*x^3-3", 3);                       // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6-3*x^3-2", 3);                       // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6-3*x^3-1", 3);                       // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6-3*x^3+1", 3);                       // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6-3*x^3+3", 3);                       // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6-3*x^3+4", 3);                       // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6-2*x^3-4", 3);                       // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6-2*x^3-2", 3);                       // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6-2*x^3-1", 3);                       // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6-2*x^3+2", 3);                       // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6-2*x^3+3", 3);                       // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6-2*x^3+4", 3);                       // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6-2*x^3+x^2+2*x+2", 4);               // and K=4
        CheckLRSDegeneracyOrder(Px, "x^6-2*x^3+4*x^2-4*x+2", 4);             // and K=4
        CheckLRSDegeneracyOrder(Px, "x^6-2*x^3+4*x^2+4*x+2", 4);             // and K=4
        CheckLRSDegeneracyOrder(Px, "x^6-x^3-4", 3);                         // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6-x^3-3", 3);                         // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6-x^3-1", 3);                         // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6-x^3+1", 3);                         // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6-x^3+2", 3);                         // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6-x^3+3", 3);                         // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6-x^3+4", 3);                         // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6+x^3-4", 3);                         // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6+x^3-3", 3);                         // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6+x^3-1", 3);                         // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6+x^3+1", 3);                         // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6+x^3+2", 3);                         // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6+x^3+3", 3);                         // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6+x^3+4", 3);                         // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6+2*x^3-4", 3);                       // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6+2*x^3-2", 3);                       // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6+2*x^3-1", 3);                       // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6+2*x^3+2", 3);                       // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6+2*x^3+3", 3);                       // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6+2*x^3+4", 3);                       // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6+2*x^3+x^2-2*x+2", 4);               // and K=4
        CheckLRSDegeneracyOrder(Px, "x^6+2*x^3+4*x^2-4*x+2", 4);             // and K=4
        CheckLRSDegeneracyOrder(Px, "x^6+2*x^3+4*x^2+4*x+2", 4);             // and K=4
        CheckLRSDegeneracyOrder(Px, "x^6+3*x^3-3", 3);                       // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6+3*x^3-2", 3);                       // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6+3*x^3-1", 3);                       // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6+3*x^3+1", 3);                       // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6+3*x^3+3", 3);                       // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6+3*x^3+4", 3);                       // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6+4*x^3-4", 3);                       // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6+4*x^3-3", 3);                       // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6+4*x^3-2", 3);                       // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6+4*x^3+1", 3);                       // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6+4*x^3+2", 3);                       // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6+x^4-4*x^3+x^2-2*x+4", 3);           // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6+x^4-2*x^3+x^2-x+1", 3);             // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6+x^4+x^2-3*x+3", 6);                 // and K=6
        CheckLRSDegeneracyOrder(Px, "x^6+x^4+x^2+3*x+3", 6);                 // and K=6
        CheckLRSDegeneracyOrder(Px, "x^6+x^4+2*x^3+x^2+x+1", 3);             // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6+x^4+4*x^3+x^2+2*x+4", 3);           // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6+2*x^4-4*x^3+4*x^2-4*x+4", 3);       // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6+2*x^4+4*x^3+4*x^2+4*x+4", 3);       // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6+x^5-x^4+3*x^2+2*x+1", 3);           // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6+x^5-3*x^3-x^2+2*x+4", 3);           // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6+x^5+3*x^3+2*x^2-x+1", 3);           // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6+x^5+2*x^4-3*x^3-x+1", 3);           // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6+x^5+2*x^4+3*x^3+3*x^2+2*x+4", 3);   // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6+x^5+3*x^4-4*x^3+3*x^2-2*x+1", 3);   // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6+2*x^5+2*x^4-4*x^3+x^2-2*x+2", 4);   // and K=4
        CheckLRSDegeneracyOrder(Px, "x^6+2*x^5+2*x^4-2*x^3+4*x^2+4*x+2", 4); // and K=4
        CheckLRSDegeneracyOrder(Px, "x^6+2*x^5+2*x^4+4*x+4", 3);             // and K=3
        CheckLRSDegeneracyOrder(Px, "x^6+2*x^5+2*x^4+2*x^3+4*x^2+4*x+2", 4); // and K=4
        CheckLRSDegeneracyOrder(Px, "x^6+2*x^5+2*x^4+4*x^3+x^2-2*x+2", 4);   // and K=4
        CheckLRSDegeneracyOrder(Px, "x^6+3*x^5+2*x^4-3*x^3+4*x^2-3*x+3", 6); // and K=6
        CheckLRSDegeneracyOrder(Px, "x^6+3*x^5+4*x^4+3*x^3-2*x^2-3*x+3", 6); // and K=6

        // DEG=8

        CheckLRSDegeneracyOrder(Px, "x^8+x^7-x^5-x^4-x^3+x+1", 3);         // and K=3  i.e. two roots have ratio a 3-th root of unity
        CheckLRSDegeneracyOrder(Px, "x^8-x^6-2*x^5+x^3+x+1", 3);           // and K=3
        CheckLRSDegeneracyOrder(Px, "x^8-x^6+2*x^5-x^3-x+1", 3);           // and K=3
        CheckLRSDegeneracyOrder(Px, "x^8-2*x^5-x^4+x^2+x+1", 3);           // and K=3
        CheckLRSDegeneracyOrder(Px, "x^8-2*x^5+x^4+x^2-x+1", 3);           // and K=3
        CheckLRSDegeneracyOrder(Px, "x^8+2*x^5-x^4+x^2-x+1", 3);           // and K=3
        CheckLRSDegeneracyOrder(Px, "x^8+2*x^5+x^4+x^2+x+1", 3);           // and K=3
        CheckLRSDegeneracyOrder(Px, "x^8+x^7-x^5-x^4-x^3+x+1", 3);         // and K=3 (already found scanning through height=1)
        CheckLRSDegeneracyOrder(Px, "x^8+x^7+x^6+2*x^5-2*x^3+x^2-x+1", 3); // and K=3

        // Palindromic examples:
        CheckLRSDegeneracyOrder(Px, "x^8+x^7+2*x^6-3*x^5-x^4-3*x^3+2*x^2+x+1", 3);       // and K=3
        CheckLRSDegeneracyOrder(Px, "x^8+2*x^7+2*x^6-2*x^5-2*x^4-2*x^3+2*x^2+2*x+1", 4); // and K=4
        CheckLRSDegeneracyOrder(Px, "x^8+3*x^7+2*x^6-3*x^5+3*x^4-3*x^3+2*x^2+3*x+1", 6); // and K=6
        CheckLRSDegeneracyOrder(Px, "x^8+3*x^7+3*x^6+2*x^4+3*x^2+3*x+1", 6);             // and K=6
        CheckLRSDegeneracyOrder(Px, "x^8+x^7+3*x^6-4*x^5+2*x^4-4*x^3+3*x^2+x+1", 3);     // and K=3
        CheckLRSDegeneracyOrder(Px, "x^8+2*x^7+2*x^6-4*x^5-x^4-4*x^3+2*x^2+2*x+1", 4);   // and K=4
        CheckLRSDegeneracyOrder(Px, "x^8+2*x^7+3*x^6-2*x^5-4*x^4-2*x^3+3*x^2+2*x+1", 3); // and K=3
        CheckLRSDegeneracyOrder(Px, "x^8+3*x^7+4*x^6+3*x^5+3*x^4+3*x^3+4*x^2+3*x+1", 6); // and K=6
        CheckLRSDegeneracyOrder(Px, "x^8+4*x^7+8*x^6+8*x^5+7*x^4+8*x^3+8*x^2+4*x+1", 4); // and K=4

        CheckLRSDegeneracyOrder(Px, "x^8+x^7-x^5-x^4-x^3+x+1", 3);                  // and K=3,5
        CheckLRSDegeneracyOrder(Px, "x^8+2*x^7+2*x^5+9*x^4+8*x^3+5*x^2+3*x+1 ", 3); // and K=3,5 (no others like this of height <= 9)

        //-----------------------------------------------------------------------------

        CheckLRSDegeneracyOrder(Px, "x^10+x^6-x^5+x^2+x+1", 3); // and K=3
        CheckLRSDegeneracyOrder(Px, "x^10+x^6+x^5+x^2-x+1", 3); // and K=3

        // Own stupid tests

        CheckLRSDegeneracyOrder(Px, "(x-3)*(x+3)", 2);
        CheckLRSDegeneracyOrder(Px, "x^2-1", 2);
        CheckLRSDegeneracyOrder(Px, "x^101-1", 101);
        CheckLRSDegeneracyOrder(Px, "x^101+1", 101);

        // Very huge degenerate polynomial (n=102307); disable Mut test for this
        /* string hugePolyString;
        ifstream hugePolyIn;
        hugePolyIn.open("degenerate-large");
        getline(hugePolyIn, hugePolyString);
        hugePolyIn.close();
        CheckLRSDegeneracyOrder(Px, hugePolyString, true); */

        // Tests from Cipu's paper

        CheckLRSDegeneracyOrder(Px, "x^12-6*x^11+23*x^10-73*x^9+191*x^8-405*x^7+766*x^6-1164*x^5+1368*x^4-1539*x^3+1863*x^2-1701*x+729", 13);
        CheckLRSDegeneracyOrder(Px, "-82*x^24+572*x^23-2328*x^22+7589*x^21-20626*x^20+45830*x^19-88495*x^18+143301*x^17-181934*x^16+203609*x^15-247594*x^14+266625*x^13-161067*x^12+28519*x^11-40816*x^10+100473*x^9-78511*x^8+88658*x^7-191918*x^6+258252*x^5-237528*x^4+211410*x^3-168723*x^2+75573*x-10206", 2);
        CheckLRSDegeneracyOrder(Px, "-7*x^24-14*x^23-21*x^22+x^21+58*x^20+87*x^19+32*x^18-168*x^17-252*x^16-345*x^15-18*x^14-27*x^13-27*x^12+18*x^11+27*x^10+20*x^9-32*x^8-48*x^7-26*x^6+76*x^5+114*x^4+59*x^3-186*x^2-279*x-372", 3);
        CheckLRSDegeneracyOrder(Px, "x^4+x^3-3*x^2+x+1", 0);
        CheckLRSDegeneracyOrder(Px, "x^10+x^9-x^7-x^6-x^5-x^4-x^3+x+1", 0);
        CheckLRSDegeneracyOrder(Px, "-294*x^12+567*x^11-378*x^10+1478*x^9+1924*x^8+2369*x^6+2601*x^5-1062*x^4+1290*x^3+7764*x^2+6615*x+3235", 0);
        CheckLRSDegeneracyOrder(Px, "-6370*x^24-990*x^23+2576*x^22-5779*x^21+260*x^20+1757*x^19+90*x^18-4755*x^17-7705*x^16+4115*x^15-10010*x^14-8757*x^13+4704*x^12+5576*x^11+216*x^10+5257*x^9+9215*x^8+2695*x^6-4510*x^5+3185*x^4-5330*x^3+1", 0);
    }



} // end of namespace CoCoA


//----------------------------------------------------------------------
// Use main() to handle any uncaught exceptions and warn the user about them.
int main()
{
    CoCoA::GlobalManager CoCoAFoundations;
    CoCoA::SignalWatcher MonitorInterrupt(SIGINT); // you must also call CheckForInterrupt every so often

    cout << boolalpha; // so that bools print out as true/false

    try
    {
        CoCoA::testDeg(); // test IsDegenerate function
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
