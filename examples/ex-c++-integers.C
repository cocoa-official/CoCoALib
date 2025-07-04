// Copyright (c) 2022  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"
using namespace CoCoA;

#include <limits>
using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This is an example showing the C++ integer type \"long\", which\n"
  "we recommend just for indexing; otherwise use \"BigInt\".  We  \n"
  "also give an ugly example about how to detect integer overflow.\n";

const string LongDescription =
  "This is an example showing the C++ integer type \"long\".  There  \n"
  "are other types with different bit-widths; all can be \"signed\"  \n"
  "or \"unsigned\" (just for non-negative values).  Computations     \n"
  "are fast, but may cause overflow -- something you should avoid! \n"
  "We give an example showing how awkward it can be to avoid or    \n"
  "anticipate overflow.  There are no such problems with BigInt!   \n";

//----------------------------------------------------------------------

namespace CoCoA
{

  // GENERAL COMMENT:
  //   C++ offers many different integral types.  When working
  //   with CoCoALib we recommend using just "long"
  //   (this is actually shorthand for "signed long int")
  //   The type "char" is good for characters in strings.
  //   C++ does "automatic promotion" of integral values
  //   [this is both good and bad].


  void program()
  {
    GlobalManager CoCoAFoundations;
    cout << ShortDescription << endl;

    // Integer constants such as 3 or -99 are known as "literals"
    // (you may see this name in a compiler message).
    // For literals outside the range -9999 to 9999 it is a good
    // idea to append the suffix "L" -- a lower case "l" is also
    // allowed but looks confusingly similar to the digit "1".
    // IMPORTANT: do not write any leading zeroes: 033 is not 33.

    const long SAFE_NEGATIVE = -2147483647L; // OK on all binary computers
    const long SAFE_POSITIVE =  2147483647L; // OK on all binary computers
    cout << "All values between " << SAFE_NEGATIVE
         << " and " << SAFE_POSITIVE << " are OK" << endl;

    // Integers outside the range above may not work on some
    // computers.  The actual permitted range can be found like this:
    const long MOST_NEGATIVE = numeric_limits<long>::min();
    const long MOST_POSITIVE = numeric_limits<long>::max();
    if (MOST_NEGATIVE < SAFE_NEGATIVE  ||  MOST_POSITIVE > SAFE_POSITIVE)
      cout << "Actual range of representable values is wider." << endl << endl;

    // It is likely that -MOST_POSITIVE works as expected,
    // but that -MOST_NEGATIVE does not work as expected!


    // What are longs useful for?
    // (*) size of a container (e.g. vector or string)
    // (*) index into a container
    // (*) loops indexing over elements in a container
    // (*) representing values from a limited range (e.g. years)

    // Arithmetic with longs:
    // (*) integer division (e.g.  2/3  gives  0)
    // (*) be careful that OVERFLOW does not occur
    //     [*YOU* must check -- the computer does not]
    // (*) if in doubt use BigInt (slower, but safer)
    
    // Silly example:  compute 1+2+3+...+n
    long n = 1000000L;  // Try 5000000000L as well...
    long sum = 0;
    for (long i=1; i <= n; ++i)
      sum += i;
    cout << "Sum of integers 1 up to " << n << "  is " << sum << endl;
    // sum will overflow if the result is greater than MOST_POSITIVE


    // ------------------------
    // *** !!! ADVANCED !!! ***
    // ------------------------
    //
    // The rest of this example is significantly trickier.
    // We want to compute the sum directly by using the "theoretical" formula:
    //   sum(1..n)  =  (1/2)*n*(n+1)
    // We ignore overflow here; that will be discussed further below.
    long bad_formula    = (1/2)*n*(n+1); // Always zero: 1/2 uses integer division!
    long better_formula = (n*(n+1))/2;   // n*(n+1)/2  is same by C++ assoc rules.
    long best_formula   = (n%2 == 0) ? ((n+1)*(n/2)) : (n*((n+1)/2));

    cout << "bad_formula gives "       << bad_formula << endl;
    cout << "better_formula gives "    << better_formula << endl;
    cout << "best_formula gives "      << best_formula   << endl;

    // ***BEFORE READING ON***:  Why did I call them as I did?

    // bad_formula:  while CoCoALib has the ability to compute
    //     with rational values (via BigRat), a "rational constant",
    //     such as 1/2, is handled by C++ as INTEGER DIVISION.
    //     So 1/2 == 0.  Be careful!!
    
    // better_formula:  using brackets we instruct the compiler
    //      to multiply first, then divide; the problem is that
    //      the multiplication may overflow even if the answer
    //      is small enough to fit into a long.
    //      NOTE: C++ associativity rules mean that n*(n+1)/2 is interpreted
    //      the same as (n*(n+1))/2; I prefer to use brackets to be explicit.

    // best_formula:  "best" but not "most readable"!
    //      depending on whether n is even or odd we tell
    //      the compiler to divide first (but knowing that
    //      the division will be exact), and then multiply.
    //      This will overflow iff sum in the loop overflows.

    // !!CAUTION!!  C++ associativity rules dictate how operators bind,
    //              but order of evaluation of the arguments is
    //              "in parallel".  Be careful if args have side-effects!
    //              Maybe your clever expression works...
    //              ...today,
    //              ...on your current computer,
    //              ...with your current compiler.
    //              But it may give different results with a different
    //              computer/compiler.   Make sure you write unambiguous
    //              portable code!

    
    // -----------------------------------------------------------------
    // WARNING: HARD, UGLY PART BELOW (how to anticipate/avoid overflow)
    // -----------------------------------------------------------------
    // Take a deep breath...
    // [AIM: convince you to use BigInt for general integer computations]
    
    // How could we check if overflow will occur?

    // We must check beforehand (to avoid "undefined behaviour"), but the
    // checks make the code both slower and harder to read  :-(
    // Below is my code to evaluate the above formula while checking
    // for overflow (and printing out a warning if overflow would occur).

    // For simplicity we ASSUME that n >= 0.
    long safest = 0;  // will hold correct value iff no overflow occurs
    if (n == 0)
      safest = 0; // handle the case n=0
    else if (n == MOST_POSITIVE)
        cerr << "!!! Overflow would occur !!!" << endl;
    // Henceforth we know that n+1 cannot cause overflow
    else if (n%2 == 0)
    {
      // Case n is even & n > 0
      const long halfn = n/2;  // division is exact, cannot overflow, halfn > 0.
      if (MOST_POSITIVE/halfn >= n+1)  // division OK because halfn non-zero.
        safest = halfn*(n+1); // no overflow in (n+1); no overflow in product.
      else
        cerr << "!!! Overflow would occur !!!" << endl;
    }
    else
    {
      // Case n is odd, n > 0
      const long halfn1 = (n+1)/2; // (n+1) cannot overflow, division is exact, halfn1 > 0.
      if (MOST_POSITIVE/n >= halfn1) // division OK because n non-zero
        safest = halfn1*n; // no overflow in product.
      else
        cerr << "!!! Overflow would occur !!!" << endl;
    }
    // Correct answer is in "safest", or else we printed an overflow message.
    cout << "If no overflow message then correct answer is " << safest << endl;
    // Phew! (no, I did not get the code right on my first attempt)

    // NOTE: in this specific instance (n*(n+1))/2 is monotonic increasing, so
    // we could just have checked whether n is greater than the largest value
    // which works (BUT but how do you determine this value?).  The aim here
    // was to exhibit a general technique.


    // ----------------------
    // ----- CONCLUSION -----
    // ----------------------
    //   In general it is best to use (signed or unsigned) longs just for indexing,
    //   and to use BigInts for general integer arithmetic (they're slow but safe).
    //   Avoid mixing "signed" and "unsigned" values in a single formula.

    // Tricky exercise:
    //     write 4 fns for safe addition, subtraction, multiplication & division
    //     of longs which avoid overflow (and print a message if overflow would occur).

    // For masochists:
    //  - look up integral types (char, short, int, long, long long)
    //  - look up integral promotions (e.g. on cppreference)
    //  - look up about unsigned integral types (overflow is silent but "harmless")
  }

} // end of namespace CoCoA

// IGNORE THE STUFF BELOW (at least for now)

//----------------------------------------------------------------------
// Use main() to handle any uncaught exceptions and warn the user about them.
int main()
{
  try
  {
    CoCoA::program();
    return 0;
  }
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
