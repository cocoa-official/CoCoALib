#ifndef CoCoA_random_H
#define CoCoA_random_H

//   Copyright (c)  2010-2012  John Abbott, and Anna M. Bigatti

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

#include "CoCoA/BigInt.H"
#include "CoCoA/MachineInt.H"

#include "gmp.h"
#include <bitset>
//using std::bitset;
#include <iosfwd>
//using std::ostream;

namespace CoCoA
{

  class RandomSource
  {
  public:
    explicit RandomSource(const MachineInt& seed);
    explicit RandomSource(const BigInt& seed = BigInt(/*0*/));
    RandomSource(const RandomSource& copy);
    RandomSource& operator=(const RandomSource& rhs);
    ~RandomSource();
    friend bool RandomBool(RandomSource& RndSrc) noexcept;
    //    friend long RandomLong(RandomSource& RndSrc);
    friend long RandomLong(RandomSource& RndSrc, const MachineInt& lwb, const MachineInt& upb);
    friend BigInt RandomBigInt(RandomSource& RndSrc, const MachineInt& lwb, const MachineInt& upb);
    friend BigInt RandomBigInt(RandomSource& RndSrc, const MachineInt& lwb, const BigInt& upb);
    friend BigInt RandomBigInt(RandomSource& RndSrc, const BigInt& lwb, const MachineInt& upb);
    friend BigInt RandomBigInt(RandomSource& RndSrc, const BigInt& lwb, const BigInt& upb);
//    friend BigInt RandomBigIntmod(RandomSource& RndSrc, const BigInt& M);  // range is 0..M-1
    friend void reseed(RandomSource& RndSrc, const MachineInt& seed);
    friend void reseed(RandomSource& RndSrc, const BigInt& seed);
    friend std::ostream& operator<<(std::ostream& out, const RandomSource& RndSrc);
  private:  // data members
    BigInt mySeed;
    gmp_randstate_t myState;
//     long myBoolCounter;
//     long myLongCounter;
//     long myBigIntCounter;
  };


  void reseed(RandomSource& RndSrc, const MachineInt& seed);
  void reseed(RandomSource& RndSrc, const BigInt& seed);

  inline bool RandomBool(RandomSource& RndSrc) noexcept
  {
    return gmp_urandomb_ui(RndSrc.myState, 1) == 1;
  }

  // CONVENIENCE FUNCTIONS using GlobalRandomSource()
  bool   RandomBool() noexcept;
  bool   RandomBiasedBool(double P); // returns true with probability P
///  bool   RandomBiasedBool(unsigned long N, unsigned long D); // returns true with probability N/D
  long   RandomLong(const MachineInt& lwb, const MachineInt& upb);
  BigInt RandomBigInt(const MachineInt& lwb, const MachineInt& upb);
  BigInt RandomBigInt(const MachineInt& lwb, const BigInt& upb);
  BigInt RandomBigInt(const BigInt& lwb, const MachineInt& upb);
  BigInt RandomBigInt(const BigInt& lwb, const BigInt& upb);

  //////////////////////////////////////////////////////////////////


  class RandomSeqLong
  {
  public:
    RandomSeqLong(const MachineInt& lwb, const MachineInt& upb, const MachineInt& seed = 0);
    RandomSeqLong(const RandomSeqLong& copy);
    RandomSeqLong& operator=(const RandomSeqLong& rhs);
    ~RandomSeqLong();
    long operator*() const noexcept;
    RandomSeqLong& operator++() noexcept;
    RandomSeqLong operator++(int);
    long myIndex() const noexcept;
  private:
    void myGenValue() noexcept; // fills myValue with a newly generated random value.
  private: // data members
    long myLwb;
    long myUpb;
    unsigned long myRange;
    long myValue;
    long myCounter;
    unsigned long mySeed;
    gmp_randstate_t myState;
    friend std::ostream& operator<<(std::ostream& out, const RandomSeqLong& RndLong);
  };

  long NextValue(RandomSeqLong& RndLong) noexcept;
  std::ostream& operator<<(std::ostream& out, const RandomSeqLong& RndLong);

  //---------------------------------------------------------------------------
  // Inline functions.

  inline long RandomSeqLong::operator*() const noexcept
  {
    return myValue;
  }


  inline RandomSeqLong& RandomSeqLong::operator++() noexcept
  {
    ++myCounter;
    myGenValue();
    return *this;
  }


  inline long NextValue(RandomSeqLong& RndLong) noexcept
  {
    return *++RndLong; // advance then get value
  }


  //////////////////////////////////////////////////////////////////


  class RandomSeqBool
  {
  public:
    explicit RandomSeqBool(const MachineInt& seed = 0);
    RandomSeqBool(const RandomSeqBool& copy);
    RandomSeqBool& operator=(const RandomSeqBool& rhs);
    ~RandomSeqBool();
    bool operator*() const noexcept;
    RandomSeqBool& operator++() noexcept;
    RandomSeqBool operator++(int);
    long myIndex() const noexcept;
  private:
    static const int ourBufferBits = std::numeric_limits<unsigned long>::digits;
    void myFillBuffer(); // fills myBuffer and resets myBitIndex to 0
  private: // data members
    int myBitIndex;
    std::bitset<ourBufferBits> myBuffer;
    long myCounter;
    unsigned long mySeed;
    gmp_randstate_t myState;

    friend std::ostream& operator<<(std::ostream& out, const RandomSeqBool& RndBool);
  };

  bool NextValue(RandomSeqBool& RndBool) noexcept;
  bool NextBiasedBool(RandomSeqBool& RndBool, double P); // returns true with probability P
///  bool   NextBiasedBool(RandomSeqBool& RndBool, unsigned long N, unsigned long D); // returns true with probability N/D
  std::ostream& operator<<(std::ostream& out, const RandomSeqBool& RndBool);

  //---------------------------------------------------------------------------
  // Inline functions.

  inline bool RandomSeqBool::operator*() const noexcept
  {
    return myBuffer[myBitIndex];
  }


  inline RandomSeqBool& RandomSeqBool::operator++() noexcept
  {
    ++myCounter;
    if (++myBitIndex >= ourBufferBits)
      myFillBuffer();
    return *this;
  }


  inline bool NextValue(RandomSeqBool& RndBool) noexcept
  {
    return *++RndBool; // advance then get value
  }


  //////////////////////////////////////////////////////////////////


  class RandomSeqBigInt
  {
  public:
    RandomSeqBigInt(const MachineInt& lwb, const MachineInt& upb, const MachineInt& seed = 0);
    RandomSeqBigInt(const MachineInt& lwb, const BigInt& upb, const MachineInt& seed = 0);
    RandomSeqBigInt(const BigInt& lwb, const MachineInt& upb, const MachineInt& seed = 0);
    RandomSeqBigInt(const BigInt& lwb, const BigInt& upb, const MachineInt& seed = 0);
    RandomSeqBigInt(const RandomSeqBigInt& copy);
    RandomSeqBigInt& operator=(const RandomSeqBigInt& rhs);
    ~RandomSeqBigInt();
    const BigInt& operator*() const noexcept;
    RandomSeqBigInt& operator++();
    RandomSeqBigInt operator++(int);
    long myIndex() const noexcept;
  private:
    void myGenValue(); // fills myValue with a newly generated random value.
  private: // data members
    BigInt myLwb;   // morally const, changed only by assignment
    BigInt myUpb;   // morally const, changed only by assignment
    BigInt myRange; // morally const, changed only by assignment
    BigInt myValue;
    long myCounter;
    unsigned long mySeed;
    gmp_randstate_t myState;
    friend std::ostream& operator<<(std::ostream& out, const RandomSeqBigInt& RndBigInt);
  };

  const BigInt& NextValue(RandomSeqBigInt& RndBigInt);
  std::ostream& operator<<(std::ostream& out, const RandomSeqBigInt& RndBigInt);

  //---------------------------------------------------------------------------
  // Inline functions.

  inline const BigInt& RandomSeqBigInt::operator*() const noexcept
  {
    return myValue;
  }


  inline RandomSeqBigInt& RandomSeqBigInt::operator++()
  {
    ++myCounter;
    myGenValue();
    return *this;
  }


  inline const BigInt& NextValue(RandomSeqBigInt& RndBigInt)
  {
    return *++RndBigInt; // advance then get value
  }


} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/random.H,v 1.12 2022/02/18 14:11:52 abbott Exp $
// $Log: random.H,v $
// Revision 1.12  2022/02/18 14:11:52  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.11  2021/02/10 19:32:01  abbott
// Summary: Added noexcept (redmine 1572)
//
// Revision 1.10  2021/01/07 14:48:33  abbott
// Summary: Corrected copyright
//
// Revision 1.9  2019/02/28 11:42:26  abbott
// Summary: Added reseed
//
// Revision 1.8  2018/04/20 18:51:25  abbott
// Summary: Changed ctors for BigInt/BigRat from string or from MPZ/MPQ
//
// Revision 1.7  2013/05/14 14:21:59  abbott
// Improved a comment.
//
// Revision 1.6  2013/02/19 18:50:46  abbott
// Added RandomBiasedBool for (small) rational probabilities; but it is
// commented out in the header.
//
// Revision 1.5  2013/02/15 17:44:38  abbott
// Added RandomBiasedBool; changed name prob -->  NextBiasedBool (swapped args too).
//
// Revision 1.4  2012/12/05 11:03:17  abbott
// Renamed RandomLongStream   --> RandomSeqLong
//         RandomBoolStream   --> RandomSeqBool
//         RandomBigIntStream --> RandomSeqBigInt
//
// Revision 1.3  2012/12/04 20:05:23  abbott
// New unified header and source for all random generators.
//
//

#endif
