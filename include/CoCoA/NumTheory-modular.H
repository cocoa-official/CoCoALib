#ifndef CoCoA_NumTheory_modular_H
#define CoCoA_NumTheory_modular_H

//   Copyright (c)  2009,2012,2020  John Abbott,  Anna M. Bigatti

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


namespace CoCoA
{
  
  class MachineInt; // fwd decl -- defined in MachineInt.H
  class BigInt; // fwd decl -- defined in BigInt.H

  // Frankly, 8 separate cases seems a bit excessive!
  long   PowerMod(const MachineInt& base, const MachineInt& exp, const MachineInt& modulus);  // NOT noexcept - arg check
  long   PowerMod(const MachineInt& base, const BigInt& exp,     const MachineInt& modulus);  // ditto
  long   PowerMod(const BigInt& base,     const MachineInt& exp, const MachineInt& modulus);  // ditto
  long   PowerMod(const BigInt& base,     const BigInt& exp,     const MachineInt& modulus);  // ditto
  BigInt PowerMod(const MachineInt& base, const MachineInt& exp, const BigInt& modulus);
  BigInt PowerMod(const BigInt& base,     const MachineInt& exp, const BigInt& modulus);
  BigInt PowerMod(const MachineInt& base, const BigInt& exp,     const BigInt& modulus);
  BigInt PowerMod(const BigInt& base,     const BigInt& exp,     const BigInt& modulus);


  enum class InvModErrorAction { ThrowOnError, RtnZeroOnError };
  constexpr InvModErrorAction ThrowOnError = InvModErrorAction::ThrowOnError;     // safely inject
  constexpr InvModErrorAction RtnZeroOnError = InvModErrorAction::RtnZeroOnError; // safely inject
  unsigned long InvModNoArgCheck(unsigned long r, unsigned long m, const InvModErrorAction ErrorAction = ThrowOnError);  // READ THE DOC!!!  ASSUMES 0 <= r < m
  long   InvMod(const MachineInt& r, const MachineInt& m, const InvModErrorAction ErrorAction = ThrowOnError);
  long   InvMod(const BigInt& r,     const MachineInt& m, const InvModErrorAction ErrorAction = ThrowOnError);
  BigInt InvMod(const MachineInt& r, const BigInt& m, const InvModErrorAction ErrorAction = ThrowOnError);
  BigInt InvMod(const BigInt& r,     const BigInt& m, const InvModErrorAction ErrorAction = ThrowOnError);



  long   MultiplicativeOrderMod(const MachineInt& residue, const MachineInt& modulus);
  long   MultiplicativeOrderMod(const BigInt& residue,     const MachineInt& modulus);
  BigInt MultiplicativeOrderMod(const MachineInt& residue, const BigInt& modulus);
  BigInt MultiplicativeOrderMod(const BigInt& residue,     const BigInt& modulus);     
  long PrimitiveRoot(const MachineInt& p); // NOT noexcept -- arg check
  long PrimitiveRoot(const BigInt& P);


  //-------------------------------------------------------
  // Quadratic residue testing:
  // Kronecker is generalization of Jacobi which is generalization of Legendre Symbol
  long KroneckerSymbol(const MachineInt& residue, const MachineInt& modulus);  // NOT noexcept -- arg check
  long KroneckerSymbol(const BigInt& residue,     const MachineInt& modulus); 
  long KroneckerSymbol(const MachineInt& residue, const BigInt& modulus);
  long KroneckerSymbol(const BigInt& residue,     const BigInt& modulus);

} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/NumTheory-modular.H,v 1.4 2022/02/18 14:11:46 abbott Exp $
// $Log: NumTheory-modular.H,v $
// Revision 1.4  2022/02/18 14:11:46  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.3  2021/03/03 22:09:32  abbott
// Summary: New enum class (redmine 894)
//
// Revision 1.2  2021/02/10 19:32:00  abbott
// Summary: Added noexcept (redmine 1572)
//
// Revision 1.1  2020/01/26 14:14:31  abbott
// Summary: Finished splitting NumTheory into smaller pieces (redmine 1161)
//
//
//

#endif
