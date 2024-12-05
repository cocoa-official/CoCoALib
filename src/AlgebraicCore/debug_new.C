//   Copyright (c)  2005,2016  John Abbott and Anna M. Bigatti

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


#include "CoCoA/debug_new.H"

#include <cstddef>
using std::size_t;
#include <cstdlib>
using std::malloc;
#include <iostream>
using std::ostream;
using std::cerr;
using std::endl;
#include <new>
using std::bad_alloc;


namespace debug_new
{

  constexpr size_t MAX_ALLOC = 128*1024*1024; // Max request in a single alloc; 128Mbytes.
  constexpr int MARGIN_MARKER = 1234567890;
  constexpr int MARGIN_MARKER_AFTER_FREE = 1122334455;

  // Global variables; don't hope for thread-safety with this code!
  size_t TOTAL_MALLOCKED = 0;
  size_t TOTAL_FREED = 0;
  size_t CURRENT_MEM = 0;
  size_t PEAK_MEM = 0;

  size_t NUM_NEW = 0;
  size_t NEW_WATCH_POINT = 0;
  size_t NUM_DELETE = 0;
  size_t DELETE_WATCH_POINT = 0;

  bool InhibitMsgs = true;

  PrintTrace::PrintTrace(bool activate)
  {
    PreviousState = InhibitMsgs;
    // make sure cerr has been used before trying to print log mesgs on it.
    InhibitMsgs = true;
    if (activate)
      cerr << "[debug_new] ACTIVATING LOGGING OF NEW/DELETE REQUESTS (ctor for " << this << ")" << endl;
    else
      cerr << "[debug_new] SUSPENDING LOGGING OF NEW/DELETE REQUESTS (ctor for " << this << ")" << endl;

    InhibitMsgs = !activate;
  }

  PrintTrace::~PrintTrace()
  {
    if (PreviousState != InhibitMsgs)
    {
      if (InhibitMsgs)
	cerr << "[debug_new] ACTIVATING LOGGING OF NEW/DELETE REQUESTS (dtor for " << this << ")" << endl;
      else
	cerr << "[debug_new] SUSPENDING LOGGING OF NEW/DELETE REQUESTS (dtor for " << this << ")" << endl;
    }
    InhibitMsgs = PreviousState;
  }


  void InterceptNew(unsigned long nth)
  {
    if (NUM_NEW >= nth)
      cerr << "[debug_new] ERROR: InterceptNew(" << nth << ") called too late" << endl;
    NEW_WATCH_POINT = nth;
  }


  void InterceptDelete(unsigned long nth)
  {
    if (NUM_DELETE >= nth)
      cerr << "[debug_new] ERROR: InterceptDelete(" << nth << ") called too late" << endl;
    DELETE_WATCH_POINT = nth;
  }

  /* MARGIN is the number of "ints" allocated immediately before and after  */
  /* each block.  These margins are checked for integrity when the space is */
  /* freed, or when CHECKMARGINS is explicitly called.  The margins are     */
  /* invisible to the caller.  Bigger margins give better safety but waste  */
  /* more space.  The value should be at least 2.                           */
  const size_t MARGIN = 8;


  void intercepted()
  {
    cerr << "[debug_new] INTERCEPTED" << endl;
  }


  static void msg_delete_zero(ostream& out)
  {
    if (InhibitMsgs) return;
    out << "[debug_new] ERROR: DELETE: ZERO POINTER (seq=" << NUM_DELETE << ")" << endl;
    if (NUM_DELETE == DELETE_WATCH_POINT) intercepted();
  }

  static void msg_double_delete(ostream& out, void* addr)
  {
    if (InhibitMsgs) return;
    out << "[debug_new] ERROR: DELETE: probable DOUBLE DELETE on pointer " << addr
        << " (seq=" << NUM_DELETE << ")" << endl;
    if (NUM_DELETE == DELETE_WATCH_POINT) intercepted();
  }

  static void msg_delete_bad_size(ostream& out, size_t sz)
  {
    if (InhibitMsgs) return;
    out << "[debug_new] ERROR: DELETE: bad block size : "
        << sz << " > " << MAX_ALLOC << " = max alloc allowed"
        << " (seq=" << NUM_DELETE << ")" << endl;
    if (NUM_DELETE == DELETE_WATCH_POINT) intercepted();
  }

  static void msg_trouble(ostream& out, void* addr, void* block_lo, void* block_hi, size_t sz)
  {
    if (InhibitMsgs) return;
    out << "[debug_new] ERROR: CHECKMARGINS: trouble at " << addr << " for block "
        << block_lo << "-" << block_hi << " (size=" << sz << ")" << endl;
  }

  static void msg_alloc(ostream& out, const int* block, size_t intsize, size_t sz)
  {
    if (InhibitMsgs) return;
    out << "[debug_new] ALLOC " << block+MARGIN << " (" << sz << " bytes) \tseq=" << NUM_NEW
	<< " (extd block is " << block << "-" << block+(intsize+2*MARGIN)-1
	<< ") \tA=" << TOTAL_MALLOCKED/1024
	<< "k F=" << TOTAL_FREED/1024
	<< "k C=" << CURRENT_MEM/1024
	<< "k MAX=" << PEAK_MEM/1024 << "k" << endl;
    if (NUM_NEW == NEW_WATCH_POINT) intercepted();
  }

  static void msg_freed(ostream& out, const int* block, size_t intsize, size_t sz, void* ptr)
  {
    if (InhibitMsgs) return;
    out << "[debug_new] FREED " << ptr << " (" << sz << " bytes) \tseq=" << NUM_DELETE
	<< " (extd block is " << block << "-" << block+(intsize+2*MARGIN)-1
	<< ") \tA=" << TOTAL_MALLOCKED/1024
	<< "k F=" << TOTAL_FREED/1024
	<< "k C=" << CURRENT_MEM/1024
	<< "k MAX=" << PEAK_MEM/1024 << "k" << endl;
    if (NUM_DELETE == DELETE_WATCH_POINT) intercepted();
  }

  static void msg_freed_corrupted(ostream& out, const int* block, size_t sz, void* ptr)
  {
    if (InhibitMsgs) return;
    out << "[debug_new] FREED " << ptr << " (" << sz << "? bytes) \tseq=" << NUM_DELETE
	<< " (CORRUPTED extd block " << block << "-\?\?)\t<<<stats left unaltered>>>" << endl;
    if (NUM_DELETE == DELETE_WATCH_POINT) intercepted();
  }


  static int CHECKMARGINS(void *ptr)
  {
    int *block;
    size_t sz;
    size_t intsize, i;
    int num_errs = 0;

    block = (int*)((char*)ptr-MARGIN*sizeof(int));
    sz = block[0];
    if (sz > MAX_ALLOC)
    {
      cerr << "[debug_new] ERROR: CHECKMARGINS: suspiciously large block (sz=" << sz << ") at " << ptr << endl;
      return 1;
    }
    for (i=1; i < MARGIN; ++i)
      if (block[i] != MARGIN_MARKER)
      {
        msg_trouble(cerr, &block[i], ptr, ((char*)ptr+sz), sz);
        ++num_errs;
      }
    intsize = 1+sz/sizeof(int);
    for (i=MARGIN+intsize; i < intsize+2*MARGIN; ++i)
      if (block[i] != MARGIN_MARKER)
      {
        msg_trouble(cerr, &block[i], ptr, ((char*)ptr+sz), sz);
        ++num_errs;
      }
    return num_errs;
  }


} // end of namespace debug_new

using namespace debug_new;

void* operator new(size_t sz)
#if __cplusplus < 201100L
  // C++98 requires this exception specification; C++11 rejects it
  throw (std::bad_alloc)
#endif
{
  if (sz > MAX_ALLOC) throw std::bad_alloc();
  ++NUM_NEW;

  size_t intsize = 1+sz/sizeof(int);
  int* block = (int*)malloc((intsize+2*MARGIN)*sizeof(int));
  block[0] = sz;
  size_t i;
  for (i=1; i < MARGIN; ++i)  block[i] = MARGIN_MARKER;
  for (i=MARGIN; i < MARGIN+intsize; ++i) block[i] = -999999999;
  for (i=MARGIN+intsize; i < intsize+2*MARGIN; ++i) block[i] = MARGIN_MARKER;
  TOTAL_MALLOCKED += sz;
  CURRENT_MEM += sz;
  if (CURRENT_MEM > PEAK_MEM) PEAK_MEM = CURRENT_MEM;
  msg_alloc(cerr, block, intsize, sz);
  return &block[MARGIN];
}


void operator delete(void *ptr) noexcept
{
  ++NUM_DELETE;

  if (ptr == nullptr) { msg_delete_zero(cerr); return; }
  int* block = (int*)((char*)ptr-MARGIN*sizeof(int));
  size_t sz = block[0];
  if (sz == MARGIN_MARKER_AFTER_FREE) { msg_double_delete(cerr, ptr); return; }
  if (sz > MAX_ALLOC) { msg_delete_bad_size(cerr, sz); return; }

  if (CHECKMARGINS(ptr) > 0)
  {
    msg_freed_corrupted(cerr, block, sz, ptr);
    return;
  }
  TOTAL_FREED += sz;
  CURRENT_MEM -= sz;
  size_t intsize = 1+sz/sizeof(int);
  msg_freed(cerr, block, intsize, sz, ptr);
  for (size_t i=0; i < intsize+2*MARGIN; ++i) block[i] = MARGIN_MARKER_AFTER_FREE;
  free(block);
}

// Compiler says I must define this: BUG???? I hope that checking sz == req_sz suffices
void operator delete(void *ptr, std::size_t req_sz) noexcept
{
  ++NUM_DELETE;

  if (ptr == nullptr) { msg_delete_zero(cerr); return; }
  int* block = (int*)((char*)ptr-MARGIN*sizeof(int));
  size_t sz = block[0];
  if (sz == MARGIN_MARKER_AFTER_FREE) { msg_double_delete(cerr, ptr); return; }
  if (sz != req_sz) { msg_delete_bad_size(cerr, sz); return; }

  if (CHECKMARGINS(ptr) > 0)
  {
    msg_freed_corrupted(cerr, block, sz, ptr);
    return;
  }
  TOTAL_FREED += sz;
  CURRENT_MEM -= sz;
  size_t intsize = 1+sz/sizeof(int);
  msg_freed(cerr, block, intsize, sz, ptr);
  for (size_t i=0; i < intsize+2*MARGIN; ++i) block[i] = MARGIN_MARKER_AFTER_FREE;
  free(block);
}
