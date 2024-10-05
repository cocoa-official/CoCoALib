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


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/debug_new.C,v 1.11 2024/03/08 19:59:59 abbott Exp $
// $Log: debug_new.C,v $
// Revision 1.11  2024/03/08 19:59:59  abbott
// Summary: Used constexpr (redmine 1511)
//
// Revision 1.10  2022/02/18 14:12:02  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.9  2021/04/26 13:59:52  abbott
// Summary: Added new operator delete (because compiler said I had to, not sure it is correct)
//
// Revision 1.8  2021/02/10 19:40:00  abbott
// Summary: Added noexcept (sometimes instead of throw()) -- redmine 1572
//
// Revision 1.7  2021/01/07 15:23:38  abbott
// Summary: Corrected copyright
//
// Revision 1.6  2019/03/19 11:07:08  abbott
// Summary: Replaced 0 by nullptr where appropriate
//
// Revision 1.5  2016/04/26 14:33:45  abbott
// Summary: Added CPP if command to make code work with C++11
//
// Revision 1.4  2016/03/30 16:08:16  abbott
// Summary: Added a comment about C++11 incompatibility
//
// Revision 1.3  2011/03/16 15:33:19  abbott
// Changed arg type "unsigned int" into "unsigned long" for the Intercept* fns.
//
// Revision 1.2  2007/10/30 17:14:06  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.2  2006/10/19 13:56:24  cocoa
// Added #include<new> whenever needed (i.e. for placement new).
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.3  2006/03/27 12:21:25  cocoa
// Minor silly changes to reduce number of complaints from some compiler or other.
//
// Revision 1.2  2006/03/12 21:28:33  cocoa
// Major check in after many changes
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.1.1.1  2005/05/03 15:47:31  cocoa
// Imported files
//
// Revision 1.3  2005/04/19 14:06:03  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.2  2005/04/01 16:18:20  cocoa
// Friday check-in.  Fixed a bug in the ctor for GeneralQuotientRingImpl.
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.3  2004/11/03 17:58:56  cocoa
// -- improved log messages
// -- minor tidying
//
// Revision 1.2  2004/01/28 15:34:19  cocoa
// Now allow larger blocks before CheckMargins gets suspicious
// (previously 10^6, now 10^7).
//
// Revision 1.1.1.1  2003/09/24 12:55:43  cocoa
// Imported files
//
// Revision 1.6  2002/03/08 18:29:24  abbott
// Added new functions debug_new::InterceptNew and debug_new::InterceptDelete
// (plus accompanying paraphernalia).
//
// Revision 1.5  2002/03/08 15:22:19  abbott
// Added some printing to the ctor of a verbose_obj: this ensures that cerr
// is used at least once before logging is turned on (otherwise buffering
// inside the cerr object causes an infinite loop).  I have also put a message
// in the dtor just for symmetry.  I doubt this fix is really portable, but
// it works for gcc-3 currently (which is what counts right now!).
//
// Revision 1.4  2002/02/18 16:28:09  abbott
// Output now sent to cerr rather than cout -- I think cout may be buffered,
// and this can cause an infinite loop when writing "ALLOC' messages.
//
// Revision 1.3  2002/01/31 11:31:04  abbott
// Added the debug_new::verbose_obj mechanism for turning on and off printing.
//
// Revision 1.2  2002/01/11 13:37:40  abbott
// Changed "Ord" into "seq" to be more compatible with the documentation
// (and the behaviour of MemPool in verbose mode).
//
// Revision 1.1  2001/11/07 19:02:37  abbott
// Initial revision
//

