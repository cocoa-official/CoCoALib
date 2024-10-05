//   Copyright (c)  2005  John Abbott,  Anna M Bigatti

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


#ifndef CoCoA_SocketStream_H
#define CoCoA_SocketStream_H

#include "CoCoA/config.H"

#if defined (_WIN32) || defined (_WIN64)
#include <winsock2.h>
#else
#include <netdb.h>
#endif
// using htons, gethostbyname, etc
#include <streambuf>
//using std::streambuf
#include <iostream>
//using std::iostream

// #include <netinet/in.h> // for macosx


namespace CoCoA
{

  class SocketStreambuf: public::std::streambuf
  {
  private:
    int mySocketFD;
    static const int ourMaxBacklog = 5; // for the call to listen
    static const std::size_t ourUngetSize = 4;                 // bytes
    static const std::size_t ourMaxChunkSize = (1<<16) - 100;  // bytes
    static const std::size_t ourInputBufferSize = ourMaxChunkSize + ourUngetSize; // bytes
    char myInputBuffer[ourInputBufferSize];
    size_t myPacketCounter;
    size_t mySentBytes;
  public:
    explicit SocketStreambuf(unsigned short PortNum);
    SocketStreambuf(const char hostname[], unsigned short PortNum);
    ~SocketStreambuf();
  protected:
    virtual int_type overflow(int_type c);
    virtual std::streamsize xsputn(const char* s, std::streamsize num);
    virtual int_type underflow();
  };


  /*-----------------------------------------------------------------*/
  /** \include SocketStream.txt  */
  /*-----------------------------------------------------------------*/
  class SocketStream: public std::iostream
  {
  protected:
    SocketStreambuf myStreambuf;
  public:
    explicit SocketStream(unsigned short PortNum);
    SocketStream(const char hostname[], unsigned short PortNum);
  };

} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/server/SocketStream.H,v 1.5 2024/03/25 08:04:26 abbott Exp $
// $Log: SocketStream.H,v $
// Revision 1.5  2024/03/25 08:04:26  abbott
// Summary: Modifications to work  on MinGW (& 32 bitters) (redmine 1661)
//
// Revision 1.4  2022/02/22 20:40:05  abbott
// Summary: Updated copyright message (redmine 1555)
//
// Revision 1.3  2007/10/30 17:14:06  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.2  2007/09/24 14:31:10  abbott
// Changed MaxChunkSize into a class level static constant ourMaxChunkSize.
//
// Revision 1.1.1.1  2007/03/09 15:16:12  abbott
// Imported files
//
// Revision 1.2  2007/02/12 18:39:58  bigatti
// -- added myPacketCounter and mySentBytes
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.1  2006/01/17 10:23:08  cocoa
// Updated DivMask; many consequential changes.
// A few other minor fixes.
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.1.1.1  2005/05/03 15:47:31  cocoa
// Imported files
//
// Revision 1.2  2005/04/19 14:06:04  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.5  2004/11/11 15:28:38  cocoa
// -- added \include *.txt for doxygen
//
// Revision 1.4  2004/07/16 10:54:04  cocoa
// -- added #include <netinet/in.h>  line needed for compiling on macosx
//
// Revision 1.3  2004/06/29 17:10:22  cocoa
// Partially tidied use of "protected" and "private" in various
// base classes.  Checking in at the end of the day -- it works,
// and I wouldn't want it to be lost next time point's disk
// misbehaves.
//

#endif
