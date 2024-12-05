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

#endif
