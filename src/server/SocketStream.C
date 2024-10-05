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


#include "SocketStream.H"
#include "CoCoA/error.H"
#include "GlobalIO.H"

#include <unistd.h>
// using fork
#include <cstdio>
// to bring the preprocessor macro EOF into "scope"
#include <cstring>
using std::memmove;
#include <sys/types.h>
// using fork
#if defined (_WIN32) || defined (_WIN64)
#include <winsock2.h>
#include <ws2tcpip.h>
#else
#include <sys/socket.h>
// using socket, bind, listen, connect, etc.
//#include <netinet/in.h> // useful???
#include <netinet/in.h>
// using struct sockaddr_in
#include <netdb.h>
// using htons, gethostbyname, etc
#endif
#include <signal.h>
// using signal and SIG_IGN

//#include <streambuf>
using std::streambuf;
//#include <iostream>
using std::iostream;
#include <algorithm>
using std::min;

#if defined (_WIN32) || defined (_WIN64)
#define SIGCHLD 17
#endif


namespace CoCoA
{

  // Static constants -- for some reason the compiler needs this.
  const std::size_t SocketStreambuf::ourMaxChunkSize;


  SocketStreambuf::SocketStreambuf(unsigned short PortNum)
  {
    myPacketCounter=0;
    mySentBytes=0;

    const int ListenFD = socket(PF_INET, SOCK_STREAM, 0);
    if (ListenFD < 0)
      CoCoA_THROW_ERROR("listening socket creation failed", "SocketStreambuf ctor for server end");

    struct sockaddr_in ListenAddress;
    ListenAddress.sin_family = AF_INET;
    ListenAddress.sin_port = htons(PortNum);
    ListenAddress.sin_addr.s_addr = (int)htonl(INADDR_ANY);

    {
      // Experimental patch -- should allow immediate reconnection if server is killed.
#if defined (_WIN32) || defined (_WIN64)
      const char dummy=1;
#else
      int dummy=1;
#endif
      setsockopt(ListenFD, SOL_SOCKET, SO_REUSEADDR, &dummy, sizeof(int));
    }

    const int BindOK = bind(ListenFD, (struct sockaddr*)&ListenAddress, sizeof(ListenAddress));
    if (BindOK < 0)
      CoCoA_THROW_ERROR("socket bind failed (port probably already in use)", "SocketStreambuf ctor for server end");
	
    listen(ListenFD, ourMaxBacklog);

    // Must catch child termination signals, o/w the children become zombies
    signal(SIGCHLD, SIG_IGN);

    struct sockaddr_in ClientAddress;
    socklen_t AddrSize = sizeof(ClientAddress);
    while (true)
    {
      mySocketFD = accept(ListenFD, (struct sockaddr*)&ClientAddress, &AddrSize);
      if (mySocketFD < 0) 
        CoCoA_THROW_ERROR("accept failed (how can this happen?)", "SocketStreambuf ctor for server end");
      if (fork() == 0) break; // ctor succeeds only in the child
      close(mySocketFD);
    }

    setg(myInputBuffer+ourUngetSize, myInputBuffer+ourUngetSize, myInputBuffer+ourUngetSize);
  }


  SocketStreambuf::SocketStreambuf(const char hostname[], unsigned short PortNum)
  {
    myPacketCounter=0;
    mySentBytes=0;

    /* IPv6 wants PF_INET6 */
    mySocketFD = socket(PF_INET, SOCK_STREAM, 0);
    if (mySocketFD < 0)
      CoCoA_THROW_ERROR("socket creation failed", "SocketStreambuf ctor for client end");

    struct sockaddr_in ServerAddress;
    ServerAddress.sin_family = AF_INET; /* IPv6 AF_INET6 */
    ServerAddress.sin_port = htons(PortNum);
	
    /* We need to translate the hostname to an IP address 
     * gethostbyname() queries /etc/hosts, DNS, NIS, ... in the correct
     * order. IPv6 needs gethostinfo(), not available on older systems
     */
    struct hostent *HostInfo;
    HostInfo = gethostbyname(hostname);
    if (HostInfo == NULL) /* DNS resolver needs its own error checking */
    {
      CoCoA_THROW_ERROR("HostInfo failed", "SocketStreambuf ctor for client end");
// 		fprintf(stderr, "Error: cocoa_connect: %s\n", hstrerror(h_errno));
// 		switch(h_errno)
// 		{
// 			case HOST_NOT_FOUND: return PARAM_ERR;
// 			case NO_ADDRESS:
// 			case NO_RECOVERY: return DNS_ERR;
// 			case TRY_AGAIN: return DNS_WARN;
// 			default: return UNKN_ERR;
// 		}
    }
	
    ServerAddress.sin_addr = *(struct in_addr *)HostInfo->h_addr;
	
    /* Connect using the socket mySocketFD to the server specified by ServerAddress struct */
    const int ConnectOK = connect(mySocketFD, (struct sockaddr*)&ServerAddress, sizeof(ServerAddress));
    if (ConnectOK < 0)
      CoCoA_THROW_ERROR("connect failed (how can this happen?)", "SocketStreambuf ctor for client end");

    setg(myInputBuffer+ourUngetSize, myInputBuffer+ourUngetSize, myInputBuffer+ourUngetSize);
  }

  SocketStreambuf::~SocketStreambuf()
  {
    // FIXME: add proper logging - mabshoff 2007-03-06
    //GlobalLogput() << "[Log] System=SocketStream Step=TerminateStream" << std::endl;
    //GlobalLogput() << "[Log] System=SocketStream myPacketCounter=" << myPacketCounter << " mySentBytes=" 
    //               << mySentBytes << " BytesPerPacket=" << (((double)mySentBytes)/((double)myPacketCounter)) << std::endl;
  }

  std::streamsize SocketStreambuf::xsputn(const char *OutgoingMesg, std::streamsize MesgLen)
  {
    ++myPacketCounter;
    mySentBytes += MesgLen;

#ifdef CoCoASTREAMDEBUG
    GlobalLogput() << "[Log] System=SocketStream MesgLen=" << MesgLen << " Message=|";
    for ( std::streamsize i=0 ; i<MesgLen ; ++i)  GlobalLogput() << OutgoingMesg[i];
    GlobalLogput() << "|" << std::endl;
    //                   << OutgoingMesg
#endif

    const char* FirstByte = &OutgoingMesg[0];
    size_t BytesRemaining = MesgLen;
    while (BytesRemaining > 0)
    {
      const size_t PktSize = std::min(ourMaxChunkSize, BytesRemaining);
      const int SendOK = send(mySocketFD, FirstByte, PktSize, 0);
      if (SendOK < 0) CoCoA_THROW_ERROR("failed to send", "SocketStreambuf::xsputn"); //??? EOF???
      BytesRemaining -= PktSize;
      FirstByte += PktSize;
    }
    return MesgLen;
  }


  SocketStreambuf::int_type SocketStreambuf::overflow(int_type c)
  {
    if (c == EOF) return c;
    const char ch = c;
    const int SendOK = send(mySocketFD, &ch, 1, 0);
    if (SendOK < 0) return EOF;//??? CoCoA_THROW_ERROR("failed to send", "SocketStreambuf::overflow");
    return c;
  }


  SocketStreambuf::int_type SocketStreambuf::underflow()
  {
    if (gptr() < egptr()) return *gptr();
    size_t RetainChars = gptr() - eback();
    if (RetainChars > ourUngetSize) RetainChars = ourUngetSize;
    memmove(myInputBuffer+(ourUngetSize-RetainChars), gptr()-RetainChars, RetainChars);
    const ssize_t BytesRecvd = recv(mySocketFD, myInputBuffer+ourUngetSize, ourInputBufferSize-ourUngetSize, 0);
    if (BytesRecvd == -1)                                  //???
      CoCoA_THROW_ERROR("receive failed","SocketStreambuf::underflow");  //??? should give EOF???
    if (BytesRecvd == 0)
      return EOF;
    setg(myInputBuffer+(ourUngetSize-RetainChars), myInputBuffer+ourUngetSize, myInputBuffer+(ourUngetSize+BytesRecvd));
    return *gptr();
  }


/////////////////////////////////////////////////////////////////////////////


  SocketStream::SocketStream(unsigned short PortNum):
      std::iostream(&myStreambuf),
      myStreambuf(PortNum)
  {}

  SocketStream::SocketStream(const char hostname[], unsigned short PortNum):
      std::iostream(&myStreambuf),
      myStreambuf(hostname, PortNum)
  {}

} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/server/SocketStream.C,v 1.12 2024/03/25 08:04:26 abbott Exp $
// $Log: SocketStream.C,v $
// Revision 1.12  2024/03/25 08:04:26  abbott
// Summary: Modifications to work  on MinGW (& 32 bitters) (redmine 1661)
//
// Revision 1.11  2022/02/22 20:40:04  abbott
// Summary: Updated copyright message (redmine 1555)
//
// Revision 1.10  2020/06/17 15:49:31  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.9  2014/05/15 12:31:34  abbott
// Summary: Now using new files server/GlobalIO.HC (previously in CoCoA/io.H)
// Author: JAA
//
// Revision 1.8  2010/02/02 19:04:45  abbott
// Added include directive so that preprocessor macro EOF is "in scope".
//
// Revision 1.7  2007/10/30 17:14:06  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.6  2007/09/25 16:32:30  abbott
// Several minor changes to silence gcc-4.3:
//    more #includes,
//    and fixed a template problemm in RegisterServerOps.C
//
// Revision 1.5  2007/09/24 14:31:10  abbott
// Changed MaxChunkSize into a class level static constant ourMaxChunkSize.
//
// Revision 1.4  2007/06/04 12:56:50  abbott
// Added Arri's "reuse addr" suggestion as it seems to work as desired.
//
// Revision 1.3  2007/05/03 10:38:43  abbott
// Added const to various local variables.
//
// Revision 1.2  2007/04/21 19:08:44  bigatti
// -- fixed debugging printing message for strings which are not NULL
//    terminated (Ubuntu)
//
// Revision 1.1.1.1  2007/03/09 15:16:12  abbott
// Imported files
//
// Revision 1.3  2007/03/07 11:34:02  bigatti
// -- commented out logging info
//
// Revision 1.2  2007/02/12 18:39:58  bigatti
// -- added myPacketCounter and mySentBytes
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.2  2006/03/14 15:01:49  cocoa
// Improved the implementation of ring member fns for computing powers.
// Should keep Intel C++ compiler quieter too.
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
// Revision 1.3  2005/04/20 15:40:48  cocoa
// Major change: modified the standard way errors are to be signalled
// (now via a macro which records filename and line number).  Updated
// documentation in error.txt accordingly.
//
// Improved the documentation in matrix.txt (still more work to be done).
//
// Revision 1.2  2005/04/19 14:06:03  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.4  2004/11/11 14:39:58  cocoa
// -- minor changes for doxygen
//
// Revision 1.3  2004/06/29 17:10:22  cocoa
// Partially tidied use of "protected" and "private" in various
// base classes.  Checking in at the end of the day -- it works,
// and I wouldn't want it to be lost next time point's disk
// misbehaves.
//
