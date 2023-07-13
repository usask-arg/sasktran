
#if !defined(NXBASE_THREADHINCLUDE)
#define NXBASE_THREADHINCLUDE

#include <mutex>
#include <condition_variable>
#include <thread>
#ifdef NX_WINDOWS
  #include "module/system/multithread/nxsyncobjects.h"
  #include "module/system/multithread/nxmemoryfifo.h"
#endif

#include "module/system/multithread/nxworkerthread.h"			// nxworkerthread.h, added May 28, 2009, support for boost::threads, (developed with boost 1.39.0 )
#endif

