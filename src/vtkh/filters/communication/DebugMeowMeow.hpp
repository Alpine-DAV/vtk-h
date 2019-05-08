#ifndef VTK_H_DEBUG_MEOW_MEOW_HPP
#define VTK_H_DEBUG_MEOW_MEOW_HPP

#define DEBUG_LOG 1
#define DEBUG_PRINT 1

#include "util.hpp"
extern ofstream dbg;

#ifdef TRACE_DEBUG
#define DBG(x) dbg<<x
#else
#define DBG(x)
#endif

#ifdef  DEBUG_PRINT
#define DPRINT(x) std::cout<<x
#else
#define DPRINT(x)
#endif

#endif
