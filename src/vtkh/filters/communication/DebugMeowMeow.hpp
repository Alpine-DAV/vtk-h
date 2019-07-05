#ifndef VTK_H_DEBUG_MEOW_MEOW_HPP
#define VTK_H_DEBUG_MEOW_MEOW_HPP

#include <vtkh/filters/util.hpp>

extern std::ofstream dbg;
extern std::ofstream wdbg;

#ifdef TRACE_DEBUG
#define DBG(x) dbg<<x
#define WDBG(x) wdbg<<x
#else
#define DBG(x)
#define WDBG(x)
#endif

#endif
