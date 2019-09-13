#ifndef VTKH_EXPORTS_H
#define VTKH_EXPORTS_H

#if __GNUC__ >= 4 && (defined(VTKH_COMPILING_FLAG))
  #define VTKH_API __attribute__ ((visibility("default")))
#else
  #define VTKH_API /* hidden by default */
#endif

#endif
