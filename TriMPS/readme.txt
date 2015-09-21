----------------------------------------------------------------
Package name: TriMPS
----------------------------------------------------------------
The essential code of our paper "Efficient triangulation of Poisson-disk sampled point sets" is provided in this package.
----------------------------------------------------------------

Usage depends on the sampling mode (uniform sampling or adaptive sampling):
   uniform: trimps circle.line U 0.025
  adaptive: trimps circle.line A 0.007 12.0 1.0
please see the code for the details. Note that the boundary file should be in the ./data folder.
 
Output: the triangulation results in the format of .ps

----------------------------------------------------------------
Installation:

Prerequisites:
    CGAL library (http://www.cgal.org);
Build Tool:
    CMake (http://www.cmake.org)
Tested Compilers:
    Microsoft Visual C++ 9.0
    GNU G++

Generating Makefile:

For both Windows and Linux users:
    you can generate Visual C++ project file or Makefile by running the following
command on TriMPS directory
    cmake-gui .

For Linux user:
    A stand-alone Makefile is provided also.

For Windows user:
    If you have not any experience on CGAL installation, please
    follow the instruction in INSTALL_WIN.txt.

----------------------------------------------------------------
