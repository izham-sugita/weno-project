# weno-project
Research on WENO implementation

This is a simple WENO implementation project. The purpose is to demonstrate the WENO implementation as barebone as
possible.

The file advec-weno.cpp is the barebone example of linear advection simulation
using 4th-order Runge-Kutta and 5th-order WENO. The interpolation was used;
instead of the reconstruction method (yes, the difference is subtle; but
significant for conservation law). 

#update 2019/12/27
File : advec-weno-2D.cpp is an implementation of two-dimensional advection equation using
5th-order WENO reconstruction. WENO reconstruction produces higher resolution results compared
to interpolation-based. Currently working on the central WENO to further simplify WENO
implementation. The final target is to build a GPU kernel for central WENO.
