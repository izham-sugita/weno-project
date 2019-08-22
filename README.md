# weno-project
Research on WENO implementation

This is a simple WENO implementation project. The purpose is to demonstrate the WENO implementation as barebone as
possible.

The file advec-weno.cpp is the barebone example of linear advection simulation
using 4th-order Runge-Kutta and 5th-order WENO. The interpolation was used;
instead of the reconstruction method (yes, the difference is subtle; but
significant for conservation law). 
