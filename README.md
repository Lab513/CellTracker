CellTracker
===========

Hough Transform to segment and track yeast cells





The tracker can either use the internal matlab integer linear programming (ILP) solver, or it can use the open source lp_solve ILP solver.

# Installation of the lp_solve matlab bindings

1. install the [lp_solve libary](http://sourceforge.net/projects/lpsolve/)

2. downlaod the lp_solve matlab bindings (can be found [here](http://sourceforge.net/projects/lpsolve/files/lpsolve/). Make sure to install the version corresponding to your lp_solve version.

3. compile the lp_solve matlab bindings
  * for example on ubuntu:
   1. go to the directory in which the matlab bindings have been downloaded and start matlab
   2. edit the Makefile.m or Makfile64.m (depending on your system architecture). Set the `lpsolvepath` variable to the path where the lp_solve code files are residing (if you installed the lp_solve development files via apt-get '/usr/include/lpsolve`)
   3. exectue the Makefile-script
   4. copy the mxlpsolve.mexa32 or mxlpsolve.mexa64 file to the CellTracker directory
   5. add the directory containing the lp_solve library (liblpsolve55.so) to the LD_LIBRARY_PATH (if you are using bash on ubuntu: `export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib/lp_solve`)
