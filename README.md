# mpi-project
Parallel Programming of Domain Decomposition Algorithm

This project applies domain decomposition to solve an ODE 
in parallel using MPI. To compile the source:

1) Install tcl and gcc if not already.
2) Refer to the readme file for m6378-mpi.tar.gz, 
install, and test using the example.
3) Edit p2.sh and specify which c source file to compile 
(e.g. SRCFILE=$HOMEDIR/p2_v2-time.c).
4) Run p2run.sh and gcc should compile the source to p2.

To run:
1) In a terminal, m6378start
2) m6378run  number-of-processes  p2  local-grid-size  number-of-iterations  *number-of-outer-iterations*
   (*v2 source code only)

To end:
1) m6378stop

Computed solutions for a single process should be consistent with those of parallel processes.
