p_lasmerge 

Background:

The LAStools lasmerge application and supporting LASlib library were extended 
with MPI to allow the application to be run in parallel on a compute 
cluster. The goal was an application that would scale to arbitrarily large 
input, limited only by the amount of disk space needed to store the input and 
output file. No intermediate files are generated and individual process
memory needs are not determined by the size of the input or output.

The algorithm_and_results.pdf slide presentation contains a description 
of the implementation and test results of merging 785 LAS files on into one
111GB las file on a compute cluster.

Dependencies:
An MPI implementation must be installed. OpenMPI 1.6 and 1.8 are known to work
and were used in development and testing. mpic++ must be found in your PATH. 

Install:

git clone https://github.com/jwend/p_lasmerge
cd p_lasmerge
make

The p_lasmerge executable is in the bin directory.

Test:

mpirun -n 4 bin/p_lasmerge -i data/test*.las -o merged.las

Limitations and Supported Features:

p_lasmerge works only with LAS version 1.0, 1.1, and 1.2 and produces only 
corresponding LAS file output. Version 1.2 was tested most extensively with 
up to the 111 GB file size input and output. This implementation does not 
support filtering or transforming the points. See:
https://github.com/jwend/p_las2las if those functions are needed after a run
of p_lasmerge.  








