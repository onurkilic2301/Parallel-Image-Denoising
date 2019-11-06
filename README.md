# Parallel Image Denoising 
This program is implemented to run on UNIX systems. 
To be able to run the program user should install OpenMPI version 1.4.4. 
The program is run through a command line interface such as bash. 
To start the program first the user should compile the program. 
To compile the program, open the terminal and change directory to the source code’s directory. 

Then to compile the program run the command:<br />
$ mpicc -g main.c -o pr

To compile the second approach program run the command:<br />
$ mpicc -g task2.c -o pr

To run the program, run the command:<br />
$ mpiexec -n \<num-of-process> ./pr \<input-file> \<output-file> \<beta-value> \<pi-value>
  
Parameters:<br />
  \<num-of-processor>: Number of processes to run the program in parallel.<br />
  \<input-file>: Path to the text file including the matrix of the image of size 200x200 to be denoised by the program.<br />
  \<output-file>: Path to the text file to write the output of the program which is the denoised version of image.<br />
  \<beta-value>: β parameter from the Ising Model which is a number between 0 and 1. <pi-value>: π parameter from the Ising Model which is a number between 0 and 1.<br />
  The program terminates itself after execution.<br />
  
  Check Report for more explanation.
