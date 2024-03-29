/*
Student Name: Onur KILIÇOĞLU
Student Number: 2015400012
Compile Status: Compiling
Program Status: Working
Notes: - 
*/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h> 
int main(int argc, char** argv) {
   /* Intializes random number generator */
   srand(time(NULL));
  	// Initialize the MPI environment
	MPI_Init(NULL, NULL);
	// Find out rank, size
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	//N is the row and column size of image matrix		
	int N = 200;
	double beta,pi,gamma;
	//Read inputs from command line
	beta = atof(argv[3]);
	pi = atof(argv[4]);
	//Checking Input values
	if(beta < 0 || beta > 1){
		printf("Beta value should be between 0 and 1\n");
		exit(1);
	}
	if(pi <= 0 ||pi >= 1){
		printf("Pi value should be between 0 and 1\n");
		exit(1);
	}
	// Gamma is initialized to use in the algorithm
	gamma = 1.0/2*log((1-pi)/pi);
	// Iteration count for the algorithm 
	int iterationCount = 500000/(world_size -1);
	if (world_rank == 0){
		// Allocating memory for picture
		// Picture holds the image matrix 
		int ** picture = (int **)malloc(sizeof(int *) * N);
		for(int i=0;i<N;i++)
			picture[i] = (int *)malloc(sizeof(int) * N);
		//Getting Input
		FILE *fptr;
		fptr = fopen(argv[1],"rb");
		if (fptr == NULL){
		    printf("Error opening file!\n");
		    exit(1);
		}
		char * line = NULL;
	    size_t len = 0;
	   	ssize_t read;
		for(int i=0;i<N;i++){
			read = getline(&line, &len, fptr);
			char* token = strtok(line, " ");
			picture[i][0] =atoi(token);
			for(int j=1;j < N ;j++) {
		    	token = strtok(NULL, " ");
		    	picture[i][j] = atoi(token);
			}
		}
		free(line);
		fclose(fptr);
		//Sending blocks of size (N/p)x200 of pictures to processors 
		for(int i = 0; i < world_size-1; i++)
			for(int j = 0; j < N/(world_size - 1); j++)
				MPI_Send(picture[i*N/(world_size - 1)+j], N, MPI_INT, i+1, 0, MPI_COMM_WORLD);
		//Recieving blocks of size (N/p)x200 of pictures to processors 
		for(int i = 0; i < world_size-1; i++)
			for(int j = 0; j < N/(world_size - 1); j++)
				MPI_Recv(picture[N/(world_size - 1)*i+j], N, MPI_INT, i+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		// Writing denoised image to the output file
		FILE *f = fopen(argv[2], "w");
		if (f == NULL){
		    printf("Error opening file!\n");
		    exit(1);
		}
		// Writing output and deallocating memory
		for(int i=0;i < N;i++){
			for(int j=0; j < N; j++)
				fprintf(f, "%d ", picture[i][j]);
			fprintf(f, "\n");
			free(picture[i]);
		}
		fclose(f);
		free(picture);
	} 
	else if (world_rank > 0) {
		//Allocating memory for X of the picture used by processor
		//Allocating memory for Z of the picture used by processor
		int ** X = (int **)malloc(sizeof(int*) * N/(world_size - 1));
		int ** Z = (int **)malloc(sizeof(int*) * N/(world_size - 1));
		for(int i=0;i<N/(world_size - 1);i++){
			X[i] = (int *)malloc(sizeof(int) * N);
			Z[i] = (int *)malloc(sizeof(int) * N);
		}
		//Allocating memory for getting the bordering pixels of the block at the top
		//Allocating memory for getting the bordering pixels of the block at the bottom 
		int* top = (int *)malloc(sizeof(int) * N);
		int* bottom = (int *)malloc(sizeof(int) * N);
		//Processors receive initial block
		for(int i=0;i < N/(world_size - 1);i++)
			MPI_Recv(X[i], N, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		//Initializing Z(0)
		for(int i = 0;i < N/(world_size - 1); i++){
			for (int j = 0; j < N; j++)
				Z[i][j] = X[i][j];
		}
		//Sending top row
		if(world_rank - 1 > 0)
			MPI_Send(Z[0], N, MPI_INT, world_rank-1, 0, MPI_COMM_WORLD);
		//Sending bottom row
		if(world_rank + 1 < world_size)
			MPI_Send(Z[N/(world_size - 1) - 1], N, MPI_INT, world_rank+1, 1, MPI_COMM_WORLD);
		//Receiving top row
		if(world_rank - 1 > 0)
			MPI_Recv(top, N, MPI_INT, world_rank-1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		//Receiving bottom row
		if(world_rank + 1 < world_size)
			MPI_Recv(bottom, N, MPI_INT, world_rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		int toSend[2];
		int toReceive[2];
		//Calculating acceptance probability
		double alpha;
		double sum;
		for(int it=0;it < iterationCount; it++){
			//Selecting a random pixel
			int i = rand()%(N/(world_size - 1));
			int j = rand()%N;
			alpha = -2*gamma*Z[i][j]*X[i][j];
			// Calculating sum of products of neighboring pixels with the random pixel
			sum = 0;
			// Checking if the pixel is at the top row 
			if(i == 0){
				// The block is on the top therefore no pixel at the top border
				if(world_rank == 1){
					// The pixel is on the leftmost column 
					if(j == 0){
						sum += Z[i][j]*Z[i+1][j];
						sum += Z[i][j]*Z[i+1][j+1];
						sum += Z[i][j]*Z[i][j+1];
					}
					// The pixel is on the rightmost column 
					else if(j == N-1){
						sum += Z[i][j]*Z[i][j-1];
						sum += Z[i][j]*Z[i+1][j-1];
						sum += Z[i][j]*Z[i+1][j];
					}
					// The pixel is in the middle and has neighbors on its right and left
					else{
						sum += Z[i][j]*Z[i][j-1];
						sum += Z[i][j]*Z[i+1][j-1];
						sum += Z[i][j]*Z[i+1][j];
						sum += Z[i][j]*Z[i+1][j+1];
						sum += Z[i][j]*Z[i][j+1];
					}
				}
				// The block is not on the top therefore pixels at the top border can be used 
				else{
					// The pixel is on the leftmost column 
					if(j == 0){
						sum += Z[i][j]*Z[i+1][j];
						sum += Z[i][j]*Z[i+1][j+1];
						sum += Z[i][j]*Z[i][j+1];
						sum += Z[i][j]*top[j];
						sum += Z[i][j]*top[j+1];

					}
					// The pixel is on the rightmost column 
					else if(j == N-1){
						sum += Z[i][j]*Z[i][j-1];
						sum += Z[i][j]*Z[i+1][j-1];
						sum += Z[i][j]*Z[i+1][j];
						sum += Z[i][j]*top[j];
						sum += Z[i][j]*top[j-1];
					}
					// The pixel is in the middle and has neighbors on its right and left
					else{
						sum += Z[i][j]*Z[i][j-1];
						sum += Z[i][j]*Z[i+1][j-1];
						sum += Z[i][j]*Z[i+1][j];
						sum += Z[i][j]*Z[i+1][j+1];
						sum += Z[i][j]*Z[i][j+1];
						sum += Z[i][j]*top[j-1];
						sum += Z[i][j]*top[j];
						sum += Z[i][j]*top[j+1];
					}
				}
			
			}
			// Checking if the pixel is at the bottom row 
			else if(i == N/(world_size - 1) -1){
				// The block is at the bottom therefore no pixel at the bottom border
				if(world_rank == world_size-1){
					// The pixel is on the leftmost column 
					if(j == 0){
						sum += Z[i][j]*Z[i-1][j];
						sum += Z[i][j]*Z[i-1][j+1];
						sum += Z[i][j]*Z[i][j+1];
					}
					// The pixel is on the rightmost column
					else if(j == N-1){
						sum += Z[i][j]*Z[i][j-1];
						sum += Z[i][j]*Z[i-1][j-1];
						sum += Z[i][j]*Z[i-1][j];
					}
					// The pixel is in the middle and has neighbors on its right and left
					else{
						sum += Z[i][j]*Z[i][j-1];
						sum += Z[i][j]*Z[i-1][j-1];
						sum += Z[i][j]*Z[i-1][j];
						sum += Z[i][j]*Z[i-1][j+1];
						sum += Z[i][j]*Z[i][j+1];
					}
				}
				// The block is not at the bottom therefore pixels at the bottom border can be used 
				else{
					// The pixel is on the leftmost column 
					if(j == 0){
						sum += Z[i][j]*Z[i-1][j];
						sum += Z[i][j]*Z[i-1][j+1];
						sum += Z[i][j]*Z[i][j+1];
						sum += Z[i][j]*bottom[j];
						sum += Z[i][j]*bottom[j+1];
					}
					// The pixel is on the rightmost column 
					else if(j == N-1){
						sum += Z[i][j]*Z[i][j-1];
						sum += Z[i][j]*Z[i-1][j-1];
						sum += Z[i][j]*Z[i-1][j];
						sum += Z[i][j]*bottom[j-1];
						sum += Z[i][j]*bottom[j];
					}
					// The pixel is in the middle and has neighbors on its right and left
					else{
						sum += Z[i][j]*Z[i][j-1];
						sum += Z[i][j]*Z[i-1][j-1];
						sum += Z[i][j]*Z[i-1][j];
						sum += Z[i][j]*Z[i-1][j+1];
						sum += Z[i][j]*Z[i][j+1];
						sum += Z[i][j]*bottom[j-1];
						sum += Z[i][j]*bottom[j];
						sum += Z[i][j]*bottom[j+1];
					}
				}
			}
			// The pixel is not on the border of an adjacent block
			else{
				// The pixel is on the leftmost column
				if(j == 0){
					sum += Z[i][j]*Z[i-1][j];
					sum += Z[i][j]*Z[i-1][j+1];
					sum += Z[i][j]*Z[i][j+1];
					sum += Z[i][j]*Z[i+1][j];
					sum += Z[i][j]*Z[i+1][j+1];
				}
				// The pixel is on the rightmost column
				else if(j == N-1){
					sum += Z[i][j]*Z[i][j-1];
					sum += Z[i][j]*Z[i-1][j-1];
					sum += Z[i][j]*Z[i-1][j];
					sum += Z[i][j]*Z[i+1][j-1];
					sum += Z[i][j]*Z[i+1][j];
				}
				// The pixel is in the middle and has all neighbors
				else{
					sum += Z[i][j]*Z[i][j-1];
					sum += Z[i][j]*Z[i-1][j-1];
					sum += Z[i][j]*Z[i-1][j];
					sum += Z[i][j]*Z[i-1][j+1];
					sum += Z[i][j]*Z[i][j+1];
					sum += Z[i][j]*Z[i+1][j-1];
					sum += Z[i][j]*Z[i+1][j];
					sum += Z[i][j]*Z[i+1][j+1];
				}
			}
			alpha += -2 * beta * sum;
			// Flipping the pixel with the probability of the acceptance probability
			if(alpha > log((double)rand() / (double)RAND_MAX))
				Z[i][j] *= -1;
			//Sending top row
			if(world_rank - 1 > 0){
				toSend[0] = Z[0][j]; 
				toSend[1] = j;
				MPI_Send(&toSend, 2, MPI_INT, world_rank-1, 0, MPI_COMM_WORLD);
			}
			//Sending bottom row
			if(world_rank + 1 < world_size){
				toSend[0] = Z[N/(world_size-1)-1][j]; 
				toSend[1] = j;
				MPI_Send(&toSend, 2, MPI_INT, world_rank+1, 1, MPI_COMM_WORLD);
			}
			//Receiving top row
			if(world_rank - 1 > 0){
				MPI_Recv(&toReceive, 2, MPI_INT, world_rank-1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				top[toReceive[1]] = toReceive[0];
			}
			//Receiving bottom row
			if(world_rank + 1 < world_size){
				MPI_Recv(&toReceive, 2, MPI_INT, world_rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				bottom[toReceive[1]] = toReceive[0];
			}
		}
		// Sending the back to the master processor
		for(int i=0;i<N/(world_size - 1);i++){
			MPI_Send(Z[i], N, MPI_INT, 0, 0, MPI_COMM_WORLD);
			free(X[i]);
			free(Z[i]);
		}
		free(X);
		free(Z);
		free(top);
		free(bottom);

	}

	MPI_Finalize();
}