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
	int N = 200;
	double beta,pi,gamma;
	beta = atof(argv[3]);
	pi = atof(argv[4]);
 	//TODO get input
	gamma = 1.0/2*log((1-pi)/pi);
	int iterationCount = 500000/(world_size -1);
	if (world_rank == 0){
		// Allocating memory for picture
		int ** picture = (int **)malloc(sizeof(int *) * N);
		for(int i=0;i<N;i++)
			picture[i] = (int *)malloc(sizeof(int) * N);
		//getting Input
		FILE *fptr;
		fptr = fopen(argv[1],"rb");
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
		//Sending frames of size (N/p)x200 of pictures to processors 
		for(int i = 0; i < world_size-1; i++)
			for(int j = 0; j < N/(world_size - 1); j++)
				MPI_Send(picture[i*N/(world_size - 1)+j], N, MPI_INT, i+1, 0, MPI_COMM_WORLD);
		//Recieving frames of size (N/p)x200 of pictures to processors 
		for(int i = 0; i < world_size-1; i++)
			for(int j = 0; j < N/(world_size - 1); j++)
				MPI_Recv(picture[N/(world_size - 1)*i+j], N, MPI_INT, i+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		FILE *f = fopen(argv[2], "w");
		if (f == NULL){
		    printf("Error opening file!\n");
		    exit(1);
		}
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
		int ** X = (int **)malloc(sizeof(int*) * N/(world_size - 1));
		int ** Z = (int **)malloc(sizeof(int*) * N/(world_size - 1));
		for(int i=0;i<N/(world_size - 1);i++){
			X[i] = (int *)malloc(sizeof(int) * N);
			Z[i] = (int *)malloc(sizeof(int) * N);
		}
		int* top = (int *)malloc(sizeof(int) * N);
		int* bottom = (int *)malloc(sizeof(int) * N);
		//Processors receive initial frame
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
			sum = 0;
			if(i == 0){
				if(world_rank == 1){
					if(j == 0){
						sum += Z[i][j]*Z[i+1][j];
						sum += Z[i][j]*Z[i+1][j+1];
						sum += Z[i][j]*Z[i][j+1];
					}
					else if(j == N-1){
						sum += Z[i][j]*Z[i][j-1];
						sum += Z[i][j]*Z[i+1][j-1];
						sum += Z[i][j]*Z[i+1][j];
					}
					else{
						sum += Z[i][j]*Z[i][j-1];
						sum += Z[i][j]*Z[i+1][j-1];
						sum += Z[i][j]*Z[i+1][j];
						sum += Z[i][j]*Z[i+1][j+1];
						sum += Z[i][j]*Z[i][j+1];
					}
				}
				else{
					if(j == 0){
						sum += Z[i][j]*Z[i+1][j];
						sum += Z[i][j]*Z[i+1][j+1];
						sum += Z[i][j]*Z[i][j+1];
						sum += Z[i][j]*top[j];
						sum += Z[i][j]*top[j+1];

					}
					else if(j == N-1){
						sum += Z[i][j]*Z[i][j-1];
						sum += Z[i][j]*Z[i+1][j-1];
						sum += Z[i][j]*Z[i+1][j];
						sum += Z[i][j]*top[j];
						sum += Z[i][j]*top[j-1];
					}
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
			else if(i == N/(world_size - 1) -1){
				if(world_rank == world_size-1){
					if(j == 0){
						sum += Z[i][j]*Z[i-1][j];
						sum += Z[i][j]*Z[i-1][j+1];
						sum += Z[i][j]*Z[i][j+1];
					}
					else if(j == N-1){
						sum += Z[i][j]*Z[i][j-1];
						sum += Z[i][j]*Z[i-1][j-1];
						sum += Z[i][j]*Z[i-1][j];
					}
					else{
						sum += Z[i][j]*Z[i][j-1];
						sum += Z[i][j]*Z[i-1][j-1];
						sum += Z[i][j]*Z[i-1][j];
						sum += Z[i][j]*Z[i-1][j+1];
						sum += Z[i][j]*Z[i][j+1];
					}
				}
				else{
					if(j == 0){
						sum += Z[i][j]*Z[i-1][j];
						sum += Z[i][j]*Z[i-1][j+1];
						sum += Z[i][j]*Z[i][j+1];
						sum += Z[i][j]*bottom[j];
						sum += Z[i][j]*bottom[j+1];
					}
					else if(j == N-1){
						sum += Z[i][j]*Z[i][j-1];
						sum += Z[i][j]*Z[i-1][j-1];
						sum += Z[i][j]*Z[i-1][j];
						sum += Z[i][j]*bottom[j-1];
						sum += Z[i][j]*bottom[j];
					}
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
			else{
				if(j == 0){
					sum += Z[i][j]*Z[i-1][j];
					sum += Z[i][j]*Z[i-1][j+1];
					sum += Z[i][j]*Z[i][j+1];
					sum += Z[i][j]*Z[i+1][j];
					sum += Z[i][j]*Z[i+1][j+1];
				}
				else if(j == N-1){
					sum += Z[i][j]*Z[i][j-1];
					sum += Z[i][j]*Z[i-1][j-1];
					sum += Z[i][j]*Z[i-1][j];
					sum += Z[i][j]*Z[i+1][j-1];
					sum += Z[i][j]*Z[i+1][j];
				}
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
			//printf("%f\n", alpha);
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