#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h> 

void send(int ** Z,int N,int n,int world_rank, int world_size){
	//Send bottom row
	if(world_rank + n < world_size)
		MPI_Send(Z[N/n -1], N/n, MPI_INT, world_rank + n, 0, MPI_COMM_WORLD);
	//Send top row
	if(world_rank - n > 0)
		MPI_Send(Z[0], N/n, MPI_INT, world_rank - n, 1, MPI_COMM_WORLD);
	//Send right column;
	if(world_rank % n  != 0)
		for(int i=0; i < N/n; i++)
			MPI_Send(&Z[i][0], 1, MPI_INT, world_rank + 1, 2, MPI_COMM_WORLD);
	//Send Left Column 
	if(world_rank % n  != 1)
		for(int i=0; i < N/n; i++)
			MPI_Send(&Z[i][N/n-1], 1, MPI_INT, world_rank - 1, 3, MPI_COMM_WORLD);
	//Send top left corner
	if(world_rank % n  != 1 && world_rank - n > 0)
		MPI_Send(&Z[0][0], 1, MPI_INT, world_rank -n - 1, 4, MPI_COMM_WORLD);
	//Send top right corner
	if(world_rank % n  != 0 && world_rank - n > 0)
		MPI_Send(&Z[0][N/n-1], 1, MPI_INT, world_rank -n + 1, 5, MPI_COMM_WORLD);
	//Send bottom left corner
	if(world_rank % n  != 1 && world_rank + n < world_size)
		MPI_Send(&Z[N/n-1][0], 1, MPI_INT, world_rank +n - 1, 6, MPI_COMM_WORLD);
	//Send bottom right corner
	if(world_rank % n  != 0 && world_rank + n < world_size)
		MPI_Send(&Z[N/n-1][N/n-1], 1, MPI_INT, world_rank +n +1, 7, MPI_COMM_WORLD);
}
void receive(int* top,int* bottom,int* left,int* right,int* topLeft,int* topRight,int* bottomLeft,int* bottomRight,int N,int n,int world_rank,int world_size){
	//Receive top row
	if(world_rank - n > 0)
		MPI_Recv(top, N/n, MPI_INT, world_rank - n, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	//Receive bottom row
	if(world_rank + n < world_size)
		MPI_Recv(bottom, N/n, MPI_INT, world_rank + n, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	//Receive left column;
	if(world_rank % n  != 1)
		for(int i=0; i < N/n; i++)
			MPI_Recv(&left[i], 1, MPI_INT, world_rank -1 , 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	//Receive right Column 
	if(world_rank % n  != 0)
		for(int i=0; i < N/n; i++)
			MPI_Recv(&right[i], 1, MPI_INT, world_rank + 1, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	//Receive bottom right corner
	if(world_rank % n  != 0 && world_rank + n < world_size)
		MPI_Recv(bottomRight, 1, MPI_INT, world_rank +n+ 1, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	//Receive bottom left corner
	if(world_rank % n  != 1 && world_rank + n < world_size) 
		MPI_Recv(bottomLeft, 1, MPI_INT, world_rank +n- 1, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);	
	//Receive top right corner
	if(world_rank % n  != 0 && world_rank - n > 0)
		MPI_Recv(topRight, 1, MPI_INT, world_rank-n+1, 6, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	//Receive top left corner
	if(world_rank % n  != 1 && world_rank - n > 0)
		MPI_Recv(topLeft, 1, MPI_INT, world_rank-n-1, 7, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}
void send2(int ** Z,int N,int n,int world_rank, int world_size,int i,int j){
	int toSend[2];
	//Send bottom row
	if(world_rank + n < world_size){
		toSend[0] = Z[N/n-1][j];
		toSend[1] = j;
		MPI_Send(&toSend, 2, MPI_INT, world_rank + n, 0, MPI_COMM_WORLD);
	}
	//Send top row
	if(world_rank - n > 0){
		toSend[0] = Z[0][j];
		toSend[1] = j;
		MPI_Send(&toSend, 2, MPI_INT, world_rank - n, 1, MPI_COMM_WORLD);
	}
	//Send right column;
	if(world_rank % n  != 0){
		toSend[0] = Z[i][N/n-1];
		toSend[1] = i;
		MPI_Send(&toSend, 1, MPI_INT, world_rank + 1, 2, MPI_COMM_WORLD);
	}
	//Send Left Column 
	if(world_rank % n  != 1){
		toSend[0] = Z[i][0];
		toSend[1] = i;
		MPI_Send(&toSend, 1, MPI_INT, world_rank - 1, 3, MPI_COMM_WORLD);
	}
	//Send top left corner
	if(world_rank % n  != 1 && world_rank - n > 0)
		MPI_Send(&Z[0][0], 1, MPI_INT, world_rank -n - 1, 4, MPI_COMM_WORLD);
	//Send top right corner
	if(world_rank % n  != 0 && world_rank - n > 0)
		MPI_Send(&Z[0][N/n-1], 1, MPI_INT, world_rank -n + 1, 5, MPI_COMM_WORLD);
	//Send bottom left corner
	if(world_rank % n  != 1 && world_rank + n < world_size)
		MPI_Send(&Z[N/n-1][0], 1, MPI_INT, world_rank +n - 1, 6, MPI_COMM_WORLD);
	//Send bottom right corner
	if(world_rank % n  != 0 && world_rank + n < world_size)
		MPI_Send(&Z[N/n-1][N/n-1], 1, MPI_INT, world_rank +n +1, 7, MPI_COMM_WORLD);
}
void receive2(int* top,int* bottom,int* left,int* right,int* topLeft,int* topRight,int* bottomLeft,int* bottomRight,int N,int n,int world_rank,int world_size,int i,int j){
	int toReceive[2];
	//Receive top row
	if(world_rank - n > 0){
		MPI_Recv(&toReceive, 2, MPI_INT, world_rank - n, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		top[toReceive[1]] =toReceive[0];
	}
	//Receive bottom row
	if(world_rank + n < world_size){
		MPI_Recv(&toReceive, 2,MPI_INT, world_rank + n, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		bottom[toReceive[1]] =toReceive[0];
	}
	//Receive left column;
	if(world_rank % n  != 1){
		MPI_Recv(&toReceive, 2, MPI_INT, world_rank -1 , 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		left[toReceive[1]] =toReceive[0];
	}
	//Receive right Column 
	if(world_rank % n  != 0){
		MPI_Recv(&toReceive, 2, MPI_INT, world_rank + 1, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		right[toReceive[1]] =toReceive[0];
	}
	//Receive bottom right corner
	if(world_rank % n  != 0 && world_rank + n < world_size)
		MPI_Recv(bottomRight, 1, MPI_INT, world_rank +n+ 1, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	//Receive bottom left corner
	if(world_rank % n  != 1 && world_rank + n < world_size) 
		MPI_Recv(bottomLeft, 1, MPI_INT, world_rank +n- 1, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);	
	//Receive top right corner
	if(world_rank % n  != 0 && world_rank - n > 0)
		MPI_Recv(topRight, 1, MPI_INT, world_rank-n+1, 6, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	//Receive top left corner
	if(world_rank % n  != 1 && world_rank - n > 0)
		MPI_Recv(topLeft, 1, MPI_INT, world_rank-n-1, 7, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}
int main(int argc, char** argv) {
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
	int n = (int)sqrt(world_size-1); 
	gamma = 1.0/2*log((1-pi)/pi);
	int iterationCount = 500000/(world_size -1);
	if (world_rank == 0){
		// Allocating memory for picture
		int *** picture = (int ***)malloc(sizeof(int **) * (world_size-1));
		for(int i=0;i < (world_size-1);i++){
			picture[i] = (int **)malloc(sizeof(int*) * N/n);
			for(int j =0;j < N/n; j++)
				picture[i][j] = (int *)malloc(sizeof(int) * N/n);
		}
		//TODO get Input
		FILE *fptr;
		fptr = fopen(argv[1],"rb");
		char * line = NULL;
	    size_t len = 0;
	   	ssize_t read;
		for(int i=0;i<N;i++){
			//TODO control
			if((read = getline(&line, &len, fptr)) != -1)
				;
			
			char* token = strtok(line, " ");
			//TODO error check
			picture[i/(N/n)][i%(N/n)][0] =atoi(token);
			for(int j=1;j < N ;j++) {
		    	token = strtok(NULL, " ");
		    	picture[i/(N/n)*n+ j/(N/n)][i%(N/n)][j%(N/n)] = atoi(token);
			}
		}
		free(line);
		fclose(fptr);
		

		//Sending frames of size (N/p)x200 of pictures to processors 
		for(int i = 0; i < world_size-1; i++)
			for(int j = 0; j < N/n; j++)
				MPI_Send(picture[i][j], N/n, MPI_INT, i+1, 0, MPI_COMM_WORLD);
		//Recieving frames of size (N/p)x200 of pictures to processors 
		for(int i = 0; i < world_size-1; i++)
			for(int j = 0; j < N/n; j++)
				MPI_Recv(picture[i][j], N/n, MPI_INT, i+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		FILE *f = fopen(argv[2], "w");
		if (f == NULL){
		    printf("Error opening file!\n");
		    exit(1);
		}
		for(int i=0;i < N;i++){
			for(int j=0; j < N; j++){
				fprintf(f, "%d ", picture[i/(N/n)*n+ j/(N/n)][i%(N/n)][j%(N/n)]);
			}
			fprintf(f, "\n");
		}
		fclose(f);
		for(int i=0;i < n;i++){
			for(int j =0;j < N/n; j++)
				free(picture[i][j]);
			free(picture[i]);
		}
		free(picture);
	}
	else if (world_rank > 0) {
		//Initializing random number generaror.
		time_t t;
		srand((int)time(&t) % world_rank);
		//Allocating memory for X of the picture used by processor
		int ** X = (int **)malloc(sizeof(int*) * N/n);
		int ** Z = (int **)malloc(sizeof(int*) * N/n);
		for(int i=0;i<N/n;i++){
			X[i] = (int *)malloc(sizeof(int) * N/n);
			Z[i] = (int *)malloc(sizeof(int) * N/n);
		}
		int* top = (int *)malloc(sizeof(int) * N/n);
		int* bottom = (int *)malloc(sizeof(int) * N/n);
		int* right = (int *)malloc(sizeof(int) * N/n);
		int* left = (int *)malloc(sizeof(int) * N/n);
		int topLeft, topRight, bottomLeft, bottomRight;
		//Processors receive initial frame
		for(int i=0;i < N/n;i++)
			MPI_Recv(X[i], N/n, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		//Initializing Z(0)
		for(int i = 0;i < N/n; i++){
			for (int j = 0; j < N/n; j++){
				Z[i][j] = X[i][j];
			}
		}
		if(world_size >2){
			send(Z,N,n,world_rank,world_size);
			receive(top,bottom,left,right,&topLeft,&topRight,&bottomLeft,&bottomRight,N,n,world_rank,world_size);
		}
		//Calculating acceptance probability
		double alpha;
		double sum;
		for(int it=0;it < iterationCount; it++){
			//Selecting a random pixel
			int i = rand()%(N/n);
			int j = rand()%(N/n);
			alpha = -2*gamma*Z[i][j]*X[i][j];
			sum = 0;
			if(i == 0){
				if(world_rank - n <= 0){
					if(j == 0){
						sum += Z[i][j]*Z[i+1][j];
						sum += Z[i][j]*Z[i+1][j+1];
						sum += Z[i][j]*Z[i][j+1];
						if(world_rank %n != 1){
							sum += Z[i][j]*left[i];
							sum += Z[i][j]*left[i+1];
						}
					}
					else if(j == N-1){
						sum += Z[i][j]*Z[i][j-1];
						sum += Z[i][j]*Z[i+1][j-1];
						sum += Z[i][j]*Z[i+1][j];
						if(world_rank %n != 0){
							sum += Z[i][j]*right[i];
							sum += Z[i][j]*right[i+1];
						}
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
						if(world_rank %n != 1){
							sum += Z[i][j]*left[i];
							sum += Z[i][j]*left[i+1];
							sum += topLeft;
						}

					}
					else if(j == N-1){
						sum += Z[i][j]*Z[i][j-1];
						sum += Z[i][j]*Z[i+1][j-1];
						sum += Z[i][j]*Z[i+1][j];
						sum += Z[i][j]*top[j];
						sum += Z[i][j]*top[j-1];
						if(world_rank %n != 0){
							sum += Z[i][j]*right[i];
							sum += Z[i][j]*right[i+1];
							sum += topRight;
						}
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
			else if(i == N/n -1){
				if(world_rank + n > world_size){
					if(j == 0){
						sum += Z[i][j]*Z[i-1][j];
						sum += Z[i][j]*Z[i-1][j+1];
						sum += Z[i][j]*Z[i][j+1];
						if(world_rank %n != 1){
							sum += Z[i][j]*left[i-1];
							sum += Z[i][j]*left[i];
						}
					}
					else if(j == N-1){
						sum += Z[i][j]*Z[i][j-1];
						sum += Z[i][j]*Z[i-1][j-1];
						sum += Z[i][j]*Z[i-1][j];
						if(world_rank %n != 0){
							sum += Z[i][j]*right[i-1];
							sum += Z[i][j]*right[i];
						}
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
						if(world_rank %n != 1){
							sum += Z[i][j]*left[i-1];
							sum += Z[i][j]*left[i];
							sum += bottomLeft;
						}
					}
					else if(j == N-1){
						sum += Z[i][j]*Z[i][j-1];
						sum += Z[i][j]*Z[i-1][j-1];
						sum += Z[i][j]*Z[i-1][j];
						sum += Z[i][j]*bottom[j-1];
						sum += Z[i][j]*bottom[j];
						if(world_rank %n != 0){
							sum += Z[i][j]*right[i-1];
							sum += Z[i][j]*right[i];
							sum += bottomRight;
						}
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
					if(world_rank %n != 1){
						sum += Z[i][j]*left[i-1];
						sum += Z[i][j]*left[i];
						sum += Z[i][j]*left[i+1];
					}
				}
				else if(j == N-1){
					sum += Z[i][j]*Z[i][j-1];
					sum += Z[i][j]*Z[i-1][j-1];
					sum += Z[i][j]*Z[i-1][j];
					sum += Z[i][j]*Z[i+1][j-1];
					sum += Z[i][j]*Z[i+1][j];
					if(world_rank %n != 0){
						sum += Z[i][j]*right[i-1];
						sum += Z[i][j]*right[i];
						sum += Z[i][j]*right[i+1];
					}
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
			if(alpha > log((double)rand() / (double)RAND_MAX))
				Z[i][j] *= -1;
			if(world_size >2){
				send2(Z,N,n,world_rank,world_size,i,j);
				receive2(top,bottom,left,right,&topLeft,&topRight,&bottomLeft,&bottomRight,N,n,world_rank,world_size,i,j);
			}
		}
		for(int i=0;i<N/n;i++){
			MPI_Send(Z[i], N/n, MPI_INT, 0, 0, MPI_COMM_WORLD);
			free(X[i]);
			free(Z[i]);
		}
		free(X);
		free(Z);
		free(top);
		free(bottom);
		free(right);
		free(left);
	}

	MPI_Finalize();
}