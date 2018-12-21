#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h> 

/* send function takes the matrix of the block of image that is processed by a block, matrix row size N,
* number of blocks in a row n, the rank of the slave processor and the world size and sends initial  
* version of all the bordering pixels to the necessary processor 
*/
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
			MPI_Send(&Z[i][N/n-1], 1, MPI_INT, world_rank + 1, 2, MPI_COMM_WORLD);
	//Send Left Column 
	if(world_rank % n  != 1)
		for(int i=0; i < N/n; i++)
			MPI_Send(&Z[i][0], 1, MPI_INT, world_rank - 1, 3, MPI_COMM_WORLD);
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
/* receive function takes pointers to the arrays that hold the information about bordering pixels top,bottom,left,right,
* the pointers to the variables that hold the information about bordering pixels topLeft,topRight,bottomLeft,bottomRight,
* matrix row size N, number of blocks in a row n, the rank of the slave processor and the world size
* and receives the initial version of all bordering pixels from the necessary processor
*/
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

/* send2 function takes the matrix of the block of image that is processed by a block, matrix row size N,
* number of blocks in a row n, the rank of the slave processor, the world size and indices of the most recently changed pixel 
* and sends the version of the possibly changed bordering pixel to the necessary processor with the index of the possible change
*/
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

/* receive2 function takes pointers to the arrays that hold the information about bordering pixels top,bottom,left,right,
* the pointers to the variables that hold the information about bordering pixels topLeft,topRight,bottomLeft,bottomRight,
* matrix row size N, number of blocks in a row n, the rank of the slave processor and the world size
* and receives the version of the possibly changed bordering pixel from the necessary processor with the index of the possible change
* then changes the pixel at the index to the most recent one
*/
void receive2(int* top,int* bottom,int* left,int* right,int* topLeft,int* topRight,int* bottomLeft,int* bottomRight,int N,int n,int world_rank,int world_size){
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
	//N is the row and column size of image matrix
	int N = 200;
	double beta,pi,gamma;
	//Read inputs from command line
	beta = atof(argv[3]);
	pi = atof(argv[4]);
	int n = (int)sqrt(world_size-1);
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
		// Picture contains blocks to be sent to slave processors
		int *** picture = (int ***)malloc(sizeof(int **) * (world_size-1));
		for(int i=0;i < (world_size-1);i++){
			picture[i] = (int **)malloc(sizeof(int*) * N/n);
			for(int j =0;j < N/n; j++)
				picture[i][j] = (int *)malloc(sizeof(int) * N/n);
		}
		// Getting Input
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
			// Putting each pixel of the image in the right block 
			picture[i/(N/n)][i%(N/n)][0] =atoi(token);
			for(int j=1;j < N ;j++) {
		    	token = strtok(NULL, " ");
		    	picture[i/(N/n)*n+ j/(N/n)][i%(N/n)][j%(N/n)] = atoi(token);
			}
		}
		free(line);
		fclose(fptr);
		

		//Sending blocks of size (N/n)x(N/n) of pictures to processors 
		for(int i = 0; i < world_size-1; i++)
			for(int j = 0; j < N/n; j++)
				MPI_Send(picture[i][j], N/n, MPI_INT, i+1, 0, MPI_COMM_WORLD);
		//Recieving blocks of size (N/n)x(N/n) of pictures to processors 
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
			// Calculating sum of products of neighboring pixels with the random pixel
			sum = 0;
			// Checking if the pixel is at the top row 
			if(i == 0){
				// The block is on the top therefore no pixel at the top border
				if(world_rank - n <= 0){
					// The pixel is on the leftmost column 
					if(j == 0){
						sum += Z[i][j]*Z[i+1][j];
						sum += Z[i][j]*Z[i+1][j+1];
						sum += Z[i][j]*Z[i][j+1];
						// If the block is not on the leftmost part of picture pixels from left border can be used 
						if(world_rank %n != 1){
							sum += Z[i][j]*left[i];
							sum += Z[i][j]*left[i+1];
						}
					}
					// The pixel is on the rightmost column 
					else if(j == N-1){
						sum += Z[i][j]*Z[i][j-1];
						sum += Z[i][j]*Z[i+1][j-1];
						sum += Z[i][j]*Z[i+1][j];
						// If the block is not on the rightmost part of picture pixels from right border can be used
						if(world_rank %n != 0){
							sum += Z[i][j]*right[i];
							sum += Z[i][j]*right[i+1];
						}
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
				// The block is not on the top therefore pixel at the top border can be used 
				else{
					// The pixel is on the leftmost column 
					if(j == 0){
						sum += Z[i][j]*Z[i+1][j];
						sum += Z[i][j]*Z[i+1][j+1];
						sum += Z[i][j]*Z[i][j+1];
						sum += Z[i][j]*top[j];
						sum += Z[i][j]*top[j+1];
						// If the block is not on the leftmost part of picture pixels from left border can be used
						if(world_rank %n != 1){
							sum += Z[i][j]*left[i];
							sum += Z[i][j]*left[i+1];
							sum += Z[i][j]*topLeft;
						}

					}
					// The pixel is on the rightmost column 
					else if(j == N-1){
						sum += Z[i][j]*Z[i][j-1];
						sum += Z[i][j]*Z[i+1][j-1];
						sum += Z[i][j]*Z[i+1][j];
						sum += Z[i][j]*top[j];
						sum += Z[i][j]*top[j-1];
						// If the block is not on the rightmost part of picture pixels from right border can be used
						if(world_rank %n != 0){
							sum += Z[i][j]*right[i];
							sum += Z[i][j]*right[i+1];
							sum += Z[i][j]*topRight;
						}
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
			else if(i == N/n -1){
				// The block is at the bottom therefore no pixel at the bottom border
				if(world_rank + n > world_size){
					// The pixel is on the leftmost column 
					if(j == 0){
						sum += Z[i][j]*Z[i-1][j];
						sum += Z[i][j]*Z[i-1][j+1];
						sum += Z[i][j]*Z[i][j+1];
						// If the block is not on the leftmost part of picture pixels from left border can be used
						if(world_rank %n != 1){
							sum += Z[i][j]*left[i-1];
							sum += Z[i][j]*left[i];
						}
					}
					// The pixel is on the rightmost column 
					else if(j == N-1){
						sum += Z[i][j]*Z[i][j-1];
						sum += Z[i][j]*Z[i-1][j-1];
						sum += Z[i][j]*Z[i-1][j];
						// If the block is not on the rightmost part of picture pixels from right border can be used
						if(world_rank %n != 0){
							sum += Z[i][j]*right[i-1];
							sum += Z[i][j]*right[i];
						}
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
						// If the block is not on the leftmost part of picture pixels from left border can be used
						if(world_rank %n != 1){
							sum += Z[i][j]*left[i-1];
							sum += Z[i][j]*left[i];
							sum += Z[i][j]*bottomLeft;
						}
					}
					// The pixel is on the rightmost column 
					else if(j == N-1){
						sum += Z[i][j]*Z[i][j-1];
						sum += Z[i][j]*Z[i-1][j-1];
						sum += Z[i][j]*Z[i-1][j];
						sum += Z[i][j]*bottom[j-1];
						sum += Z[i][j]*bottom[j];
						// If the block is not on the rightmost part of picture pixels from right border can be used
						if(world_rank %n != 0){
							sum += Z[i][j]*right[i-1];
							sum += Z[i][j]*right[i];
							sum += Z[i][j]*bottomRight;
						}
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
			// The pixel is not on the border of an adjacent block from top or bottom
			else{
				// The pixel is on the leftmost column 
				if(j == 0){
					sum += Z[i][j]*Z[i-1][j];
					sum += Z[i][j]*Z[i-1][j+1];
					sum += Z[i][j]*Z[i][j+1];
					sum += Z[i][j]*Z[i+1][j];
					sum += Z[i][j]*Z[i+1][j+1];
					// If the block is not on the leftmost part of picture pixels from left border can be used
					if(world_rank %n != 1){
						sum += Z[i][j]*left[i-1];
						sum += Z[i][j]*left[i];
						sum += Z[i][j]*left[i+1];
					}
				}
				// The pixel is on the rightmost column 
				else if(j == N-1){
					sum += Z[i][j]*Z[i][j-1];
					sum += Z[i][j]*Z[i-1][j-1];
					sum += Z[i][j]*Z[i-1][j];
					sum += Z[i][j]*Z[i+1][j-1];
					sum += Z[i][j]*Z[i+1][j];
					// If the block is not on the rightmost part of picture pixels from right border can be used
					if(world_rank %n != 0){
						sum += Z[i][j]*right[i-1];
						sum += Z[i][j]*right[i];
						sum += Z[i][j]*right[i+1];
					}
				}
				// The pixel is in the middle and has neighbors on its right and left
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
			// Sending possible changes to the adjacent processors if exists
			if(world_size >2){
				send2(Z,N,n,world_rank,world_size,i,j);
				receive2(top,bottom,left,right,&topLeft,&topRight,&bottomLeft,&bottomRight,N,n,world_rank,world_size);
			}

		}
		// Sending the back to the master processor
		// Deallocating memory
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