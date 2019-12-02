#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <time.h>
#include <string.h>
#include <math.h>

#include <definitions.h>
#include <utilities.h>

void parseArgv(int argc, char *argv[], char *DataFilePath, int *HypercubeDimension);
void checkArgv(int Size, int HypercubeDimension, char *DataFilePath);
void notifyErrorAndClose(MPI_Comm comm, short ErrorCode);
void printUsage();

void neightborsHypercube(int HypercubeDimension, int Rank, int *Neightbors);

void dispatchDataToNodes(char *DataFilePath, int Size, MPI_Comm comm, int *RootData);

void calcMaxHypercubeNetwork(int Rank, int HypercubeDimension,int *Neightbors, int LocalData, int *Max); 

int main(int argc, char *argv[])
{
	int Rank, Size, Data, Max, HypercubeDimension;
	clock_t ExecutionTime;
	int *Neightbors;
	MPI_Status status;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
	MPI_Comm_size(MPI_COMM_WORLD, &Size);

	if(Rank==ROOT)
	{
		char DataFilePath[MAX_DATA_BLOCK];
		int HypercubeDimension;
		
		parseArgv(argc, argv, DataFilePath, &HypercubeDimension);
		
		checkArgv(Size, HypercubeDimension, DataFilePath);

		short StatusCode = STATUS_OK;
		
		if(MPI_Bcast(&StatusCode, 1, MPI_SHORT, ROOT, MPI_COMM_WORLD)!=MPI_SUCCESS)
		{
			fprintf(stderr,"Fail in processes notification\n");
			exit(EXIT_FAILURE);
		}
		
		dispatchDataToNodes(DataFilePath, Size, MPI_COMM_WORLD, &Data);
		
	}
	else
	{
		short StatusCode;
		MPI_Bcast(&StatusCode, 1, MPI_SHORT, ROOT, MPI_COMM_WORLD);
		
		if(StatusCode<0)
		{
			fprintf(stderr, "Stop process %d\n",Rank);
			MPI_Finalize();
			exit(EXIT_FAILURE);
		}
		
		MPI_Recv(&Data, 1, MPI_INT, ROOT, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	}
	
	HypercubeDimension=atoi(argv[2]);
	
	Neightbors=malloc(HypercubeDimension*sizeof(int));
	
	neightborsHypercube(HypercubeDimension, Rank, Neightbors);
	
	ExecutionTime=clock();	
	calcMaxHypercubeNetwork(Rank,HypercubeDimension, Neightbors, Data, &Max);
	ExecutionTime=clock()-ExecutionTime;
	
	if(Rank==ROOT)
	{
		printf("Execution Time of The Algorithm: %f\n", (ExecutionTime/(double)CLOCKS_PER_SEC));
		printf("Max of the network: %d\n", Max);
	}
	
	free(Neightbors);
	MPI_Finalize();
	
	return EXIT_SUCCESS;
}
	
void parseArgv(int argc, char *argv[], char *DataFilePath, int *HypercubeDimension)
{
	
	if(argc!=3)
	{
		printUsage();
		notifyErrorAndClose(MPI_COMM_WORLD, WRONG_ARGUMENT_ERROR_CODE);
	}
	
	strcpy(DataFilePath,argv[1]);
	*HypercubeDimension=atoi(argv[2]);
}

void checkArgv(int Size, int HypercubeDimension, char *DataFilePath)
{
	
	if(Size!=(pow(2,HypercubeDimension)))
	{
		fprintf(stderr,"Error, wrong Hypercube Dimension Size\n");
		notifyErrorAndClose(MPI_COMM_WORLD, WRONG_NUMBER_OF_LAUNCHED_PROCESSES);
	}
	
	if(isFile(DataFilePath)==-1)
	{
		fprintf(stderr, "Error, in the specified file path\n"); 
		notifyErrorAndClose(MPI_COMM_WORLD, WRONG_FILE_PATH);
	}
	
	if(countDataNumber(DataFilePath)<Size)
	{
		fprintf(stderr, "Error, not enought data for Hypercube\n");
		notifyErrorAndClose(MPI_COMM_WORLD, NOT_ENOUGHT_DATA_FOR_PROCESSES);
	}
	
}

void printUsage()
{
	printf("Correct use of the program: mpirun -n <Number of processes to launch> ./Hypercube <Data file path> <Hypercube Dimension>\n");
}

void notifyErrorAndClose(MPI_Comm comm, short ErrorCode)
{
	if(MPI_Bcast(&ErrorCode, 1, MPI_SHORT, ROOT, comm)!=MPI_SUCCESS)
	{
		fprintf(stderr,"Fail in processes notification\n");
	}
	MPI_Finalize();
	exit(EXIT_FAILURE);
}

void dispatchDataToNodes(char *DataFilePath, int Size, MPI_Comm comm, int *RootData)
{
	FILE *DataFile;

	int DataToSend,i;
	
	if((DataFile = fopen(DataFilePath,"r"))==NULL)
	{
		fprintf(stderr, "Error in open file %s : %s\n",DataFilePath,strerror(errno));
		exit(EXIT_FAILURE);
	}
	
	fscanf(DataFile,"%d,",RootData);
	
	for(i=1; i<Size; i++)
	{
		fscanf(DataFile,"%d,",&DataToSend);
		MPI_Send(&DataToSend, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
	}
	
	fclose(DataFile);
}

void calcMaxHypercubeNetwork(int Rank, int HypercubeDimension,int *Neightbors, int LocalData, int *Max)
{
	int i,ReceivedData;
	*Max=LocalData;

	for(i=0; i<HypercubeDimension; i++)
	{
		MPI_Send(Max, 1, MPI_INT, Neightbors[i], 0, MPI_COMM_WORLD);
		MPI_Recv(&ReceivedData, 1, MPI_INT, Neightbors[i], MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		
		if(ReceivedData>*Max)
		{
			*Max=ReceivedData;
		}
	} 
	
}
	
void neightborsHypercube(int HypercubeDimension, int Rank, int *Neightbors)
{	
	int i;
	for(i=0; i<HypercubeDimension; i++)
	{
		Neightbors[i] = Rank ^ (int) pow(2,i);
	}
}
	

