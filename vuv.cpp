/*
PRL Project 3
@author Michal Ormos
@email xormos00@stud.fit.vutbr.cz
@date April 2019
*/

#include <stdio.h>
#include <iostream>
#include <mpi.h>
#include <math.h>
#include <vector>
#include <string>
#include <fstream>
#include <time.h>
#include <algorithm>
#include <iterator>
using namespace std; 


#define DEBUG 0
#define TIME 0	
#define TAG 0

/**
	Function for just displaying content of vector
	@param vector of int numbers
	@return display content to stdout
*/
void display_vector(const vector<int> &v)
{
    copy(v.begin(), v.end(),
        ostream_iterator<int>(cout, " "));
}
void display_array(int arr[], int size)
{
	cout << "0:" << arr[0] << "\n";
	for (int i=1; i < size; i++) {
	 	cout << i << ":" << arr[i] << "\n";
	}
	cout << endl;
}
int getElementIndex(int array[], int size, int value)
{ 
	int index = 0;
	while ( index < size && array[index] != value ) ++index;
	return ( index == size ? -1 : index );
}

/**
Algorithm: - level of vertex
1) 	for each e do in parallel
		if e is forward edge then weight(e) = -1
						     else weight(e) = +1
		end if
	end for
2) 	weight = SuffixSums(Etour, weight)
3) 	for each e do in parallel
		if e = (u, v) is forward edge then level(v) = weight(e) + 1
		endif
	endfor
	level(root) = 0
*/
int main(int argc, char *argv[]) {
	const int EOS=-1;
	int processesCount = 0;
	int processId = 0;
	string input;
	// Initialize the MPI
    MPI_Init(&argc, &argv);
    MPI_Status stat;
    // Get number of process
    MPI_Comm_size(MPI_COMM_WORLD,&processesCount); //get number of processes
    MPI_Comm_rank(MPI_COMM_WORLD,&processId); //get actual process id

	input = argv[argc-1];
	int mySuccessor, myPredecessor = 0;
	int isPredParent, isSuccParent = 0;
	int isLeftChildPred, isLeftChildSucc = 0;
	int isRightChildPred, isRightChildSucc = 0;
	int inverseEdge = 0;
	int weight = 0;
	int suffixSum = 0;
	vector<int> adjacencyList;

	if(processesCount == 1) {
		cout << input << ":0\n";
		MPI_Finalize(); 
		return 0;
	}

    /************************ Adjency List Creation **************************/
	// calculate my neighbours and my inverse edge
	if (processId % 2 == 0) {
		mySuccessor = processId/4;
		myPredecessor = processId/2+1;
		inverseEdge = processId + 1;
		weight = -1;
	} else {
		mySuccessor = processId/2+1;
		myPredecessor = processId/4;
		inverseEdge = processId - 1;
		weight = 1;
	}
	// calculate my list ajdency list to lnow, where i am
	if (myPredecessor != 0) {
		adjacencyList.push_back(((myPredecessor-1)*2)+1);
		adjacencyList.push_back((myPredecessor-1)*2);
	}

	int recountMyPredecessor = (((myPredecessor+1)*4)-2);
	if (recountMyPredecessor <= processesCount) {
		adjacencyList.push_back(myPredecessor*4);
		adjacencyList.push_back((myPredecessor*4)+1);
	}

	recountMyPredecessor = ((myPredecessor+1)*4);
	if (recountMyPredecessor <= processesCount) {
		adjacencyList.push_back((myPredecessor*4)+2);
		adjacencyList.push_back((myPredecessor*4)+3);
	}	

	// get next Etour edge for this edge
	int nextEdgeEtour = 0;
	for (int i = 0; i < adjacencyList.size(); i = i + 2) {
		if (adjacencyList[i] == inverseEdge ) {
			if ( i >= (adjacencyList.size()-2)) {
				nextEdgeEtour = adjacencyList[0];
			} else {
				nextEdgeEtour = adjacencyList[i+2];
			}
			break;
		}
	}
	if (DEBUG) {
		cout << "\tprocesID: " << processId << " mySuccessor: " << mySuccessor << " myPredecessor: " << myPredecessor << " weight: " << weight << " inverse edge: " << inverseEdge << " next etour: " << nextEdgeEtour << "\n";
		display_vector(adjacencyList);
	}

    /**************************** Euler's Path *******************************/
	// create Euler's Path
	int eTour[processesCount];
	eTour[processId] = nextEdgeEtour;

	MPI_Allgather(&eTour[processId], 1, MPI_INT, eTour, 1, MPI_INT, MPI_COMM_WORLD);
	// wait for all processors to get here
	MPI_Barrier(MPI_COMM_WORLD);
	if (DEBUG && processId == 0) {
		display_array(eTour, (sizeof(eTour)/sizeof(*eTour)));
	}

	/************************ Suffix Sum ************************ SEQUENTIAL */
	int lastEdgeIndex = getElementIndex(eTour, processesCount, getElementIndex(eTour, processesCount, 0));
	int suffixFromSender;
	if(getElementIndex(eTour, processesCount, processId) == eTour[lastEdgeIndex]){ // root element
		MPI_Recv(&suffixFromSender, 1, MPI_INT, eTour[processId], 0, MPI_COMM_WORLD, &stat);
		suffixSum = suffixFromSender + weight;
	} else if(processId == eTour[lastEdgeIndex]){ // last element of eTour
		suffixSum = weight;
		MPI_Send(&suffixSum, 1, MPI_INT, getElementIndex(eTour, processesCount, processId), 0, MPI_COMM_WORLD);
	} else {
		MPI_Recv(&suffixFromSender, 1, MPI_INT, eTour[processId], 0, MPI_COMM_WORLD, &stat);
		suffixSum = suffixFromSender + weight;
		MPI_Send(&suffixSum, 1, MPI_INT, getElementIndex(eTour, processesCount, processId), 0, MPI_COMM_WORLD);
	}
	suffixSum = suffixSum * weight;
	if (DEBUG) {
		cout << "My ID: " << processId << " Suffix Sum: " << suffixSum << "\n";
	}
	MPI_Barrier(MPI_COMM_WORLD); // wait for all processors to get here

	/****************** Gather Suffixes and print the tree *******************/
	if (processId == 0) {
		int tmpSuffix = 0;
		int results[processesCount];
		int i = 0;
		results[0] = suffixSum;
		for (i = 1; i < processesCount; i++) {
			MPI_Recv(&tmpSuffix, 1, MPI_INT, i, TAG, MPI_COMM_WORLD, &stat);
			results[i] = tmpSuffix;
		}
		// print nodes by value (first the greatest ones)
		int z = 0;
		cout << input[0] << ":0";
		for (i=1; i<input.length(); i++) {
			if (i == 1) {
				z = i;
			} else if ( i == 2 ) {
				z = i+1;
			} else {
				z = i*2-1;
			}
			cout << "," << input[i] << ":" << results[z];
		}
		cout << "\n";
	} else {
		MPI_Send(&suffixSum, 1, MPI_INT, 0, TAG, MPI_COMM_WORLD);
		MPI_Send(&myPredecessor, 1, MPI_INT, 0, TAG, MPI_COMM_WORLD);
	}
    // releive all alocated processes
    MPI_Finalize();
    // finish programm
    return 0;    
}	