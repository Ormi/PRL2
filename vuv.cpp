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

/**
	Function for just displaying content of array
	@param array of int numbers
	@param size of array 
	@return display content to stdout
*/
void display_array(int arr[], int size)
{
	int i = 0;
	cout << "0:" << arr[0] << "\n";
	for (i=1; i < size; i++) {
	 	cout << i << ":" << arr[i] << "\n";
	}
	cout << endl;
}

/**
	Function for just getting index of desired element in array
	@param array of int numbers
	@param size of array 
	@param val of which index we are looking for
	@return desired index or -1 if not found
*/
int get_index_of_element(int arr[], int size, int val)
{
	int index = 0; 
	for (index = 0; index < size; index++) {
		if ( arr[index] == val ) {
			return index;
		}
	}
	return -1;
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
	int processesCount = 0;
	int processId = 0;
	// Initialize the MPI
    MPI_Init(&argc, &argv);
    MPI_Status stat;
    // Get number of process
    MPI_Comm_size(MPI_COMM_WORLD,&processesCount); //get number of processes
    MPI_Comm_rank(MPI_COMM_WORLD,&processId); //get actual process id

    // get input string from script
	string input;
	input = argv[argc-1];

	// prepare variables for tree structure calculating
	int mySuccessor, myPredecessor = 0;
	int isPredParent, isSuccParent = 0;
	int isLeftChildPred, isLeftChildSucc = 0;
	int isRightChildPred, isRightChildSucc = 0;
	int inverseEdge = 0;
	int weight = 0;
	int suffixSum = 0;
	vector<int> adjacencyList;

	// If there is only one item in tree
	// Just print it out
	if(processesCount == 1) {
		cout << input << ":0\n";
		MPI_Finalize(); 
		return 0;
	}

    /************************ Adjency List Creation **************************/
	// calculate my neighbours and my inverse edges
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
	// calculate my ajdency list to know, where am I
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
	// each node has 2 edges, iterate i+2
	// if first item of list of vertex
	// else next existuje
	int nextEdgeEtour, i = 0;
	int adjacencyListSize = adjacencyList.size();
	for (i = 0; i < adjacencyListSize; i = i+2) {
		if (adjacencyList[i] == inverseEdge ) {
			if ( i >= (adjacencyListSize - 2)) {
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
	// Processor are building the eTour array together ( MPI_Allgather )
	// And wait until are processes get there ( MPI_Barrier )
	int eTour[processesCount];
	eTour[processId] = nextEdgeEtour;

	MPI_Allgather(&eTour[processId], 1, MPI_INT, eTour, 1, MPI_INT, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	if (DEBUG && processId == 0) {
		display_array(eTour, (sizeof(eTour)/sizeof(*eTour)));
	}

	/****************************** Suffix Sum ********************************/
	int highestEdgeIndex = get_index_of_element(eTour, processesCount, get_index_of_element(eTour, processesCount, 0));
	int suffixFromSender = 0;
	int elementIndex = get_index_of_element(eTour, processesCount, processId);

	// if root element count for root
	// else if last element of etour
	// else other elements
	if(elementIndex == eTour[highestEdgeIndex]) {
		MPI_Recv(&suffixFromSender, 1, MPI_INT, eTour[processId], 0, MPI_COMM_WORLD, &stat);
		suffixSum = suffixFromSender + weight;
	} else if(processId == eTour[highestEdgeIndex]) {
		suffixSum = weight;
		MPI_Send(&suffixSum, 1, MPI_INT, elementIndex, 0, MPI_COMM_WORLD);
	} else {
		MPI_Recv(&suffixFromSender, 1, MPI_INT, eTour[processId], 0, MPI_COMM_WORLD, &stat);
		suffixSum = suffixFromSender + weight;
		MPI_Send(&suffixSum, 1, MPI_INT, elementIndex, 0, MPI_COMM_WORLD);
	}

	suffixSum = suffixSum * weight;

	if (DEBUG) {
		cout << "My ID: " << processId << " Suffix Sum: " << suffixSum << "\n";
	}

	MPI_Barrier(MPI_COMM_WORLD);

	/****************** Gather Suffixes and print the tree *******************/
	if (processId == 0) {
		int finalSuffix = 0;
		int finalWeights[processesCount];
		int i = 0;
		finalWeights[0] = suffixSum;
		for (i = 1; i < processesCount; i++) {
			MPI_Recv(&finalSuffix, 1, MPI_INT, i, TAG, MPI_COMM_WORLD, &stat);
			finalWeights[i] = finalSuffix;
		}
		// print first node of A:0
		// gather other nodes by their back edges
		int z = 0;
		cout << input[0] << ":0";
		for (i = 1; i<input.length(); i++) {
			if (i == 1) {
				z = i;
			} else if ( i == 2 ) {
				z = i+1;
			} else {
				z = i*2-1;
			}
			cout << "," << input[i] << ":" << finalWeights[z];
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