/*

Copyright (c) 2014, Nicolas Kourtellis (extensions)

Copyright (c) 2012, Tharaka Alahakoon, Rahul Tripathi, Nicolas Kourtellis

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

*/

#include "readgml.h"
#include "betweenness.h"

void user_inputs(f64 &epsilon, f64 &c_thr, f64 &sup, NETWORK *network, char *argv[]); 
void Print (f64 CB[], f64 RCB [], f64 ACB [], NETWORK *network, FILE *OutPtr ); 

/* 
 * Main function 
 */
int main (int argc, char *argv[]) {
	
	ui64 i, nvertices = 0, nedges = 0;
	f64 brandes_time = 0, randBrandes_time = 0, AdapSample_time = 0, epsilon, c_thr, sup; 
	f64 *CB, *RCB, *ACB; 
	FILE *InPtr, *OutPtr; 
	NETWORK network; 
	
	// Check command-line arguments 
	if (argc < 6) { 
		cout << "Usage: ./rand-brandes_adap-sampl_centrality <infile.gml> <outfile.csv> ";
		cout << "<epsilon for randomized betweennness> <c-threshold for adaptive sampling> <number of pivots for adaptive sampling>" << endl;
		exit(1);
	}

	// Open the gml file 
	InPtr = fopen(argv[1], "r"); 
	if (InPtr == NULL) {
		cout << "Unable to open the input file" << endl; 
		return 0; 
	}
	
	// Read the gml file and create the network 
	if (read_network(&network, InPtr) != 0) {
		cout << "Error creating the network" << endl; 
		return 0; 
	}
	
	// Close the gml file
	fclose(InPtr); 
	
	// Display the number of vertices and number of edges 
	nvertices = network.nvertices; 
	nedges = network.nedges; 
	cout << "Number of vertices =  " << nvertices << "  and number of edges == " << nedges << endl; 
	cout << "Maximum weight = " << network.MAX_Weight << " and " << "Minimum weight = " << network.MIN_Weight << endl; 
	cout << "Network is directed(1) OR undirected(0) = " << network.directed << endl; 
	
	// Allocate memory for CB == Betweenness Centrality and NOV == Number of Visits for K-Path Centrality
	CB = (f64 *) calloc (nvertices, sizeof(f64)); 
	RCB = (f64 *) calloc (nvertices, sizeof(f64)); 
	ACB = (f64 *) calloc (nvertices, sizeof(f64)); 
	
	if ((CB == NULL) || (RCB == NULL) || (ACB == NULL)) { 
		cout << "Allocating memory failed" << endl; 
		free(CB); 
		free(RCB); 
		free(ACB); 
		free_network(&network); 
		return 0; 
	}
	
	// Set the centralities of all vertices to zero 
	for (i=0; i < nvertices; i++) {
		CB[i] = 0; 
		RCB[i] = 0; 
		ACB[i] = 0; 
	}
	
	// User inputs - epsilon and c-thr and sup
	user_inputs(epsilon, c_thr, sup, &network, argv); 
	
	// Open the output file 
	OutPtr = fopen(argv[2], "w"); 
	if (OutPtr == NULL) {
		cout << "Unable to open the output file" << endl; 
		return 0; 
	}
	
	//Compute and print betweenness centrality
	BrandesAlgorithm(CB, &network, brandes_time);

	//Compute and print randomized approximate betweenness centrality
	Rand_BrandesAlgorithm(RCB, &network, epsilon, randBrandes_time);
	
	//Compute and print Adaptive randomized sampling algorithm for betweenness centrality
	Adaptive_Sampling_Algorithm(ACB, &network, c_thr, sup, AdapSample_time);
	
	//Write file header
	fprintf(OutPtr, "Input file name:," );
	fprintf(OutPtr, "%s", argv[1] );
	fprintf(OutPtr, ",nvertices:,%ld,nedges:,%ld,directed:,%ld,", (ui64) network.nvertices, (ui64) network.nedges, (ui64) network.directed );
	fprintf(OutPtr, "max_weight:,%f,min_weight:,%f\n", network.MAX_Weight, network.MIN_Weight );
	fprintf(OutPtr, "epsilon:,%f,c-threshold:,%f,pivots:,%f\n", epsilon, c_thr, sup);
	fprintf(OutPtr, "Brandes time:,%f,RandBrandes time:,%f,AdapSample time:,%f\n", brandes_time, randBrandes_time, AdapSample_time ); 
	fprintf(OutPtr, "Vertex,Brandes,RandBrandes,AdaptiveSample\n"); 
	
	Print(CB, RCB, ACB, &network, OutPtr);  
	
	//Free memory
	free(CB);
	free(RCB); 
	free(ACB); 
	free_network(&network);
	fclose(OutPtr); 
	cout << "Done" << endl; 
	return 0;
	
} // End Main 

/* 
 * User inputs - alpha, plength, and epsilon 
 */ 
void user_inputs(f64 &epsilon, f64 &c_thr, f64 &sup, NETWORK *network, char *argv[]) {
	
	ui64 numV, numE;
	numV = network->nvertices; 
	numE = network->nedges; 
	
	// epsilon for randomized betweenness
	epsilon = atof(argv[3]);
	if ((epsilon < 0) || (epsilon > 1)) {
		epsilon = 0.01;
		cout << "Using the default value of epsilon = " << epsilon << endl; 
	}
	
	// c-threshold for adaptive sampling
	c_thr = atof(argv[4]); 
	if (c_thr < 2) {
		c_thr = 5;
		cout << "Using the default value of c-threshold = " << c_thr << endl; 
	}

	// sup for number of pivots in adaptive sampling
	sup = atof(argv[5]);
	if (sup < 20) {
		sup = 20;
		cout << "Using the default value of sup = " << sup << endl; 
	}

	cout << "epsilon = " << epsilon << " and c-threshold = " << c_thr << " and sup =" << sup << endl;

	return;
} // End user_inputs


/*
 * Print All centrality values
 */ 
void Print (f64 CB[], f64 RCB [], f64 ACB [], NETWORK *network, FILE *OutPtr) {

	ui64 i; 
	
	for ( i = 0; i < (ui64) network->nvertices; i++ ) {
		fprintf(OutPtr, "%lu,%f,%f,%f\n", i, CB[i], RCB[i], ACB[i]); 
	}
	return;
}
