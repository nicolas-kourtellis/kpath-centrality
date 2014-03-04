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

// This main program allows the parsing of user input values to call the kpath centrality.
// Various checks are performed on the user input but kind handling is expected.

#include "readgml.h"
#include "betweenness.h"
#include "kpath.h"

void user_inputs(f64 &alpha, ui64 &plength, NETWORK *network, char *argv[]); 
void Print (f64 CB[], f64 NOV [ ], NETWORK *network, FILE *OutPtr ); 

/* 
 * Main function 
 */
int main (int argc, char *argv[]) {
	
	ui64 i, plength = 0, nvertices = 0, nedges = 0;
	f64 alpha, brandes_time = 0, kpath_time = 0;
	f64 *CB, *NOV;
	FILE *InPtr, *OutPtr; 
	NETWORK network; 
	
	// Check command-line arguments 
	if (argc < 5) { 
		cout << "Usage: ./kpath_centrality <infile.gml> <outfile.csv> <k-path alpha> <k-path length> " << endl;
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
	NOV = (f64 *) calloc (nvertices, sizeof(f64)); 
	
	if ((CB == NULL) || (NOV == NULL)) { 
		cout << "Allocating memory failed" << endl; 
		free(CB); 
		free(NOV); 
		free_network(&network); 
		return 0; 
	}
	
	// Set the centralities of all vertices to zero 
	for (i=0; i < nvertices; i++) {
		CB[i] = 0; 
		NOV[i] = 0; 
	}
	
	// User inputs - alpha, plength, and epsilon and c-thr and sup
	user_inputs(alpha, plength, &network, argv); 
	
	// Open the output file 
	OutPtr = fopen(argv[2], "w"); 
	if (OutPtr == NULL) {
		cout << "Unable to open the output file" << endl; 
		return 0; 
	}
	
	//Compute and print betweenness centrality
	BrandesAlgorithm(CB, &network, brandes_time);

	//Compute and print k-path centrality
	kpathcentrality(NOV, &network, alpha, plength, kpath_time); 
	
	//Write file header
	fprintf(OutPtr, "Input file name:," );
	fprintf(OutPtr, "%s", argv[1] );
	fprintf(OutPtr, ",nvertices:,%ld,nedges:,%ld,directed:,%ld,", (ui64) network.nvertices, (ui64) network.nedges, (ui64) network.directed );
	fprintf(OutPtr, "max_weight:,%f,min_weight:,%f\n", network.MAX_Weight, network.MIN_Weight );
	fprintf(OutPtr, "alpha:,%f,plength:%ld\n", alpha, plength);
	fprintf(OutPtr, "Brandes time:,%f,kpath time:,%f\n", brandes_time, kpath_time); 
	fprintf(OutPtr, "Vertex,Brandes,KPath\n"); 
	
	Print(CB, NOV, &network, OutPtr);
	
	//Free memory
	free(CB);
	free(NOV);
	free_network(&network);
	fclose(OutPtr); 
	cout << "Done" << endl; 
	return 0;
	
} // End Main 

/* 
 * User inputs - alpha, plength 
 */ 
void user_inputs(f64 &alpha, ui64 &plength, NETWORK *network, char *argv[]) {
	
	ui64 numV, numE;
	numV = network->nvertices; 
	numE = network->nedges; 
	
	// k-path alpha
	alpha = atof(argv[3]); 
	if ((alpha < -0.5) || (alpha > 0.5)) {
		alpha = 0;
		cout << "Using the default value of alpha = " << alpha << endl; 
	}
	
	// k-path length
	plength = atol(argv[4]); 
	if ((plength <= 0) || (plength > numV)) {
		plength = (ui64) log((f64) (numV + numE));
		if (log ((f64) (numV + numE)) - plength >= 0.5) 
			plength += 1; 
		cout << "Using the default value of k-path length = " << plength << endl; 
	}
	
	cout << "alpha = " << alpha << " and k-path length = " << plength << endl;

	return;
} // End user_inputs


/*
 * Print All centrality values
 */ 
void Print (f64 CB[], f64 NOV [ ], NETWORK *network, FILE *OutPtr) {

	ui64 i; 
	
	for ( i = 0; i < (ui64) network->nvertices; i++ ) {
		fprintf(OutPtr, "%lu,%f,%f\n", i, CB[i], NOV[i]); 
	}
	
	return; 
}
