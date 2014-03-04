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

#include "kpath.h"

/* 
 * K-Path Centrality for weighted graphs
 */ 
void kpathcentrality_Weighted( f64 NOV[ ], NETWORK *network, f64 alpha, ui64 plength, f64 &time_dif ) {
	
	ui64 i, j, k, x, nloops;
	ui64 nvertices, degree, randL, *Explored;
	f64 randWeight, Weight, TotInvWeight;
	stack <ui64> S;  
	time_t start, end;
	
	//Start time before k-path Centrality Algorithm
	time ( &start );
	
	//Get user inputs and calculate number of loops
	nvertices = (ui64) network->nvertices;
	nloops = (ui64) (2 * plength * plength * pow((f64)nvertices,(1 -(2*alpha))) * log((f64)nvertices) ) + 1;
	
	//Allocate memory
	Explored = ( ui64* ) calloc ( nvertices, sizeof ( ui64 ) );
	if ( Explored == NULL ) {
		cout << "Allocating memory for Explored list failed." << endl;
		return;
	}
	
	// Set all vertices to be unexplored and set NOV to zero 
	for ( j = 0; j < nvertices; j++ ) { 
		Explored[ j ] = 0;
		NOV[j] = 0; 
	}
	
	//Generate a random seed using time
	srand ( time ( NULL ) );
	
	// k-path Centrality Algorithm for weighted graphs
	for ( i = 0; i < nloops; i++ ) {
		
		//pick a random vertex as the source vertex
		x = (ui64) ((((f64) rand())/((f64) RAND_MAX + 1.0))*nvertices);
		
		// Get the degree of the randomly chosen vertex x
		degree = (ui64) network->vertex[ x ].degree;
		
		// As long as the degress is zero, keep randomly choosing the vertex x
		while( degree == 0 ) {
			x = (ui64) ((((f64) rand())/((f64) RAND_MAX + 1.0))*nvertices);
			degree = (ui64) network->vertex[ x ].degree;
		}
		
		Explored[ x ] = 1;
		S.push(x); 
		
		// Pick a random length less or equal to path length
		randL = (ui64) ((((f64) rand())/((f64) RAND_MAX + 1.0))*plength) + 1;

		// Inner loop less or equal to path length
		for ( j = 0; j < randL; j++ ) {
			TotInvWeight = 0;
			
			// Add all edge weights that lead to unexplored vertices
			for ( k = 0; k < degree; k++ ) {
				if ( Explored[ network->vertex[ x ].edge[ k ].target ] == 0 ) {
					Weight = network->vertex[ x ].edge[ k ].weight;
					if ( Weight != 0 )
						TotInvWeight += ( 1 / Weight );
				}
			}

			// If all edges lead to explored vertices break from the inner loop
			if ( TotInvWeight == 0 )
				break;
			
			/* 
			 * Randomly pick an edge out of the remaining unexplored edges with probability 
			 * inversely proportional to its edge weight
			 */
			randWeight = (((f64) rand())/((f64) RAND_MAX + 1.0))*TotInvWeight;
			TotInvWeight = 0;
			for ( k = 0; k < degree; k++ ) {
				if ( Explored[ network->vertex[ x ].edge[ k ].target ] == 0 ) {
					Weight = network->vertex[ x ].edge[ k ].weight;
					if ( Weight != 0 ) {
						TotInvWeight += ( 1 / Weight );
						if ( TotInvWeight > randWeight )
							break;
					}
				}
			}

			// Set the target vertex as the new source vertex
			x = network->vertex[ x ].edge[ k ].target;
			
			// Set the degree to the new vertex degree
			degree = (ui64) network->vertex[ x ].degree;
			 
			//mark the new vertex as explored and increase the number of visites
			Explored[ x ] = 1;
			NOV[ x ] += 1;
			S.push(x); 
			
		} // End For loop for path length
		
		while(!S.empty()) {
			x = S.top(); 
			Explored[x] = 0; 
			S.pop();

			/* if message traversal stops in less than l edges, reset count values to the old ones */
			if(j < randL )
				NOV[ x ] -= 1;
		}
		
	} // End For loop for number of iterations
	
	//End time after k-path Centrality Algorithm and Time difference
	time ( &end );
	time_dif = difftime ( end, start );
	cout << "It took " << time_dif << " seconds to calculate k-path Centrality on a weighted graph" << endl;
	
	//Approximate value
	for ( i = 0; i < nvertices; i++ ) 
		NOV[ i ] = (NOV[ i ] * plength * nvertices ) / nloops;		
	
	//Dealocate memory
	free ( Explored );
	return; 
	
} // End Program

/* 
 * K-Path Centrality for unweighted graphs
 */ 
void kpathcentrality_Unweighted( f64 NOV[ ], NETWORK *network, f64 alpha, ui64 plength, f64 &time_dif ) {
	
	ui64 i, j, k, x, count, randCount, nloops;
	ui64 nvertices, degree, randL, *Explored;
	stack <ui64> S;  
	time_t start, end;
	
	//Start time before k-path Centrality Algorithm
	time ( &start );
	
	//Get user inputs and calculate number of loops
	nvertices = (ui64) network->nvertices;
	nloops = (ui64) (2 * plength * plength * pow((f64)nvertices,(1 -(2*alpha))) * log((f64)nvertices) ) + 1;
	 
	//Allocate memory
	Explored = ( ui64* ) calloc ( nvertices, sizeof ( ui64 ) );
	if ( Explored == NULL ) {
		cout << "Allocating memory for Explored list failed." << endl;
		return;
	}

	// Set all vertices to be unexplored and set NOV to zero 
	for ( j = 0; j < nvertices; j++ ) { 
		Explored[ j ] = 0;
		NOV[j] = 0; 
	}
	
	//Generate a random seed using time
	srand ( time ( NULL ) );

	/* k-path Centrality Algorithm for unweighted graphs */ 
	for ( i = 0; i < nloops; i++ ) {
	
		//pick a random vertex as the source vertex
		x = (ui64) ((((f64) rand())/((f64) RAND_MAX + 1.0))*nvertices);
		
		// Get the degree of the randomly chosen vertex x
		degree = (ui64) network->vertex[ x ].degree;
			
		// As long as the degress is zero, keep randomly choosing the vertex x
		while( degree == 0 ) {
			x = (ui64) ((((f64) rand())/((f64) RAND_MAX + 1.0))*nvertices);
			degree = (ui64) network->vertex[ x ].degree;
		}
		
		Explored[ x ] = 1;
		S.push(x); 
		
		// Pick a random length less or equal to path length
		randL = (ui64) ((((f64) rand())/((f64) RAND_MAX + 1.0))*plength) + 1;

		// Inner loop less or equal to path length
		for ( j = 0; j < randL; j++ ) {
			count = 0;
			
			// Add all edge weights that lead to unexplored vertices
			for ( k = 0; k < degree; k++ )
				if ( Explored[ network->vertex[ x ].edge[ k ].target ] == 0 )
					count++;
			
			// If all edges lead to explored vertices break from the inner loop
			if ( count == 0 )
				break;
			
			// Randomly pick an edge out of the remaining unexplored edges
			randCount = (ui64) ((((f64) rand())/((f64) RAND_MAX + 1.0))*count);
			count = 0;
			for ( k = 0; k < degree; k++ ) {
				if ( Explored[ network->vertex[ x ].edge[ k ].target ] == 0 ) {
					count++;
					if ( count > randCount )
						break;
				}
			}
			
			// Set the target vertex as the new source vertex
			x = network->vertex[ x ].edge[ k ].target;
					
			// Set the degree to the new vertex degree
			degree = (ui64) network->vertex[ x ].degree;
			
			//mark the new vertex as explored and increase the number of visites
			Explored[ x ] = 1;
			NOV[ x ] += 1;
			S.push(x); 
			
		} // End For loop for path length
		
		while(!S.empty()) {
			x = S.top(); 
			Explored[x] = 0; 
			S.pop(); 

			/* if message traversal stops in less than l edges, reset count values to the old ones */
			if(j < randL )
				NOV[ x ] -= 1;
		}		
		
	} // End For loop for number of iterations

	//End time after k-path Centrality Algorithm and Time difference
	time ( &end );
	time_dif = difftime ( end, start );
	cout << "It took " << time_dif << " seconds to calculate k-path Centrality on a unweighted graph" << endl;

	//Approximate value
	for ( i = 0; i < nvertices; i++ ) 
		 NOV[ i ] = (NOV[ i ] * plength * nvertices ) / nloops;

	//Dealocate memory
	free ( Explored );
	
	return; 

} // End Program


/* 
 * K-Path Centrality - Choose between weighted or unweighted graphs
 */ 
void kpathcentrality( f64 NOV[ ], NETWORK *network, f64 alpha, ui64 plength, f64 &time_dif ) {
	if ( ( network->MAX_Weight != 1 ) || ( network->MIN_Weight != 1 ) )
		kpathcentrality_Weighted ( NOV, network, alpha, plength, time_dif );
	else
		kpathcentrality_Unweighted ( NOV, network, alpha, plength, time_dif );
	return; 
}
