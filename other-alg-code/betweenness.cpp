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

#include "betweenness.h"

/* 
 * Brandes Algorithm for weighted graphs
 */ 
void BrandesAlgorithm_Weighted(f64 CB[], NETWORK *network, f64 &time_dif) { 
	
	ui64 i, j, u, v;
	ui64 nvertices = (ui64) network->nvertices;	// The number of vertices in the network
	f64 u_distance, v_distance, edgeWeight;		// Variables to store distance estimates or edge weights
	
	time_t start, end;							// Time variables
	
	vector<ui64> sigma;							// sigma is the number of shortest paths
	vector<f64> delta;							// A vector storing dependency of the source vertex on all other vertices
	vector< vector <ui64> > PredList;			// A list of predecessors of all vertices 
	
	stack <ui64> S;								// A stack containing vertices in the order found by Dijkstra's Algorithm
	
	FibHeap PQueue;								// A priority queue storing vertices
	FibHeapNode nodeTemp;						// A particular node stored in the priority queue
	FibHeapNode *nodePtr;						// Pointer to a vertex element stored in the priority queue
	vector<FibHeapNode *> nodeVector;			// A vector of all priority queue elements 
	
	// Set the start time of Brandes' Algorithm
	time(&start); 

	nodeVector.assign ( nvertices, NULL );
	for (i=0; i < nvertices; i++) {
		nodeVector[i] = new FibHeapNode(); 
		nodeVector[i]->Set_vertexPosition(i); 
		//Set all Nodes distance to unsigned long max ULONG_MAX that is assumed to be infinity
		nodeVector[i]->Set_key(ULONG_MAX); 
	}
	
	// Compute Betweenness Centrality for every vertex i
	for (i=0; i < nvertices; i++) {
		
		/* Initialize */ 
		PredList.assign(nvertices, vector <ui64> (0, 0)); 
		sigma.assign(nvertices, 0); 
		sigma[i] = 1; 
		delta.assign(nvertices, 0); 
		
		nodeVector[i]->Set_key(0); 
		PQueue.Insert(nodeVector[i]); 
		
		// While the priority queue is nonempty 
		while (PQueue.GetNumNodes() != 0) {
			// Get the element in the priority queue with the minimum key 
			nodePtr = PQueue.ExtractMin(); 
			// Get the vertex corresponding to the queue element with the minimum key
			u = nodePtr->Get_vertexPosition(); 
			// Push u onto the stack S. Needed later for betweenness computation
			S.push(u); 
			// Shortest path distance from source i to vertex u
			u_distance = nodeVector[u]->Get_key(); 
			// Iterate over all the neighbors of u 
			for (j=0; j < (ui64) network->vertex[u].degree; j++) { 
				// Get the neighbor v of vertex u
				v = (ui64) network->vertex[u].edge[j].target; 
				// Get the weight of the edge (u,v) 
				edgeWeight = (f64) network->vertex[u].edge[j].weight; 
				// If v's shortest path distance estimate has not been set yet, then 
				// set the distance estimate of v and store v in the priority queue
				if (nodeVector[v]->Get_key() == ULONG_MAX) {
					nodeVector[v]->Set_key(u_distance + edgeWeight); 
					PQueue.Insert(nodeVector[v]); 
				}
				// Get the current shortest path distance estimate of v
				v_distance = nodeVector[v]->Get_key(); 
				
				/* Relax and Count */ 
				if (v_distance == u_distance + edgeWeight) { 
					sigma[v] += sigma[u]; 
					PredList[v].push_back(u); 
				}
				if (v_distance > u_distance + edgeWeight) {
					sigma[v] = sigma[u]; 
					PredList[v].clear(); 
					PredList[v].push_back(u); 
					nodeTemp.Set_vertexPosition(v); 
					nodeTemp.Set_key(u_distance + edgeWeight); 
					if (PQueue.DecreaseKey(nodeVector[v], nodeTemp) != 0) 
						cout << "Error decreasing the node key" << endl; 
				}
				
			} // End For 
			
		} // End While
					

		/* Accumulation */ 
		while (!S.empty()) { 
			u = S.top(); 
			S.pop(); 
			for (j=0; j < PredList[u].size(); j++) {
				delta[PredList[u][j]] += ((f64) sigma[PredList[u][j]]/sigma[u]) * (1+delta[u]); 
			}
			if (u != i) 
				CB[u] += delta[u]; 
		}
		
		// Clear data for the next run
		PredList.clear(); 
		sigma.clear(); 
		delta.clear(); 
		for (j=0; j < nvertices; j++) 
			nodeVector[j]->Set_key(ULONG_MAX); 
		
	} // End For 
	
	// End time after Brandes' algorithm and the time difference
	time(&end); 
	time_dif = difftime(end, start); 
	cout << "It took " << time_dif << " seconds to calculate Betweenness Centrality in an weighted graph" << endl; 
	
	// Deallocate memory 
	for (i=0; i < nvertices; i++) 
		delete nodeVector[i]; 
	return; 
	
} // End of BrandesAlgorithm_Weighted 

/* 
 * Brandes Algorithm for unweighted graphs
 */ 
void BrandesAlgorithm_Unweighted(f64 CB[], NETWORK *network, f64 &time_dif) { 
	
	ui64 i, j, u, v;
	ui64 nvertices = (ui64) network->nvertices;	// The number of vertices in the network
	
	time_t start, end;							// Time variables
	
	vector<ui64> d;								// A vector storing shortest distance estimates
	vector<ui64> sigma;							// sigma is the number of shortest paths
	vector<f64> delta;							// A vector storing dependency of the source vertex on all other vertices
	vector< vector <ui64> > PredList;			// A list of predecessors of all vertices 
	
	queue <ui64> Q;								// A priority queue soring vertices
	stack <ui64> S;								// A stack containing vertices in the order found by Dijkstra's Algorithm
	
	// Set the start time of Brandes' Algorithm
	time(&start); 
	
	// Compute Betweenness Centrality for every vertex i
	for (i=0; i < nvertices; i++) {
		
		/* Initialize */ 
		PredList.assign(nvertices, vector <ui64> (0, 0)); 
		d.assign(nvertices, ULONG_MAX); 
		d[i] = 0; 
		sigma.assign(nvertices, 0); 
		sigma[i] = 1; 
		delta.assign(nvertices, 0); 
		Q.push(i); 
		
		// Use Breadth First Search algorithm 
		while (!Q.empty()) {
			// Get the next element in the queue
			u = Q.front(); 
			Q.pop(); 
			// Push u onto the stack S. Needed later for betweenness computation
			S.push(u); 
			// Iterate over all the neighbors of u 
			for (j=0; j < (ui64) network->vertex[u].degree; j++) { 
				// Get the neighbor v of vertex u
				v = (ui64) network->vertex[u].edge[j].target; 
				
				/* Relax and Count */
				if (d[v] == ULONG_MAX) { 
					 d[v] = d[u] + 1; 
					 Q.push(v); 
				} 
				if (d[v] == d[u] + 1) {
					sigma[v] += sigma[u]; 
					PredList[v].push_back(u); 
				}
			} // End For
			
		} // End While 
					
		/* Accumulation */ 
		while (!S.empty()) { 
			u = S.top(); 
			S.pop(); 
			for (j=0; j < PredList[u].size(); j++) {
				delta[PredList[u][j]] += ((f64) sigma[PredList[u][j]]/sigma[u]) * (1+delta[u]); 
			}
			if (u != i) 
				CB[u] += delta[u]; 
		}
		
		// Clear data for the next run
		PredList.clear(); 
		d.clear(); 
		sigma.clear(); 
		delta.clear(); 

	} // End For 
	
	// End time after Brandes' algorithm and the time difference
	time(&end); 
	time_dif = difftime(end, start); 
	cout << "It took " << time_dif << " seconds to calculate Betweenness Centrality in an unweighted graph" << endl; 
	
	return; 
	
} // End of BrandesAlgorithm_Unweighted 

/* 
 * Brandes' Algorithm - Choose between weighted or unweighted graphs 
 */ 
void BrandesAlgorithm(f64 CB[], NETWORK *network, f64 &time_dif) {
	if ((network->MAX_Weight != 1) || (network->MIN_Weight  != 1) ) 
		BrandesAlgorithm_Weighted(CB, network, time_dif); 
	else 
		BrandesAlgorithm_Unweighted(CB, network, time_dif); 
	return;
} 

/* 
 * Randomized Brandes Algorithm for weighted graphs
 */ 
void Rand_BrandesAlgorithm_Weighted(f64 RCB[], NETWORK *network, f64 epsilon, f64 &time_dif) { 

	ui64 i, j, u, v, numSample, randvx;
	ui64 nvertices = (ui64) network->nvertices;	// The number of vertices in the network
	f64 u_distance, v_distance, edgeWeight;		// Variables to store distance estimates or edge weights
	
	time_t start, end;							// Time variables
	
	vector<ui64> sigma;							// sigma is the number of shortest paths
	vector<f64> delta;							// A vector storing dependency of the source vertex on all other vertices
	vector< vector <ui64> > PredList;			// A list of predecessors of all vertices 
	
	vector<ui64> SampleVertex; 
	vector<ui64>::iterator it;						// An iterator of vector elements
	
	stack <ui64> S;								// A stack containing vertices in the order found by Dijkstra's Algorithm
	
	FibHeap PQueue;								// A priority queue storing vertices
	FibHeapNode nodeTemp;						// A particular node stored in the priority queue
	FibHeapNode *nodePtr;						// Pointer to a vertex element stored in the priority queue
	vector<FibHeapNode *> nodeVector;			// A vector of all priority queue elements 
	
	// Set the start time of Randomized Brandes' Algorithm
	time(&start); 

	nodeVector.assign ( nvertices, NULL );
	for (i=0; i < nvertices; i++) {
		nodeVector[i] = new FibHeapNode(); 
		nodeVector[i]->Set_vertexPosition(i); 
		//Set all Nodes distance to unsigned long max ULONG_MAX that is assumed to be infinity
		nodeVector[i]->Set_key(ULONG_MAX); 
	}
	
	// Generate random seed 
	srand((unsigned)time(NULL)); 
	
	numSample = (ui64) (2 * (log ((f64) nvertices) / (epsilon * epsilon)));
	
	SampleVertex.resize(numSample); 
	
	for (i=0; i < numSample; i++) {
		
		// Generate a random vertex
		randvx = (ui64) ((((f64) rand())/((f64) RAND_MAX + 1.0))*nvertices); 
		
		// Insert the randomly sampled vertex
		SampleVertex.push_back(randvx); 
	}
	
	// Compute Randomized Betweenness Centrality using sampled vertices
	for (it= SampleVertex.begin(); it < SampleVertex.end(); it++) {
		
		i = *it; 
	
		/* Initialize */ 
		PredList.assign(nvertices, vector <ui64> (0, 0)); 
		sigma.assign(nvertices, 0); 
		sigma[i] = 1; 
		delta.assign(nvertices, 0); 
		
		nodeVector[i]->Set_key(0); 
		PQueue.Insert(nodeVector[i]); 
		
		// While the priority queue is nonempty 
		while (PQueue.GetNumNodes() != 0) {
			// Get the element in the priority queue with the minimum key 
			nodePtr = PQueue.ExtractMin(); 
			// Get the vertex corresponding to the queue element with the minimum key
			u = nodePtr->Get_vertexPosition(); 
			// Push u onto the stack S. Needed later for betweenness computation
			S.push(u); 
			// Shortest path distance from source i to vertex u
			u_distance = nodeVector[u]->Get_key(); 
			// Iterate over all the neighbors of u 
			for (j=0; j < (ui64) network->vertex[u].degree; j++) { 
				// Get the neighbor v of vertex u
				v = (ui64) network->vertex[u].edge[j].target; 
				// Get the weight of the edge (u,v) 
				edgeWeight = (f64) network->vertex[u].edge[j].weight; 
				// If v's shortest path distance estimate has not been set yet, then 
				// set the distance estimate of v and store v in the priority queue
				if (nodeVector[v]->Get_key() == ULONG_MAX) {
					nodeVector[v]->Set_key(u_distance + edgeWeight); 
					PQueue.Insert(nodeVector[v]); 
				}
				// Get the current shortest path distance estimate of v
				v_distance = nodeVector[v]->Get_key(); 
				
				/* Relax and Count */ 
				if (v_distance == u_distance + edgeWeight) { 
					sigma[v] += sigma[u]; 
					PredList[v].push_back(u); 
				}
				if (v_distance > u_distance + edgeWeight) {
					sigma[v] = sigma[u]; 
					PredList[v].clear(); 
					PredList[v].push_back(u); 
					nodeTemp.Set_vertexPosition(v); 
					nodeTemp.Set_key(u_distance + edgeWeight); 
					if (PQueue.DecreaseKey(nodeVector[v], nodeTemp) != 0) 
						cout << "Error decreasing the node key" << endl; 
				}
				
			} // End For 
			
		} // End While
		
		
		/* Accumulation */ 
		while (!S.empty()) { 
			u = S.top(); 
			S.pop(); 
			for (j=0; j < PredList[u].size(); j++) {
				delta[PredList[u][j]] += ((f64) sigma[PredList[u][j]]/sigma[u]) * (1+delta[u]); 
			}
			if (u != i) 
				RCB[u] += delta[u]; 
		}
		
		// Clear data for the next run
		PredList.clear(); 
		sigma.clear(); 
		delta.clear(); 
		for (j=0; j < nvertices; j++) 
			nodeVector[j]->Set_key(ULONG_MAX); 
		
	} // End For 
	
	for (i=0; i < nvertices; i++) 
		RCB[i] = nvertices * (RCB[i]/numSample); 
	
	// End time after Brandes' algorithm and the time difference
	time(&end); 
	time_dif = difftime(end, start); 
	cout << "It took " << time_dif << " seconds to calculate Randomized Approximate Betweenness Centrality in a weighted graph" << endl; 
	
	// Deallocate memory 
	for (i=0; i < nvertices; i++) 
		delete nodeVector[i]; 
	return; 
	
} // End of Rand_BrandesAlgorithm_Weighted 


/* 
 * Randomized Brandes Algorithm for unweighted graphs
 */ 
void Rand_BrandesAlgorithm_Unweighted(f64 RCB[], NETWORK *network, f64 epsilon, f64 &time_dif) { 
	
	ui64 i, j, u, v, numSample, randvx;
	ui64 nvertices = (ui64) network->nvertices;	// The number of vertices in the network
	
	time_t start, end;							// Time variables
	
	vector<ui64> d;								// A vector storing shortest distance estimates
	vector<ui64> sigma;							// sigma is the number of shortest paths
	vector<f64> delta;							// A vector storing dependency of the source vertex on all other vertices
	vector< vector <ui64> > PredList;			// A list of predecessors of all vertices 
	
	vector<ui64> SampleVertex; 
	vector<ui64>::iterator it;						// An iterator of vector elements
	
	queue <ui64> Q;								// A priority queue soring vertices
	stack <ui64> S;								// A stack containing vertices in the order found by Dijkstra's Algorithm
	
	// Set the start time of Brandes' Algorithm
	time(&start); 
	
	// Generate random seed 
	srand((unsigned)time(NULL)); 
	
	numSample = (ui64) (2 * (log ((f64) nvertices) / (epsilon * epsilon))); 
	
	SampleVertex.resize(numSample); 
	
	for (i=0; i < numSample; i++) {
		
		// Generate a random vertex
		randvx = (ui64) ((((f64) rand())/((f64) RAND_MAX + 1.0))*nvertices); 
		
		// Insert the randomly sampled vertex
		SampleVertex.push_back(randvx); 
	}
	
	// Compute Randomized Betweenness Centrality using sampled vertices
	for (it= SampleVertex.begin(); it < SampleVertex.end(); it++) {
		
		i = *it; 
		
		/* Initialize */ 
		PredList.assign(nvertices, vector <ui64> (0, 0)); 
		d.assign(nvertices, ULONG_MAX); 
		d[i] = 0; 
		sigma.assign(nvertices, 0); 
		sigma[i] = 1; 
		delta.assign(nvertices, 0); 
		Q.push(i); 
		
		// Use Breadth First Search algorithm 
		while (!Q.empty()) {
			// Get the next element in the queue
			u = Q.front(); 
			Q.pop(); 
			// Push u onto the stack S. Needed later for betweenness computation
			S.push(u); 
			// Iterate over all the neighbors of u 
			for (j=0; j < (ui64) network->vertex[u].degree; j++) { 
				// Get the neighbor v of vertex u
				v = (ui64) network->vertex[u].edge[j].target; 
				
				/* Relax and Count */
				if (d[v] == ULONG_MAX) { 
					d[v] = d[u] + 1; 
					Q.push(v); 
				} 
				if (d[v] == d[u] + 1) {
					sigma[v] += sigma[u]; 
					PredList[v].push_back(u); 
				}
			} // End For
			
		} // End While 
		
		/* Accumulation */ 
		while (!S.empty()) { 
			u = S.top(); 
			S.pop(); 
			for (j=0; j < PredList[u].size(); j++) {
				delta[PredList[u][j]] += ((f64) sigma[PredList[u][j]]/sigma[u]) * (1+delta[u]); 
			}
			if (u != i) 
				RCB[u] += delta[u]; 
		}
		
		// Clear data for the next run
		PredList.clear(); 
		d.clear(); 
		sigma.clear(); 
		delta.clear(); 
		
	} // End For 
	
	for (i=0; i < nvertices; i++) 
		RCB[i] = nvertices * (RCB[i]/numSample); 
	
	// End time after Brandes' algorithm and the time difference
	time(&end); 
	time_dif = difftime(end, start); 
	cout << "It took " << time_dif << " seconds to calculate Randomized Approximate Betweenness Centrality in an unweighted graph" << endl; 	
	return; 
	
} // End of Rand_BrandesAlgorithm_Unweighted 

/* 
 * Randomized Brandes' Algorithm - Choose between weighted or unweighted graphs 
 */ 
void Rand_BrandesAlgorithm(f64 RCB[], NETWORK *network, f64 epsilon, f64 &time_dif) {
	if ((network->MAX_Weight != 1) || (network->MIN_Weight  != 1) ) 
		Rand_BrandesAlgorithm_Weighted(RCB, network, epsilon, time_dif); 
	else 
		Rand_BrandesAlgorithm_Unweighted(RCB, network, epsilon, time_dif); 
	return;
} 


/* 	 
 * Adaptive Randomized Sampling Algorithm for Weighted graphs. The cut-off on the number of samples is n/20. 
*/ 
void Adaptive_Sampling_Weighted(f64 ACB[], NETWORK *network, f64 c_thr, f64 sup, f64 &time_dif) { 
	
	ui64 i, j, u, v, numSample, randvx;
	ui64 nvertices = (ui64) network->nvertices;	// The number of vertices in the network
	ui64 count = 0; 
	f64 u_distance, v_distance, edgeWeight;		// Variables to store distance estimates or edge weights
	
	time_t start, end;							// Time variables
	
	vector<ui64> sigma;							// sigma is the number of shortest paths
	vector<f64> delta;							// A vector storing dependency of the source vertex on all other vertices
	vector< vector <ui64> > PredList;			// A list of predecessors of all vertices 
	
	vector<ui64> SampleVertex; 
	vector<ui64>::iterator it;					// An iterator of vector elements
	vector<bool> Flag; 
	
	stack <ui64> S;								// A stack containing vertices in the order found by Dijkstra's Algorithm
	
	FibHeap PQueue;								// A priority queue storing vertices
	FibHeapNode nodeTemp;						// A particular node stored in the priority queue
	FibHeapNode *nodePtr;						// Pointer to a vertex element stored in the priority queue
	vector<FibHeapNode *> nodeVector;			// A vector of all priority queue elements 
	
	// Set the start time of Randomized Brandes' Algorithm
	time(&start); 

	nodeVector.assign ( nvertices, NULL );
	for (i=0; i < nvertices; i++) {
		nodeVector[i] = new FibHeapNode(); 
		nodeVector[i]->Set_vertexPosition(i); 
		//Set all Nodes distance to unsigned long max ULONG_MAX that is assumed to be infinity
		nodeVector[i]->Set_key(ULONG_MAX); 
	}
	
	// Generate random seed 
	srand((unsigned)time(NULL)); 
	
	numSample = (ui64) (nvertices/sup); 
	
	if (numSample < 1) 
		numSample = nvertices; 
	
	SampleVertex.resize(numSample); 
	Flag.assign(nvertices, false); 
	
	for (i=0; i < numSample; i++) {
		
		// Generate a random vertex
		randvx = (ui64) ((((f64) rand())/((f64) RAND_MAX + 1.0))*nvertices); 
		
		// Insert the randomly sampled vertex
		SampleVertex.push_back(randvx); 
	}
	
	// Compute Randomized Betweenness Centrality using sampled vertices
	for (it= SampleVertex.begin(); it < SampleVertex.end(); it++) {
		
		count += 1;
		i = *it; 
		
		/* Initialize */ 
		PredList.assign(nvertices, vector <ui64> (0, 0)); 
		sigma.assign(nvertices, 0); 
		sigma[i] = 1; 
		delta.assign(nvertices, 0); 
		
		nodeVector[i]->Set_key(0); 
		PQueue.Insert(nodeVector[i]); 
		
		// While the priority queue is nonempty 
		while (PQueue.GetNumNodes() != 0) {
			// Get the element in the priority queue with the minimum key 
			nodePtr = PQueue.ExtractMin(); 
			// Get the vertex corresponding to the queue element with the minimum key
			u = nodePtr->Get_vertexPosition(); 
			// Push u onto the stack S. Needed later for betweenness computation
			S.push(u); 
			// Shortest path distance from source i to vertex u
			u_distance = nodeVector[u]->Get_key(); 
			// Iterate over all the neighbors of u 
			for (j=0; j < (ui64) network->vertex[u].degree; j++) { 
				// Get the neighbor v of vertex u
				v = (ui64) network->vertex[u].edge[j].target; 
				// Get the weight of the edge (u,v) 
				edgeWeight = (f64) network->vertex[u].edge[j].weight; 
				// If v's shortest path distance estimate has not been set yet, then 
				// set the distance estimate of v and store v in the priority queue
				if (nodeVector[v]->Get_key() == ULONG_MAX) {
					nodeVector[v]->Set_key(u_distance + edgeWeight); 
					PQueue.Insert(nodeVector[v]); 
				}
				// Get the current shortest path distance estimate of v
				v_distance = nodeVector[v]->Get_key(); 
				
				/* Relax and Count */ 
				if (v_distance == u_distance + edgeWeight) { 
					sigma[v] += sigma[u]; 
					PredList[v].push_back(u); 
				}
				if (v_distance > u_distance + edgeWeight) {
					sigma[v] = sigma[u]; 
					PredList[v].clear(); 
					PredList[v].push_back(u); 
					nodeTemp.Set_vertexPosition(v); 
					nodeTemp.Set_key(u_distance + edgeWeight); 
					if (PQueue.DecreaseKey(nodeVector[v], nodeTemp) != 0) 
						cout << "Error decreasing the node key" << endl; 
				}
				
			} // End For 
			
		} // End While
		
		
		/* Accumulation */ 
		while (!S.empty()) { 
			u = S.top(); 
			S.pop(); 
			for (j=0; j < PredList[u].size(); j++) {
				delta[PredList[u][j]] += ((f64) sigma[PredList[u][j]]/sigma[u]) * (1+delta[u]); 
			}
			
			if ((u != i) && (!Flag[u])) {
				ACB[u] += delta[u]; 
				if (ACB[u] > c_thr * nvertices) {
					ACB[u] = nvertices * (ACB[u]/count);
					Flag[u] = true;
				}
			} // End If 
			
		} // End While
		
		// Clear data for the next run
		PredList.clear(); 
		sigma.clear(); 
		delta.clear(); 
		for (j=0; j < nvertices; j++) 
			nodeVector[j]->Set_key(ULONG_MAX); 
		
	} // End For 
	
	for (i=0; i < nvertices; i++) {
		if (!Flag[i]) {
			ACB[i] = nvertices * (ACB[i]/numSample); 
		}
	}
	
	// End time after Brandes' algorithm and the time difference
	time(&end); 
	time_dif = difftime(end, start); 
	cout << "It took " << time_dif << " seconds to calculate Adaptive Sampling Based Approximate Centrality Values in a weighted graph" << endl;
	
	// Deallocate memory 
	for (i=0; i < nvertices; i++) 
		delete nodeVector[i]; 
	return; 
	
} // End of Adaptive_Sampling_Weighted 


/* 
 * Adaptive Randomized Sampling Algorithm for unweighted graphs. The cut-off on the number of samples is n/sup 
 */ 
void Adaptive_Sampling_Unweighted(f64 ACB[], NETWORK *network, f64 c_thr, f64 sup, f64 &time_dif) { 
	
	ui64 i, j, u, v, numSample, randvx;
	ui64 nvertices = (ui64) network->nvertices;	// The number of vertices in the network
	ui64 count = 0; 
	
	time_t start, end;							// Time variables
	
	vector<ui64> d;								// A vector storing shortest distance estimates
	vector<ui64> sigma;							// sigma is the number of shortest paths
	vector<f64> delta;							// A vector storing dependency of the source vertex on all other vertices
	vector< vector <ui64> > PredList;			// A list of predecessors of all vertices 
	
	vector<ui64> SampleVertex; 
	vector<ui64>::iterator it;						// An iterator of vector elements
	vector<bool> Flag; 
	
	queue <ui64> Q;								// A priority queue soring vertices
	stack <ui64> S;								// A stack containing vertices in the order found by Dijkstra's Algorithm
	
	// Set the start time of Brandes' Algorithm
	time(&start); 
	
	// Generate random seed 
	srand((unsigned)time(NULL)); 
	
	numSample = (ui64) (nvertices/sup); 
	if (numSample < 1) 
		numSample = nvertices; 
	
	SampleVertex.resize(numSample); 
	Flag.assign(nvertices, false); 
	
	for (i=0; i < numSample; i++) {
		
		// Generate a random vertex
		randvx = (ui64) ((((f64) rand())/((f64) RAND_MAX + 1.0))*nvertices);
		
		// Insert the randomly sampled vertex
		SampleVertex.push_back(randvx); 
	}
	
	// Compute Randomized Betweenness Centrality using sampled vertices
	for (it= SampleVertex.begin(); it < SampleVertex.end(); it++) {
		
		count += 1; 
		i = *it; 
		
		/* Initialize */ 
		PredList.assign(nvertices, vector <ui64> (0, 0)); 
		d.assign(nvertices, ULONG_MAX); 
		d[i] = 0; 
		sigma.assign(nvertices, 0); 
		sigma[i] = 1; 
		delta.assign(nvertices, 0); 
		Q.push(i); 
		
		// Use Breadth First Search algorithm 
		while (!Q.empty()) {
			// Get the next element in the queue
			u = Q.front(); 
			Q.pop(); 
			// Push u onto the stack S. Needed later for betweenness computation
			S.push(u); 
			// Iterate over all the neighbors of u 
			for (j=0; j < (ui64) network->vertex[u].degree; j++) { 
				// Get the neighbor v of vertex u
				v = (ui64) network->vertex[u].edge[j].target; 
				
				/* Relax and Count */
				if (d[v] == ULONG_MAX) { 
					d[v] = d[u] + 1; 
					Q.push(v); 
				} 
				if (d[v] == d[u] + 1) {
					sigma[v] += sigma[u]; 
					PredList[v].push_back(u); 
				}
			} // End For
			
		} // End While 
		
		/* Accumulation */ 
		while (!S.empty()) { 
			u = S.top(); 
			S.pop(); 
			for (j=0; j < PredList[u].size(); j++) {
				delta[PredList[u][j]] += ((f64) sigma[PredList[u][j]]/sigma[u]) * (1+delta[u]); 
			}
			
			if ((u != i) && (!Flag[u])) {
				ACB[u] += delta[u]; 
				if (ACB[u] > c_thr * nvertices) {
					ACB[u] = nvertices * (ACB[u]/count);
					Flag[u] = true;
				}
			} // End If 
			
		} // End While
		
		// Clear data for the next run
		PredList.clear(); 
		d.clear(); 
		sigma.clear(); 
		delta.clear(); 
		
	} // End For 
	
	for (i=0; i < nvertices; i++) {
		if (!Flag[i]) {
			ACB[i] = nvertices * (ACB[i]/numSample); 
		}
	}
	
	// End time after Brandes' algorithm and the time difference
	time(&end); 
	time_dif = difftime(end, start); 
	cout << "It took " << time_dif << " seconds to calculate Adaptive Sampling Based Approximate Centrality Values in an unweighted graph" << endl; 
	
	return; 
	
} // End of Adaptive_Sampling_Unweighted 


/* 
 *  Adaptive Sampling Based Randomized Approximation Algorithm - Choose between weighted or unweighted graphs 
 */ 
void Adaptive_Sampling_Algorithm(f64 ACB[], NETWORK *network, f64 c_thr, f64 sup, f64 &time_dif) {
	if ((network->MAX_Weight != 1) || (network->MIN_Weight  != 1) ) 
		Adaptive_Sampling_Weighted(ACB, network, c_thr, sup, time_dif); 
	else 
		Adaptive_Sampling_Unweighted(ACB, network, c_thr, sup, time_dif); 
	return;
} 


	

