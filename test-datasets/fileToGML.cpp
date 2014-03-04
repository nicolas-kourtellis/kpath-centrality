#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
using namespace std;

//Main
int main(int argc, char *argv[]) {
	int i, j, k;
	int nvertices = 0, nedges = 0;
	vector < int > nodes;
	FILE *ptrRead, *ptrWrite;
	
	// Check command-line arguments 
	if (argc < 3) { 
		printf( "Usage: ./FileToGML <infile.txt> <outfile.gml>\nNote: the infile.txt must have lines with 3 space-separated columns:\nnode1 node2 weight\n" );
		return 0; 
	}

	//Open the read file
	ptrRead = fopen( argv[1], "r" );
	if( ptrRead == NULL ) {
		printf( "Unable to open the data file.\n" );
		return 0;
	}

	//Open the write file
	ptrWrite = fopen( argv[2], "w" );
	if( ptrRead == NULL ) {
		printf( "Unable to open the write file.\n" );
		return 0;
	}
	fprintf( ptrWrite, "GML file created from file " );
	fprintf( ptrWrite, "%s", argv[1] );
	fprintf( ptrWrite, " to gml file\ngraph\n[\n" );

	//Count number of vertices and edges
	while( fscanf( ptrRead, "%i %i %i\n", &i, &j, &k) != EOF ) {
		if( nvertices < i )
			nvertices = i;
		if( nvertices < j )
			nvertices = j;
		nedges++;
	}
	nodes.assign( nvertices + 1, 0 );

	//Write all node id's to the gml file
	rewind( ptrRead );
	while( fscanf( ptrRead, "%i %i %i\n", &i, &j, &k) != EOF ) {
		if( nodes[i] == 0 ) {
			fprintf( ptrWrite, "  node\n  [\n    id %i\n  ]\n", i );
			nodes[i] = 1;
		}
		if( nodes[j] == 0 ) {
			fprintf( ptrWrite, "  node\n  [\n    id %i\n  ]\n", j );
			nodes[j] = 1;
		}
	}

	//Write all edges to the gml file
	rewind( ptrRead );
	while( fscanf( ptrRead, "%i %i %i\n", &i, &j, &k) != EOF ) {
		fprintf( ptrWrite, "  edge\n  [\n    source %i\n", i );
		fprintf( ptrWrite, "    target %i\n    value %i\n  ]\n", j, k );
	}

	//Close all files and release allocated memory
	fprintf( ptrWrite, "]\n" );
	fclose( ptrRead );
	fclose( ptrWrite );
	nodes.clear( );
	printf( "gml file created successfully\n" );
	return 0;
}

