// Header file for the parser for GML files
//
// Mark Newman  11 AUG 06

#ifndef _READGML_H
#define _READGML_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "network.h"

int read_network(NETWORK *network, FILE *stream);
void free_network(NETWORK *network);

#endif
