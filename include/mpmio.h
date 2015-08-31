/**
 *
 * mpmio.h
 *
 **/


#ifndef _MPMIO_H_
#define _MPMIO_H_

#include <stdio.h>
#include <stdlib.h>
#include <mpmstructs.h>

#define ASCII 0
#define BINARY 1
void writeParticlesVtkFile(GridData *grid_data, Node **grid, Particle *particles, const char* filename);
void writeGridVtkFile(GridData *grid_data, Node **grid, const char* filename);
void writeHeader(FILE * outfile, const char * title, int datatype );
void writeParticleData(FILE * outfile, GridData *grid_data, Node **grid, Particle *particles);
void writeGridData(FILE * outfile, GridData *grid_data, Node **grid);

#endif
