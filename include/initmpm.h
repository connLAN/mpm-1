/**
 *
 * initmpm.h
 *
 */

#ifndef _INITMPM_H_
#define _INITMPM_H_

#include<stdlib.h>
#include<stdio.h>
#include<mpmstructs.h>

Node** createGrid(GridData *grid_data, double* box, int *num_cells);
void wpCreateGrid(GridData *grid_data, double *box, int *num_cells);
void freeGrid(Node** grid, GridData *grid_data);

Particle* createMaterial(GridData *grid_data, InitialParticleData *ipart_data);
void freeMaterial();

void initializeGrid(GridData *grid_data, Node** grid);
void initializeMaterial(GridData *grid_data, InitialParticleData *ipart_data, Particle *particles);
#endif
