/**
 *
 * mpmfunctions.h
 *
 */

#ifndef _MPMFUNCTIONS_H_
#define _MPMFUNCTIONS_H_

#include<mpmstructs.h>

void clearGridValues(GridData *grid_data, Node** grid);
void wpClearGridValues(GridData *grid_data, Node** grid);
void mapToGrid(GridData *grid_data, Node** grid, Particle* particles);
void momentumToVelocityOnGrid(GridData *grid_data, Node** grid);
void updateStressAndStrain(GridData *grid_data, Node** grid, Particle* particles, double dt);
void computeForces(GridData *grid_data, Node **grid, Particle *particles);
void computeAcceleration(GridData *grid_data, Node **grid);
void updateNodalValues(GridData *grid_data, Node **grid, double dt);
void updateParticlePosition(GridData *grid_data, Node **grid, Particle *particles, double dt);
void pUpdateParticlePosition(GridData *grid_data, Node **grid, Particle *particles,
                             Particle *post_particles, Particle **shalo_parts, double dt);
#endif
