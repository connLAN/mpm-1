/**
 *
 * commfunctions.h
 *
 */

#ifndef _COMMFUNCTIONS_H_
#define _COMMFUNCTIONS_H_

#include<mpmstructs.h>
#include<mpi.h>

void decomposeGrid(int rank, int *num_proc, int *domain_nc, int *num_cells, double *domain_box, double *box,
                   int *halo_cells, GridData *grid_data);
void decomposeMaterial(double *box, double *m_dom_box, double *m_box, double *ipart_dim, int *num_part);

void allocateHaloParticles(GridData *grid_data, int sfactor, int *pp_cell,
                           Particle **rhalo_parts, Particle **shalo_parts);

void freeHaloParticles(GridData *grid_data, Particle **rhalo_parts, Particle **shalo_parts);

void gatherHaloParticles(GridData *grid_data, Particle *particles, Particle **rhalo_parts, Particle **shalo_parts);
void wpGatherHaloParticles(GridData **grid_data, Particle ***particles, Particle ***new_particles,
                           Particle ***halo_particles);

void sendRecvParticles(GridData *grid_data, Particle **rhalo_parts, Particle **shalo_parts);
void sendRecvParticles1(GridData *grid_data, Particle **rhalo_parts, Particle **shalo_parts,
                        MPI_Request **halo_req, MPI_Status **halo_stat);
void sendRecvParticles2(GridData *grid_data, MPI_Request **halo_req, MPI_Status **halo_stat);
void updateParticleList(GridData *grid_data, Particle *particles, Particle **rhalo_parts);
Particle* gatherParticles(int rank, int size, int array_size, Particle *particles, Particle *post_particles);

void wpRelocateParticles(GridData **grid, Particle *particles, Particle ***particles_new);

void getIndex(int rank, int procX, int procY, int* gridI, int* gridJ);
int getRank(int procX, int procY, int gridI, int gridJ);

#endif
