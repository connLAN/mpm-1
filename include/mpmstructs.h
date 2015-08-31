/**
 * mpmstructs.h
 *
 **/

#ifndef _MPMSTRUCTS_H_
#define _MPMSTRUCTS_H_

typedef enum {MPM_FALSE, MPM_TRUE} MpmBool;
typedef enum {BOUNDARY, HALO, INTERIOR, PERIMETER} NodeType;
typedef enum {DIRICHLET, NEUMANN, REFLEXIVE} BoundaryType;

typedef struct{
  int i;
  int j;
}IJIndex;

typedef struct{
  double vel[2];
}BoundaryCondition;

typedef struct{
  double density;
  double bulk;
  double shear;
  double E;
  double poisson;
  double idim[2];
  double velocity[2];
  double box[4];
  int num_particles[2];
  int domain_num_particles[2];
}InitialParticleData;

/**
   * grid boundary and halo cell numberining goes counter clockwise
   *
   *        3
   *       ___
   *    0 |   | 2
   *       ---
   *        1
   **/

typedef struct{
  MpmBool parallel;

  int grid_rank;
	int mem_size;
	int buf_size;
	int halo_index;
  int particle_count;
  int halo_particle_count[8];
  int phalo_size[8];
  int num_nodes[2];
  int halo_cells[4];
  int num_proc[2];
  int rank[8];
  int max_num_nodes[2];

  double gravity[2];
  double cell_dim[2];
  double box[4];  // array index: 0 = x0, 1 = y0, 2 = x1, 3 = y1

  NodeType gridBoundaries[4];
  BoundaryType gridBoundaryType[4];
  IJIndex patch_index;
}GridData;

/**
 * For a particle cell the numbers go
 * counter clockwise with 0 being in the
 * bottom left
 *
 *    3 --- 2
 *    |  p  |
 *    0 --- 1
 *
 * */
typedef struct{
  int particle_count;
  int halo_count;
  int halo_index;
  int halo_cells[4];
  int num_nodes[2];
  MpmBool halt;
  double box[4];
  double mass;
  double vol;
  double position[2];
  double velocity[2];
  double vel_grad[2][2];
  double stress[2][2];
  double strain[2][2];
  double F[2][2];
  double basis[4];
  double dbasisX[4];
  double dbasisY[4];
  double E;
  double poisson;
  IJIndex nodes[4];
  IJIndex gridID;
}Particle;

typedef struct{
  int cellID[2];
  double mass;
  double position[2];
  double int_force[2];
  double body_force[2];
  double momentum[2];
  double velocity[2];
  double acceleration[2];
  int particleIDs[30];
  int num_particles;
  NodeType node_type;
  MpmBool has_mass;
}Node;

#endif
