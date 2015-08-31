/**
 *
 * initmpm.c
 *
 */

#include<initmpm.h>

Node** createGrid(GridData *grid_data, double *box, int *num_cells)
{
  int i;

  Node** grid;

  grid_data->box[0] = box[0];
  grid_data->box[1] = box[1];
  grid_data->box[2] = box[2];
  grid_data->box[3] = box[3];

  if(grid_data->parallel == MPM_TRUE){
    grid_data->num_nodes[0] = grid_data->halo_cells[0] + grid_data->halo_cells[2] + num_cells[0] + 1; 
    grid_data->num_nodes[1] = grid_data->halo_cells[1] + grid_data->halo_cells[3] + num_cells[1] + 1; 
  }else{
    grid_data->num_nodes[0] = num_cells[0] + 1; 
    grid_data->num_nodes[1] = num_cells[1] + 1; 
  }

  grid = (Node**) malloc(grid_data->num_nodes[0] * sizeof(Node*)); 
  for(i = 0; i < grid_data->num_nodes[0]; i++){
    grid[i] = (Node*) malloc(grid_data->num_nodes[1] * sizeof(Node)); 
  }
  return grid;

} // end of function createGrid()

void wpCreateGrid(GridData *grid_data, double *box, int *num_cells)
{
  int i;

  grid_data->box[0] = box[0];
  grid_data->box[1] = box[1];
  grid_data->box[2] = box[2];
  grid_data->box[3] = box[3];

  if(grid_data->parallel == MPM_TRUE){
    grid_data->num_nodes[0] = grid_data->halo_cells[0] + grid_data->halo_cells[2] + num_cells[0] + 1; 
    grid_data->num_nodes[1] = grid_data->halo_cells[1] + grid_data->halo_cells[3] + num_cells[1] + 1; 
  }else{
    grid_data->num_nodes[0] = num_cells[0] + 1; 
    grid_data->num_nodes[1] = num_cells[1] + 1; 
  }
} // end of function createGrid()

void freeGrid(Node** grid, GridData *grid_data)
{
  int i;
  for(i = 0; i < grid_data->num_nodes[0]; i++){
    free(grid[i]);
  }
  free(grid);

} // end of function freeGrid()

Particle* createMaterial(GridData *grid_data, InitialParticleData *ipart_data){
  int array_size = ipart_data->domain_num_particles[0] * ipart_data->domain_num_particles[1];
  grid_data->particle_count = ipart_data->num_particles[0] * ipart_data->num_particles[1];

  Particle* particles = (Particle*) malloc(array_size * sizeof(Particle));

  return particles;
}

void freeMaterial(Particle* particles){
  free(particles);
}

void initializeGrid(GridData *grid_data, Node** grid)
{
  int i,j;
  double bx0 = grid_data->box[0];
  double by0 = grid_data->box[1];
  double hx = grid_data->cell_dim[0];
  double hy = grid_data->cell_dim[1];
  int haloX0 = grid_data->halo_cells[0];
  int haloY0 = grid_data->halo_cells[1];
  int haloX1 = grid_data->halo_cells[2];
  int haloY1 = grid_data->halo_cells[3];
  
  if(grid_data->parallel == MPM_TRUE){
    bx0 = grid_data->box[0] - (haloX0 * hx);
    by0 = grid_data->box[1] - (haloY0 * hy);
  }else{
    bx0 = grid_data->box[0];
    by0 = grid_data->box[1];
  }

  for(i = 0; i < grid_data->num_nodes[0]; i++){
    for(j = 0; j < grid_data->num_nodes[1]; j++){
      grid[i][j].position[0] = bx0 + i * hx;
      grid[i][j].position[1] = by0 + j * hy;
    }
  }

} // end of function initializeGrid()

void initializeMaterial(GridData *grid_data, InitialParticleData *ipart_data, Particle *particles)
{
  int i,j,index;
  int particle_count = ipart_data->num_particles[0] * ipart_data->num_particles[1];
  double phx = ipart_data->idim[0];
  double phy = ipart_data->idim[1];
  double volume = phx * phy;
  double bx0 = ipart_data->box[0] + (phx/2);
  double by0 = ipart_data->box[1] + (phy/2);
  //printf("particle count: %d, phx: %f, phy: %f\n", particle_count, phx, phy);
  index = 0;
  for(i = 0; i < ipart_data->num_particles[0]; i++){
    for(j = 0; j < ipart_data->num_particles[1]; j++){
      particles[index].particle_count = particle_count;
      particles[index].mass = ipart_data->density * volume;
      particles[index].E = ipart_data->E;
      particles[index].poisson = ipart_data->poisson;
      particles[index].vol = volume;
      particles[index].position[0] = bx0 + i * phx;
      particles[index].position[1] = by0 + j * phy;
      particles[index].velocity[0] = ipart_data->velocity[0];
      particles[index].velocity[1] = ipart_data->velocity[1];
      particles[index].stress[0][0] = 0;
      particles[index].stress[0][1] = 0;
      particles[index].stress[1][0] = 0;
      particles[index].stress[1][1] = 0;
      particles[index].strain[0][0] = 0;
      particles[index].strain[0][1] = 0;
      particles[index].strain[1][0] = 0;
      particles[index].strain[1][1] = 0;
      particles[index].F[0][0] = 1;
      particles[index].F[0][1] = 0;
      particles[index].F[1][0] = 0;
      particles[index].F[1][1] = 1;
      index++;
    }
  }
  particles[0].particle_count = index;

} // end of function initializeMaterial()
