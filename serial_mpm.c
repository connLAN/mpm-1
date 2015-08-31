/**
 *
 * serial_mpm.c
 *
 */

#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include<mpmstructs.h>
#include<initmpm.h>
#include<commfunctions.h>
#include<mpmfunctions.h>
#include<mpmio.h>

int main(int argc, char* argv[]){
  // domain variables
  int iter = 0;
  int file_index;
  int max_iter;
  double max_time = 1.0;
  double dt = .00005;
  double domain_box[4] = {0.0, 0.0, 1.0, 2.0};
  double domain_m_box[4] = {.2, 1.0, 0.8, 2.00};
  int domain_num_cells[2] = {10, 20};
  int pp_cell[2] = {2,2};
  double ipart_dim[2];
  int domain_num_particles[2];
  double test_px, test_py;
  char p_filename[40];
  char g_filename[40];

  //timing variables
  double start_time;
  double end_time;

  // patch variables
  int rank = 0;
  int size;
  int num_proc[2];
  double box[4];
  double m_box[4];
  int num_cells[2];
  int halo_cells[4] = {1,1,1,1};
  int num_particles[2];

  int *p1;
  int *p2;
  int *ph;
  

  //printf("pointer1: %d, pointer2: %d\n", test[0], &test[0][0]);

  GridData grid_data;
  Node** grid;

  InitialParticleData ipart_data;
  Particle* particles;

  max_iter = (int)(max_time/dt);

  if(argc == 4){
    domain_num_cells[0] = atoi(argv[1]);
    domain_num_cells[1] = atoi(argv[2]);
    max_iter = atoi(argv[3]);
  }

  num_proc[0] = 1;
  num_proc[1] = 1;
  rank = 0;

  grid_data.gravity[0] = 0;
  grid_data.gravity[1] = -900.8;
  grid_data.cell_dim[0] = (domain_box[2] - domain_box[0])/domain_num_cells[0];
  grid_data.cell_dim[1] = (domain_box[3] - domain_box[1])/domain_num_cells[1];
  grid_data.gridBoundaryType[0] = REFLEXIVE;
  grid_data.gridBoundaryType[1] = REFLEXIVE;
  grid_data.gridBoundaryType[2] = REFLEXIVE;
  grid_data.gridBoundaryType[3] = REFLEXIVE;

  test_px = grid_data.cell_dim[0]/pp_cell[0];
  test_py = grid_data.cell_dim[1]/pp_cell[1];

  domain_num_particles[0] = (int)((domain_m_box[2] - domain_m_box[0])/test_px);
  domain_num_particles[1] = (int)((domain_m_box[3] - domain_m_box[1])/test_py);

  ipart_dim[0] = (domain_m_box[2] - domain_m_box[0])/domain_num_particles[0];
  ipart_dim[1] = (domain_m_box[3] - domain_m_box[1])/domain_num_particles[1];


  decomposeGrid(rank, num_proc, domain_num_cells, num_cells, domain_box, box, halo_cells, &grid_data);
  decomposeMaterial(box, domain_m_box, m_box, ipart_dim, num_particles);

  ipart_data.density = 1000;
  ipart_data.bulk = 3;
  ipart_data.shear = 4;
  ipart_data.E = 40000;
  ipart_data.poisson = .30;
  ipart_data.idim[0] = ipart_dim[0];
  ipart_data.idim[1] = ipart_dim[1];
  ipart_data.velocity[0] = 0;
  ipart_data.velocity[1] = 0;
  ipart_data.box[0] = m_box[0];
  ipart_data.box[1] = m_box[1];
  ipart_data.box[2] = m_box[2];
  ipart_data.box[3] = m_box[3];
  ipart_data.num_particles[0] = num_particles[0];
  ipart_data.num_particles[1] = num_particles[1];
  ipart_data.domain_num_particles[0] = domain_num_particles[0];
  ipart_data.domain_num_particles[1] = domain_num_particles[1];

  grid = createGrid(&grid_data, box, num_cells);
  particles = createMaterial(&grid_data, &ipart_data);

  initializeGrid(&grid_data, grid);
  initializeMaterial(&grid_data, &ipart_data, particles);

  file_index = 0;

  start_time = MPI_Wtime();
  for(iter = 0; iter <= max_iter; iter++){
    clearGridValues(&grid_data, grid);
    mapToGrid(&grid_data, grid, particles);
    momentumToVelocityOnGrid(&grid_data, grid);
    computeForces(&grid_data, grid, particles);
    computeAcceleration(&grid_data, grid);
    updateNodalValues(&grid_data, grid, dt);
    updateStressAndStrain(&grid_data, grid, particles, dt);
    updateParticlePosition(&grid_data, grid, particles, dt);

    //print output
    if(iter%10 == 0){
      sprintf(p_filename, "spart%06d.vtk", file_index);
      //sprintf(g_filename, "grid_output%06d.vtk", file_index);
      writeParticlesVtkFile(&grid_data, grid, particles, p_filename);
      //writeGridVtkFile(&grid_data, grid, g_filename);
      file_index++;
    }
  }
  end_time = MPI_Wtime();

  printf("Serial, np: %d, iterations: %d, time: %f\n", 1, max_iter, end_time-start_time);
  //writeParticlesVtkFile(&grid_data, grid, particles, "sparts.vtk");
 
  freeGrid(grid, &grid_data);
  freeMaterial(particles);

  return 0;
}
