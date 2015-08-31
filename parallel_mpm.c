/**
 *
 * parallel_mpm.c
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
  int i;
  int file_index;
  int max_iter;
  double max_time = 1.0;
  double dt = .0000001;
  double domain_box[4] = {0.0, 0.0, 1.0, 2.0};
  double domain_m_box[4] = {.10, 0.2, 0.90, 2.00};
  int domain_num_cells[2] = {10, 20};
  int pp_cell[2] = {2,2};
  double ipart_dim[2];
  int domain_num_particles[2];
  double test_px, test_py;
  char p_filename[40];
  char p_filename2[40];
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
  int sfactor = 2;

  GridData grid_data;
  Node** grid;

  InitialParticleData ipart_data;
  Particle* particles;
  Particle* post_particles;
  Particle* particle_list;
  Particle* rhalo_parts[8];
  Particle* shalo_parts[8];

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  MPI_Request halo_req[8][2];
  MPI_Status halo_stat[8][2];

  max_iter = (int)(max_time/dt);

  if(rank == 0){
    if(argc == 3){
      num_proc[0] = atoi(argv[1]);
      num_proc[1] = atoi(argv[2]);
    }else if(argc == 6){
      num_proc[0] = atoi(argv[1]);
      num_proc[1] = atoi(argv[2]);
      domain_num_cells[0] = atoi(argv[3]);
      domain_num_cells[1] = atoi(argv[4]);
      max_iter = atoi(argv[5]);
    }
  }

  MPI_Bcast(num_proc,         2, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(domain_num_cells, 2, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&max_iter,        1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
 
  grid_data.gravity[0] = 0;
  grid_data.gravity[1] = -900.8;
  grid_data.cell_dim[0] = (domain_box[2] - domain_box[0])/domain_num_cells[0];
  grid_data.cell_dim[1] = (domain_box[3] - domain_box[1])/domain_num_cells[1];

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
  ipart_data.E = 90000;
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
  post_particles = createMaterial(&grid_data, &ipart_data);

  initializeGrid(&grid_data, grid);
	initializeMaterial(&grid_data, &ipart_data, particles);
  
  //allocate halo particles
  allocateHaloParticles(&grid_data, sfactor, pp_cell, rhalo_parts, shalo_parts);
  
  file_index = 0;

  start_time = MPI_Wtime();
  for(iter = 0; iter <= max_iter; iter++){
    gatherHaloParticles(&grid_data, particles, rhalo_parts, shalo_parts);
    sendRecvParticles(&grid_data, rhalo_parts, shalo_parts);
    clearGridValues(&grid_data, grid);
    mapToGrid(&grid_data, grid, particles);
    for(i = 0; i < 8; i++){
      if(grid_data.rank[i] >= 0 && rhalo_parts[i][0].particle_count > 0){
        mapToGrid(&grid_data, grid, rhalo_parts[i]);
      }
    }
    momentumToVelocityOnGrid(&grid_data, grid);
    computeForces(&grid_data, grid, particles);
    for(i = 0; i < 8; i++){
      if(grid_data.rank[i] >= 0 && rhalo_parts[i][0].particle_count > 0){
        computeForces(&grid_data, grid, rhalo_parts[i]);
      }
    }
    computeAcceleration(&grid_data, grid);

    /**
    particle_list = gatherParticles(rank, size, domain_num_particles[0] * domain_num_particles[1],
                                    particles, post_particles);
    if(rank == 0){
      //if(iter%100 == 0){
      sprintf(p_filename, "ppart%06d.vtk", file_index);
      //sprintf(p_filename2, "center%06d.vtk", file_index);
      //sprintf(g_filename, "grid_output%06d.vtk", file_index);
      //writeParticlesVtkFile(&grid_data, grid, rhalo_parts[0], p_filename);
      writeParticlesVtkFile(&grid_data, grid, particle_list, p_filename);
      //writeParticlesVtkFile(&grid_data, grid, particle_list, p_filename);
      //writeParticlesVtkFile(&grid_data, grid, particles, p_filename);
      //writeGridVtkFile(&grid_data, grid, g_filename);
      //writeGridVtkFile(&grid_data, grid, g_filename);
      //free(particle_list);
      file_index++;
      //}
    }
    **/

    updateNodalValues(&grid_data, grid, dt);
    updateStressAndStrain(&grid_data, grid, particles, dt);
    pUpdateParticlePosition(&grid_data, grid, particles, post_particles, shalo_parts, dt);
    sendRecvParticles(&grid_data, rhalo_parts, shalo_parts);
    updateParticleList(&grid_data, particles, rhalo_parts);

    /**
    particle_list = gatherParticles(rank, size, domain_num_particles[0] * domain_num_particles[1],
                                    particles, post_particles);
    if(rank == 0){
      sprintf(p_filename, "particle_output%06d.vtk", file_index);
      file_index++;
      //sprintf(g_filename, "grid_output%06d.vtk", file_index);
      writeParticlesVtkFile(&grid_data, grid, particles, p_filename);
      //writeParticlesVtkFile(&grid_data, grid, particle_list, p_filename);
      //writeParticlesVtkFile(&grid_data, grid, particles, p_filename);
      //writeGridVtkFile(&grid_data, grid, g_filename);
      free(particle_list);
    }
    **/

    MPI_Barrier(MPI_COMM_WORLD);
  }
  end_time = MPI_Wtime();

  MPI_Barrier(MPI_COMM_WORLD);
  //particle_list = gatherParticles(rank, size, domain_num_particles[0] * domain_num_particles[1],
  //                                particles, post_particles);
  if(rank == 0){
    printf("Parallel, np: %d, iterations: %d, time: %f\n", size, max_iter, end_time-start_time);
    //writeParticlesVtkFile(&grid_data, grid, particle_list, "pparts.vtk");
    //free(particle_list);
  }

  
  freeGrid(grid, &grid_data);
  freeMaterial(particles);
  freeMaterial(post_particles);
  freeHaloParticles(&grid_data, rhalo_parts, shalo_parts);

  //MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

  return 0;
}
