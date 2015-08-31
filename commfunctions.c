/**
 *
 * commfunctions.c
 *
 */

#include<commfunctions.h>
#include<stdio.h>
#include<stdlib.h>

void decomposeGrid(int rank, int *num_proc, int *domain_nc, int *num_cells, double *domain_box, double *box,
                   int *halo_cells, GridData *grid_data){

  int i;
  int remI = domain_nc[0]%num_proc[0];
  int remJ = domain_nc[1]%num_proc[1];
  int patch_index[2];

  getIndex(rank, num_proc[0], num_proc[1], &patch_index[0], &patch_index[1]); 

  for(i = 0; i < 8; i++)
    grid_data->rank[i] = -10;

  if(patch_index[0] < remI){
    num_cells[0] = (domain_nc[0]/num_proc[0]) + 1;
    box[0] = domain_box[0] + grid_data->cell_dim[0] * ((double)(patch_index[0] * num_cells[0]));
  }else{
    num_cells[0] = (domain_nc[0]/num_proc[0]);
    box[0] = domain_box[0] + grid_data->cell_dim[0] * ((double)(remI + patch_index[0] * num_cells[0]));
  }
  box[2] = box[0] + (grid_data->cell_dim[0] * ((double)num_cells[0]));

  if(patch_index[1] < remJ){
    num_cells[1] = (domain_nc[1]/num_proc[1]) + 1;
    box[1] = domain_box[1] + grid_data->cell_dim[1] * ((double)(patch_index[1] * num_cells[1]));
  }else{
    num_cells[1] = (domain_nc[1]/num_proc[1]);
    box[1] = domain_box[1] + grid_data->cell_dim[1] * ((double)(remJ + patch_index[1] * num_cells[1]));
  }
  box[3] = box[1] + (grid_data->cell_dim[1] * ((double)num_cells[1]));

  if(patch_index[0] == 0){
    grid_data->gridBoundaries[0] = BOUNDARY;
    grid_data->halo_cells[0] = 0;
  }else{
    grid_data->gridBoundaries[0] = HALO;
    grid_data->halo_cells[0] = halo_cells[0];
    grid_data->rank[0] = getRank(num_proc[0], num_proc[1], patch_index[0]-1, patch_index[1]);
  }

  if(patch_index[0] == (num_proc[0]-1)){
    grid_data->gridBoundaries[2] = BOUNDARY;
    grid_data->halo_cells[2] = 0;
  }else{
    grid_data->gridBoundaries[2] = HALO;
    grid_data->halo_cells[2] = halo_cells[2];
    grid_data->rank[4] = getRank(num_proc[0], num_proc[1], patch_index[0]+1, patch_index[1]);
  }

  if(patch_index[1] == 0){
    grid_data->gridBoundaries[1] = BOUNDARY;
    grid_data->halo_cells[1] = 0;
  }else{
    grid_data->gridBoundaries[1] = HALO;
    grid_data->halo_cells[1] = halo_cells[1];
    grid_data->rank[2] = getRank(num_proc[0], num_proc[1], patch_index[0], patch_index[1]-1);
  }

  if(patch_index[1] == (num_proc[1]-1)){
    grid_data->gridBoundaries[3] = BOUNDARY;
    grid_data->halo_cells[3] = 0;
  }else{
    grid_data->gridBoundaries[3] = HALO;
    grid_data->halo_cells[3] = halo_cells[3];
    grid_data->rank[6] = getRank(num_proc[0], num_proc[1], patch_index[0], patch_index[1]+1);
  }

  if(grid_data->gridBoundaries[0] == HALO && grid_data->gridBoundaries[1] == HALO){
    grid_data->rank[1] = getRank(num_proc[0], num_proc[1], patch_index[0]-1, patch_index[1]-1);
  }

  if(grid_data->gridBoundaries[1] == HALO && grid_data->gridBoundaries[2] == HALO){
    grid_data->rank[3] = getRank(num_proc[0], num_proc[1], patch_index[0]+1, patch_index[1]-1);
  }

  if(grid_data->gridBoundaries[2] == HALO && grid_data->gridBoundaries[3] == HALO){
    grid_data->rank[5] = getRank(num_proc[0], num_proc[1], patch_index[0]+1, patch_index[1]+1);
  }

  if(grid_data->gridBoundaries[3] == HALO && grid_data->gridBoundaries[0] == HALO){
    grid_data->rank[7] = getRank(num_proc[0], num_proc[1], patch_index[0]-1, patch_index[1]+1);
  }

  grid_data->grid_rank = rank;
  grid_data->patch_index.i = patch_index[0];
  grid_data->patch_index.j = patch_index[1];
  grid_data->num_proc[0] = num_proc[0];
  grid_data->num_proc[1] = num_proc[1];
  grid_data->parallel = MPM_TRUE;
}

void decomposeMaterial(double *box, double *m_dom_box, double *m_box, double *ipart_dim, int *num_part){
  int nump_x0, nump_y0;
  int nump_x1, nump_y1;
  MpmBool compute = MPM_FALSE;
  double phx = ipart_dim[0];
  double phy = ipart_dim[1];
  double ph2x = ipart_dim[0]/2;
  double ph2y = ipart_dim[1]/2;

  double dmX0 = m_dom_box[0];
  double dmY0 = m_dom_box[1];
  double dmX1 = m_dom_box[2];
  double dmY1 = m_dom_box[3];
  double pmX0 = dmX0 + ph2x;
  double pmY0 = dmY0 + ph2y;
  double pmX1 = dmX1 - ph2x;
  double pmY1 = dmY1 - ph2y;
  double bX0 = box[0];
  double bY0 = box[1];
  double bX1 = box[2];
  double bY1 = box[3];

  //printf("phx: %f, phy: %f\n", phx, phy);
  //printf("pmX0: %f, pmX1: %f\n", pmX0, pmX1);
  //printf("pmY0: %f, pmY1: %f\n", pmY0, pmY1);

  if(pmX0 <= bX0 && pmX1 >= bX1){
    nump_x0 = 1 + (int)((bX0 - pmX0)/phx);
    nump_x1 = 1 + (int)((bX1 - pmX0)/phx);
    if(pmY0 <= bY0 && pmY1 >= bY1){
      nump_y0 = 1 + (int)((bY0 - pmY0)/phy);
      nump_y1 = 1 + (int)((bY1 - pmY0)/phy);
      compute = MPM_TRUE;
    }else if(pmY0 > bY0 && pmY0 < bY1 && pmY1 >= bY1){
      nump_y0 = 0;
      nump_y1 = 1 + (int)((bY1 - pmY0)/phy);
      compute = MPM_TRUE;
    }else if(pmY0 <= bY0 && pmY1 > bY0 && pmY1 < bY1){
      nump_y0 = 1 + (int)((bY0 - pmY0)/phy);
      nump_y1 = 1 + (int)((pmY1 - dmY0)/phy);
      compute = MPM_TRUE;
    }else if(pmY0 > bY0 && pmY0 < bY1 && pmY1 > bY0 && pmY1 < bY1){
      nump_y0 = 0;
      nump_y1 = 1 + (int)((pmY1 - dmY0)/phy);
      compute = MPM_TRUE;
    }
  }else if(pmX0 > bX0 && pmX0 < bX1 && pmX1 >= bX1){
    nump_x0 = 0;
    nump_x1 = 1 + (int)((bX1 - pmX0)/phx);
    if(pmY0 <= bY0 && pmY1 >= bY1){
      nump_y0 = 1 + (int)((bY0 - pmY0)/phy);
      nump_y1 = 1 + (int)((bY1 - pmY0)/phy);
      compute = MPM_TRUE;
    }else if(pmY0 > bY0 && pmY0 < bY1 && pmY1 >= bY1){
      nump_y0 = 0;
      nump_y1 = 1 + (int)((bY1 - pmY0)/phy);
      compute = MPM_TRUE;
    }else if(pmY0 <= bY0 && pmY1 > bY0 && pmY1 < bY1){
      nump_y0 = 1 + (int)((bY0 - pmY0)/phy);
      nump_y1 = 1 + (int)((pmY1 - dmY0)/phy);
      compute = MPM_TRUE;
    }else if(pmY0 > bY0 && pmY0 < bY1 && pmY1 > bY0 && pmY1 < bY1){
      nump_y0 = 0;
      nump_y1 = 1 + (int)((pmY1 - dmY0)/phy);
      compute = MPM_TRUE;
    }
  }else if(pmX0 <= bX0 && pmX1 > bX0 && pmX1 < bX1){
    nump_x0 = 1 + (int)((bX0 - pmX0)/phx);
    nump_x1 = 1 + (int)((pmX1 - dmX0)/phx);
    if(pmY0 <= bY0 && pmY1 >= bY1){
      nump_y0 = 1 + (int)((bY0 - pmY0)/phy);
      nump_y1 = 1 + (int)((bY1 - pmY0)/phy);
      compute = MPM_TRUE;
    }else if(pmY0 > bY0 && pmY0 < bY1 && pmY1 >= bY1){
      nump_y0 = 0;
      nump_y1 = 1 + (int)((bY1 - pmY0)/phy);
      compute = MPM_TRUE;
    }else if(pmY0 <= bY0 && pmY1 > bY0 && pmY1 < bY1){
      nump_y0 = 1 + (int)((bY0 - pmY0)/phy);
      nump_y1 = 1 + (int)((pmY1 - dmY0)/phy);
      compute = MPM_TRUE;
    }else if(pmY0 > bY0 && pmY0 < bY1 && pmY1 > bY0 && pmY1 < bY1){
      nump_y0 = 0;
      nump_y1 = 1 + (int)((pmY1 - dmY0)/phy);
      compute = MPM_TRUE;
    }
  }else if(pmX0 > bX0  && pmX0 < bX1 && pmX1 > bX0 && pmX1 < bX1){
    nump_x0 = 0;
    nump_x1 = 1 + (int)((pmX1 - dmX0)/phx);
    if(pmY0 <= bY0 && pmY1 >= bY1){
      nump_y0 = 1 + (int)((bY0 - pmY0)/phy);
      nump_y1 = 1 + (int)((bY1 - pmY0)/phy);
      compute = MPM_TRUE;
    }else if(pmY0 > bY0 && pmY0 < bY1 && pmY1 >= bY1){
      nump_y0 = 0;
      nump_y1 = 1 + (int)((bY1 - pmY0)/phy);
      compute = MPM_TRUE;
    }else if(pmY0 <= bY0 && pmY1 > bY0 && pmY1 < bY1){
      nump_y0 = 1 + (int)((bY0 - pmY0)/phy);
      nump_y1 = 1 + (int)((pmY1 - dmY0)/phy);
      compute = MPM_TRUE;
    }else if(pmY0 > bY0 && pmY0 < bY1 && pmY1 > bY0 && pmY1 < bY1){
      nump_y0 = 0;
      nump_y1 = 1 + (int)((pmY1 - dmY0)/phy);
      compute = MPM_TRUE;
    }
  }
  if(compute == MPM_TRUE){
      m_box[0] = dmX0 + nump_x0*phx;
      m_box[1] = dmY0 + nump_y0*phy;
      m_box[2] = dmX0 + nump_x1*phx;
      m_box[3] = dmY0 + nump_y1*phy;
      num_part[0] = nump_x1 - nump_x0;
      num_part[1] = nump_y1 - nump_y0;
  }else{
    m_box[0] = 0;
    m_box[1] = 0;
    m_box[2] = 0;
    m_box[3] = 0;
    num_part[0] = 0;
    num_part[1] = 0;
  }
  //printf("m_box[0]: %f, m_box[2]: %f, num_part[0]: %d\n", m_box[0], m_box[2], num_part[0]);
  //printf("m_box[1]: %f, m_box[3]: %f, num_part[1]: %d\n", m_box[1], m_box[3], num_part[1]);
} // end of decomposeGrid()

void allocateHaloParticles(GridData *grid_data, int sfactor, int *pp_cell,
                           Particle **rhalo_parts, Particle **shalo_parts)
{
  int index, i;
  int init_c_size[4];
  int init_corner[4];
  int init_ppcell = sfactor * pp_cell[0] * pp_cell[1];
  init_c_size[0] = init_ppcell * grid_data->num_nodes[1] * grid_data->halo_cells[0];
  init_c_size[1] = init_ppcell * grid_data->num_nodes[0] * grid_data->halo_cells[1];
  init_c_size[2] = init_ppcell * grid_data->num_nodes[1] * grid_data->halo_cells[2];
  init_c_size[3] = init_ppcell * grid_data->num_nodes[0] * grid_data->halo_cells[3];
  init_corner[0] = init_ppcell * grid_data->halo_cells[0] * grid_data->halo_cells[1];
  init_corner[1] = init_ppcell * grid_data->halo_cells[1] * grid_data->halo_cells[2];
  init_corner[2] = init_ppcell * grid_data->halo_cells[2] * grid_data->halo_cells[3];
  init_corner[3] = init_ppcell * grid_data->halo_cells[3] * grid_data->halo_cells[0];
  
  index = 0;
  for(i = 0; i <= 6; i+=2){
    if(grid_data->gridBoundaries[index] == HALO){
      grid_data->phalo_size[i] = init_c_size[index];
      rhalo_parts[i] = (Particle*) malloc(grid_data->phalo_size[i] * sizeof(Particle));
      shalo_parts[i] = (Particle*) malloc(grid_data->phalo_size[i] * sizeof(Particle));
      rhalo_parts[i][0].particle_count = 0;
      shalo_parts[i][0].particle_count = 0;
    }else{
      grid_data->phalo_size[i] = 0;
    }
    index++;
  }

  index = 0;
  for(i = 1; i <= 7; i+=2){
    if(i != 7){
      if(grid_data->gridBoundaries[index] == HALO && grid_data->gridBoundaries[index+1] == HALO){
        grid_data->phalo_size[i] = init_corner[index];
        rhalo_parts[i] = (Particle*) malloc(grid_data->phalo_size[i] * sizeof(Particle));
        shalo_parts[i] = (Particle*) malloc(grid_data->phalo_size[i] * sizeof(Particle));
        rhalo_parts[i][0].particle_count = 0;
        shalo_parts[i][0].particle_count = 0;
      }else{
        grid_data->phalo_size[i] = 0;
      }
    }else{
      if(grid_data->gridBoundaries[index] == HALO && grid_data->gridBoundaries[0] == HALO){
        grid_data->phalo_size[i] = init_corner[index];
        rhalo_parts[i] = (Particle*) malloc(grid_data->phalo_size[i] * sizeof(Particle));
        shalo_parts[i] = (Particle*) malloc(grid_data->phalo_size[i] * sizeof(Particle));
        rhalo_parts[i][0].particle_count = 0;
        shalo_parts[i][0].particle_count = 0;
      }else{
        grid_data->phalo_size[i] = 0;
      }
    }
    index++;
  }
} //end of allocateHaloParticles

void gatherHaloParticles(GridData *grid_data, Particle *particles, Particle **rhalo_parts, Particle **shalo_parts){
  int i,j,index,hindex[8];
  int particle_count = particles[0].particle_count;
  int haloX0 = grid_data->halo_cells[0];
  int haloY0 = grid_data->halo_cells[1];
  int haloi00 = -10;
  int haloi01 = -10;
  int haloj0 = -10;
  int haloi1 = -10;
  int haloj1 = -10;
 
  if(grid_data->gridBoundaries[0] == HALO){
    haloi00 = grid_data->halo_cells[0];
    haloi01 = haloi00 + grid_data->halo_cells[0] - 1;
  }
  if(grid_data->gridBoundaries[1] == HALO)
    haloj0 = grid_data->halo_cells[1];
  if(grid_data->gridBoundaries[2] == HALO)
    haloi1 = grid_data->num_nodes[0] - grid_data->halo_cells[2] - 2;
  if(grid_data->gridBoundaries[3] == HALO)
    haloj1 = grid_data->num_nodes[1] - grid_data->halo_cells[3] - 2;
  
  double bx0, by0;
  double nx0, ny0;
  double hx = grid_data->cell_dim[0];
  double hy = grid_data->cell_dim[1];
  
  bx0 = grid_data->box[0] - (haloX0 * hx);
  by0 = grid_data->box[1] - (haloY0 * hy);
  
  for(i = 0; i < 8; i++){
    hindex[i] = 0;
  }
  for(index = 0; index < particle_count; index++){
    nx0 = particles[index].position[0];
    ny0 = particles[index].position[1];
    i = (int)((nx0 - bx0)/hx);
    j = (int)((ny0 - by0)/hy);
    
    if(i >= haloi00 && i <= haloi01){
      shalo_parts[0][hindex[0]] = particles[index];
      hindex[0]++;
      if(j == haloj0){
        shalo_parts[1][hindex[1]] = particles[index];
        hindex[1]++;
      }
    }
    if(j == haloj0){
      shalo_parts[2][hindex[2]] = particles[index];
      hindex[2]++;
      if(i == haloi1){
        shalo_parts[3][hindex[3]] = particles[index];
        hindex[3]++;
      }
    }
    if(i == haloi1){
      shalo_parts[4][hindex[4]] = particles[index];
      hindex[4]++;
      if(j == haloj1){
        shalo_parts[5][hindex[5]] = particles[index];
        hindex[5]++;
      }
    }
    if(j == haloj1){
      shalo_parts[6][hindex[6]] = particles[index];
      hindex[6]++;
      if(i >= haloi00 && i <= haloi01){
        shalo_parts[7][hindex[7]] = particles[index];
        hindex[7]++;
      }
    }
  }
  
  for(i = 0; i < 8; i++){
    if(grid_data->phalo_size[i] > 0){
      shalo_parts[i][0].particle_count = hindex[i];
      rhalo_parts[i][0].particle_count = 0;
      grid_data->halo_particle_count[i] = hindex[i];
    }
  }
  
} //end of gatherHaloParticles

void wpGatherHaloParticles(GridData **grid_data, Particle ***particles, Particle ***new_particles,
                           Particle ***halo_particles)
{
  int i,j,index;
  int pi,pj;
	int npx = grid_data[0][0].num_proc[0];
	int npy = grid_data[0][0].num_proc[1];
	int mem_size;
	int halo_index;
  int particle_count;
  int haloX0;
  int haloY0;
  int haloi0;
  int haloj0;
  int haloi1;
  int haloj1;
  double bx0, by0;
  double nx0, ny0;
  double hx;
  double hy;
 
	for(j = 0; j < npy; j++){
	  for(i = 0; i < npx; i++){
      halo_particles[i][j][0].particle_count = 0;
      new_particles[i][j][0].particle_count = 0;
		}
	}

	for(j = 0; j < npy; j++){
	  for(i = 0; i < npx; i++){
      particle_count = particles[i][j][0].particle_count;
      haloi0 = -10;
      haloj0 = -10;
      haloi1 = -10;
      haloj1 = -10;

      if(grid_data[i][j].gridBoundaries[0] == HALO)
        haloi0 = grid_data[i][j].halo_cells[0];
      if(grid_data[i][j].gridBoundaries[1] == HALO)
        haloj0 = grid_data[i][j].halo_cells[1];
      if(grid_data[i][j].gridBoundaries[2] == HALO)
        haloi1 = grid_data[i][j].num_nodes[0] - grid_data[i][j].halo_cells[2] - 2;
      if(grid_data[i][j].gridBoundaries[3] == HALO)
        haloj1 = grid_data[i][j].num_nodes[1] - grid_data[i][j].halo_cells[3] - 2;
      
      haloX0 = grid_data[i][j].halo_cells[0];
      haloY0 = grid_data[i][j].halo_cells[1];

      hx = grid_data[i][j].cell_dim[0];
      hy = grid_data[i][j].cell_dim[1];
      
      bx0 = grid_data[i][j].box[0] - (haloX0 * hx);
      by0 = grid_data[i][j].box[1] - (haloY0 * hy);

      for(index = 0; index < particle_count; index++){
        nx0 = particles[i][j][index].position[0];
        ny0 = particles[i][j][index].position[1];
        pi = (int)((nx0 - bx0)/hx);
        pj = (int)((ny0 - by0)/hy);
       
        if(pi == haloi0){
          halo_index = halo_particles[i-1][j][0].particle_count;
          halo_particles[i-1][j][halo_index] = particles[i][j][index];
          halo_particles[i-1][j][0].particle_count = halo_index + 1;
          if(pj == haloj0){
            halo_index = halo_particles[i-1][j-1][0].particle_count;
            halo_particles[i-1][j-1][halo_index] = particles[i][j][index];
            halo_particles[i-1][j-1][0].particle_count = halo_index + 1;
          }
        }
        if(pj == haloj0){
          halo_index = halo_particles[i][j-1][0].particle_count;
          halo_particles[i][j-1][halo_index] = particles[i][j][index];
          halo_particles[i][j-1][0].particle_count = halo_index + 1;
          if(pi == haloi1){
            halo_index = halo_particles[i+1][j-1][0].particle_count;
            halo_particles[i+1][j-1][halo_index] = particles[i][j][index];
            halo_particles[i+1][j-1][0].particle_count = halo_index + 1;
          }
        }
        if(pi == haloi1){
          halo_index = halo_particles[i+1][j][0].particle_count;
          halo_particles[i+1][j][halo_index] = particles[i][j][index];
          halo_particles[i+1][j][0].particle_count = halo_index + 1;
          if(pj == haloj1){
            halo_index = halo_particles[i+1][j+1][0].particle_count;
            halo_particles[i+1][j+1][halo_index] = particles[i][j][index];
            halo_particles[i+1][j+1][0].particle_count = halo_index + 1;
          }
        }
        if(pj == haloj1){
          halo_index = halo_particles[i][j+1][0].particle_count;
          halo_particles[i][j+1][halo_index] = particles[i][j][index];
          halo_particles[i][j+1][0].particle_count = halo_index + 1;
          if(pi == haloi0){
            halo_index = halo_particles[i-1][j+1][0].particle_count;
            halo_particles[i-1][j+1][halo_index] = particles[i][j][index];
            halo_particles[i-1][j+1][0].particle_count = halo_index + 1;
          }
        }
      }
		} //end j loop
	} //end i loop
} //end of wpGatherHaloParticles

void sendRecvParticles(GridData *grid_data, Particle **rhalo_parts, Particle **shalo_parts){
 
  int i, rank, rsize, ssize;
  int particle_size = sizeof(Particle);
  MPI_Request halo_req[8][2];
  MPI_Status halo_stat[8][2];
  
  for(i = 0; i < 8; i++){
    rank = grid_data->rank[i];
    if(rank >= 0){
      rsize = grid_data->phalo_size[i];
      ssize = shalo_parts[i][0].particle_count;
      rhalo_parts[i][0].particle_count = 0;
      MPI_Irecv(rhalo_parts[i], rsize*particle_size, MPI_BYTE, rank, 0, MPI_COMM_WORLD, &halo_req[i][0]);
      MPI_Isend(shalo_parts[i], ssize*particle_size, MPI_BYTE, rank, 0, MPI_COMM_WORLD, &halo_req[i][1]);
    }
  } 
  for(i = 0; i < 8; i++){
    rank = grid_data->rank[i];
    if(rank >= 0){
      MPI_Waitall(2, halo_req[i], halo_stat[i]);
    }
  } 
} //end of sendRecvParticles

void sendRecvParticles1(GridData *grid_data, Particle **rhalo_parts, Particle **shalo_parts,
                        MPI_Request **halo_req, MPI_Status **halo_stat)
{
  int i, rank, rsize, ssize;
  int particle_size = sizeof(Particle);
  
  for(i = 0; i < 8; i++){
    rank = grid_data->rank[i];
    if(rank >= 0){
      rsize = grid_data->phalo_size[i];
      ssize = shalo_parts[i][0].particle_count;
      rhalo_parts[i][0].particle_count = 0;
      MPI_Irecv(rhalo_parts[i], rsize*particle_size, MPI_BYTE, rank, 0, MPI_COMM_WORLD, &halo_req[i][0]);
      MPI_Isend(shalo_parts[i], ssize*particle_size, MPI_BYTE, rank, 0, MPI_COMM_WORLD, &halo_req[i][1]);
    }
  } 
} //end of sendRecvParticles1

void sendRecvParticles2(GridData *grid_data, MPI_Request **halo_req, MPI_Status **halo_stat)
{
  int i, rank;
  
  for(i = 0; i < 8; i++){
    rank = grid_data->rank[i];
    if(rank >= 0){
      MPI_Waitall(2, halo_req[i], halo_stat[i]);
    }
  } 
} //end of sendRecvParticles2

void freeHaloParticles(GridData *grid_data, Particle **rhalo_parts, Particle **shalo_parts){
  int i;
  for(i = 0; i < 8; i++){
    if(grid_data->phalo_size[i] > 0){
      free(rhalo_parts[i]);
      free(shalo_parts[i]);
    }
  }
}

void updateParticleList(GridData *grid_data, Particle *particles, Particle **rhalo_parts){
  int i,j;
  int rank;
  int index = particles[0].particle_count;
  int r_particle_count;

  for(i = 0; i < 8; i++){
    rank = grid_data->rank[i];
    if(rank >= 0){
      r_particle_count = rhalo_parts[i][0].particle_count;
      for(j = 0; j < r_particle_count; j++){
        particles[index] = rhalo_parts[i][j];
        index++;
      }
    }
  }
  particles[0].particle_count = index;


} //end of updateParticleList

Particle* gatherParticles(int rank, int size, int array_size, Particle *particles, Particle *post_particles){
  int i,j,index;
  int particle_size = sizeof(Particle);
  int particle_count = particles[0].particle_count;
  int r_particle_count;
  Particle *particle_list = NULL;

  if(rank == 0){
    particle_list = (Particle*) malloc(array_size * particle_size);
    for(index = 0; index < particle_count; index++){
      particle_list[index] = particles[index];
    }

    for(i = 1; i < size; i++){
      MPI_Recv(post_particles, array_size * particle_size, MPI_BYTE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      r_particle_count = post_particles[0].particle_count;
      for(j = 0; j < r_particle_count; j++){
        particle_list[index] = post_particles[j];
        index++;
      }
    }
    particle_list[0].particle_count = index;
  }else{
      MPI_Send(particles, particle_count * particle_size, MPI_BYTE, 0, 0, MPI_COMM_WORLD);
  }
  
  return particle_list;
}

void wpRelocateParticles(GridData **grids, Particle *particles, Particle ***particles_new)
{
  int i, j, index1, index2;
  int nx, ny;
  int particle_count = particles[0].particle_count;
  double px,py;
  GridData grid_data;
  MpmBool particle_removed;

  i = particles[0].gridID.i;
  j = particles[0].gridID.j;
  grid_data = grids[i][j];

  for(index1 = 0; index1 < particle_count; index1++){
    particle_removed = MPM_FALSE;
    px = particles[index1].position[0];
    py = particles[index1].position[1];

    if(px < grid_data.box[0] && grid_data.gridBoundaries[0] == HALO){
      if(py < grid_data.box[1] && grid_data.gridBoundaries[1] == HALO){
        index2 = particles_new[i-1][j-1][0].particle_count;
        particles_new[i-1][j-1][index2] = particles[index1];
        particles_new[i-1][j-1][0].particle_count = index2 + 1;
        particle_removed = MPM_TRUE;
      }else if(py > grid_data.box[3] && grid_data.gridBoundaries[3] == HALO){
        index2 = particles_new[i-1][j+1][0].particle_count;
        particles_new[i-1][j+1][index2] = particles[index1];
        particles_new[i-1][j+1][0].particle_count = index2 + 1;
        particle_removed = MPM_TRUE;
      }else{
        index2 = particles_new[i-1][j][0].particle_count;
        particles_new[i-1][j][index2] = particles[index1];
        particles_new[i-1][j][0].particle_count = index2 + 1;
        particle_removed = MPM_TRUE;
      }
    }else if(px > grid_data.box[2] && grid_data.gridBoundaries[2] == HALO){
      if(py < grid_data.box[1] && grid_data.gridBoundaries[1] == HALO){
        index2 = particles_new[i+1][j-1][0].particle_count;
        particles_new[i+1][j-1][index2] = particles[index1];
        particles_new[i+1][j-1][0].particle_count = index2 + 1;
        particle_removed = MPM_TRUE;
      }else if(py > grid_data.box[3] && grid_data.gridBoundaries[3] == HALO){
        index2 = particles_new[i+1][j+1][0].particle_count;
        particles_new[i+1][j+1][index2] = particles[index1];
        particles_new[i+1][j+1][0].particle_count = index2 + 1;
        particle_removed = MPM_TRUE;
      }else{
        index2 = particles_new[i+1][j][0].particle_count;
        particles_new[i+1][j][index2] = particles[index1];
        particles_new[i+1][j][0].particle_count = index2 + 1;
        particle_removed = MPM_TRUE;
      }
    }

    if(py < grid_data.box[1] && grid_data.gridBoundaries[1] == HALO){
      if(px > grid_data.box[0] && px < grid_data.box[2]){
        index2 = particles_new[i][j-1][0].particle_count;
        particles_new[i][j-1][index2] = particles[index1];
        particles_new[i][j-1][0].particle_count = index2 + 1;
        particle_removed = MPM_TRUE;
      }
    }else if(py > grid_data.box[3] && grid_data.gridBoundaries[3] == HALO){
      if(px > grid_data.box[0] && px < grid_data.box[2]){
        index2 = particles_new[i][j+1][0].particle_count;
        particles_new[i][j+1][index2] = particles[index1];
        particles_new[i][j+1][0].particle_count = index2 + 1;
        particle_removed = MPM_TRUE;
      }
    }
    
    if(particle_removed == MPM_FALSE){
      index2 = particles_new[i][j][0].particle_count;
      particles_new[i][j][index2] = particles[index1];
      particles_new[i][j][0].particle_count = index2 + 1;
    }else{
      printf("particle moved: %f, %f\n", px, py);
      printf("b0: %f, b1: %f, b2: %f, b3: %f\n", grid_data.box[0], grid_data.box[1], grid_data.box[2], grid_data.box[3]);
    }
  }
}

void getIndex(int rank, int procX, int procY, int* gridI, int* gridJ){
  int i,j;

  for(j = 0; j < procY; j++){
    for(i = 0; i < procX; i ++){
      if(rank == i + j*procX){
        *gridI = i;
        *gridJ = j;
        return;
      }
    }
  }

}

int getRank(int procX, int procY, int gridI, int gridJ){
  return gridI + gridJ * procX;
}
