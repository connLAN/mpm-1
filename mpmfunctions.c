/**
 *
 * mpmfunctions.c
 *
 **/

#include<mpmfunctions.h>
#include<stdio.h>

void clearGridValues(GridData *grid_data, Node** grid){
  int i,j;
  int num_nodesX = grid_data->num_nodes[0];
  int num_nodesY = grid_data->num_nodes[1];
  for(i = 0; i < num_nodesX; i++){
    for(j = 0; j < num_nodesY; j++){
      grid[i][j].mass = 0;
      grid[i][j].int_force[0] = 0;
      grid[i][j].int_force[1] = 0;
      grid[i][j].body_force[0] = 0;
      grid[i][j].body_force[1] = 0;
      grid[i][j].momentum[0] = 0;
      grid[i][j].momentum[1] = 0;
      grid[i][j].velocity[0] = 0;
      grid[i][j].velocity[1] = 0;
      grid[i][j].acceleration[0] = 0;
      grid[i][j].acceleration[1] = 0;
      grid[i][j].num_particles = 0;
      grid[i][j].has_mass = MPM_FALSE;
    }
  }
}

void wpClearGridValues(GridData *grid_data, Node** grid){
  int i,j;
  int num_nodesX = grid_data->num_nodes[0];
  int num_nodesY = grid_data->num_nodes[1];
  int haloX0 = grid_data->halo_cells[0];
  int haloY0 = grid_data->halo_cells[1];
  int haloX1 = grid_data->halo_cells[2];
  int haloY1 = grid_data->halo_cells[3];

  double hx = grid_data->cell_dim[0];
  double hy = grid_data->cell_dim[1];
  double bx0 = grid_data->box[0] - (haloX0 * hx);
  double by0 = grid_data->box[1] - (haloY0 * hy);


  for(i = 0; i < num_nodesX; i++){
    for(j = 0; j < num_nodesY; j++){
      grid[i][j].mass = 0;
      grid[i][j].int_force[0] = 0;
      grid[i][j].int_force[1] = 0;
      grid[i][j].body_force[0] = 0;
      grid[i][j].body_force[1] = 0;
      grid[i][j].momentum[0] = 0;
      grid[i][j].momentum[1] = 0;
      grid[i][j].velocity[0] = 0;
      grid[i][j].velocity[1] = 0;
      grid[i][j].acceleration[0] = 0;
      grid[i][j].acceleration[1] = 0;
      grid[i][j].num_particles = 0;
      grid[i][j].has_mass = MPM_FALSE;
      grid[i][j].position[0] = bx0 + i * hx;
      grid[i][j].position[1] = by0 + j * hy;
    }
  }
}

void mapToGrid(GridData *grid_data, Node** grid, Particle* particles)
{

  int i,j,index;
  int particle_count = particles[0].particle_count;
  int haloX0 = grid_data->halo_cells[0];
  int haloY0 = grid_data->halo_cells[1];
  
  double bx0, by0;
  double nx0, ny0;
  double wx, wy, dwx, dwy;
  double hx = grid_data->cell_dim[0];
  double hy = grid_data->cell_dim[1];

  if(grid_data->parallel == MPM_TRUE){
    bx0 = grid_data->box[0] - (haloX0 * hx);
    by0 = grid_data->box[1] - (haloY0 * hy);
  }else{
    bx0 = grid_data->box[0];
    by0 = grid_data->box[1];
  }

  for(index = 0; index < particle_count; index++){
    nx0 = particles[index].position[0];
    ny0 = particles[index].position[1];
    i = (int)((nx0 - bx0)/hx);
    j = (int)((ny0 - by0)/hy);

    particles[index].nodes[0].i = i;
    particles[index].nodes[0].j = j;
    particles[index].nodes[1].i = i+1;
    particles[index].nodes[1].j = j;
    particles[index].nodes[2].i = i+1;
    particles[index].nodes[2].j = j+1;
    particles[index].nodes[3].i = i;
    particles[index].nodes[3].j = j+1;

    // basis calculations
    wx = (nx0 - grid[i][j].position[0])/hx;
    wy = (ny0 - grid[i][j].position[1])/hy;
    dwx = 1/hx;
    dwy = 1/hy;

    particles[index].basis[0]   = (1-wx) * (1-wy);
    particles[index].dbasisX[0] = (-dwx) * (1-wy);
    particles[index].dbasisY[0] = (1-wx) * (-dwy);

    particles[index].basis[1]   = (wx)   * (1-wy);
    particles[index].dbasisX[1] = (dwx)  * (1-wy);
    particles[index].dbasisY[1] = (wx)   * (-dwy);

    particles[index].basis[2]   = (wx)   * (wy);
    particles[index].dbasisX[2] = (dwx)  * (wy);
    particles[index].dbasisY[2] = (wx)   * (dwy);

    particles[index].basis[3]   = (1-wx) * (wy);
    particles[index].dbasisX[3] = (-dwx) * (wy);
    particles[index].dbasisY[3] = (1-wx) * (dwy);

    //Mapping mass and momemtum to grid
    grid[i][j].has_mass     = MPM_TRUE;
    grid[i][j].mass        += particles[index].mass * particles[index].basis[0];
    grid[i][j].momentum[0] += particles[index].mass * particles[index].velocity[0] * particles[index].basis[0];
    grid[i][j].momentum[1] += particles[index].mass * particles[index].velocity[1] * particles[index].basis[0];

    grid[i+1][j].has_mass     = MPM_TRUE;
    grid[i+1][j].mass        += particles[index].mass * particles[index].basis[1];
    grid[i+1][j].momentum[0] += particles[index].mass * particles[index].velocity[0] * particles[index].basis[1];
    grid[i+1][j].momentum[1] += particles[index].mass * particles[index].velocity[1] * particles[index].basis[1];

    grid[i+1][j+1].has_mass     = MPM_TRUE;
    grid[i+1][j+1].mass        += particles[index].mass * particles[index].basis[2];
    grid[i+1][j+1].momentum[0] += particles[index].mass * particles[index].velocity[0] * particles[index].basis[2];
    grid[i+1][j+1].momentum[1] += particles[index].mass * particles[index].velocity[1] * particles[index].basis[2];

    grid[i][j+1].has_mass     = MPM_TRUE;
    grid[i][j+1].mass        += particles[index].mass * particles[index].basis[3];
    grid[i][j+1].momentum[0] += particles[index].mass * particles[index].velocity[0] * particles[index].basis[3];
    grid[i][j+1].momentum[1] += particles[index].mass * particles[index].velocity[1] * particles[index].basis[3];
  }
}

void momentumToVelocityOnGrid(GridData *grid_data, Node** grid)
{
  int i,j,index;
  int num_nodesX = grid_data->num_nodes[0];
  int num_nodesY = grid_data->num_nodes[1];

  for(i = 0; i < num_nodesX; i++){
    for(j = 0; j < num_nodesY; j++){
      if(grid[i][j].has_mass == MPM_TRUE){
        grid[i][j].velocity[0] = grid[i][j].momentum[0]/grid[i][j].mass;
        grid[i][j].velocity[1] = grid[i][j].momentum[1]/grid[i][j].mass;
      }
      if(grid_data->gridBoundaries[0] == BOUNDARY && i == 0){
        grid[i][j].velocity[0] = 0;
        grid[i][j].velocity[1] = 0;
      }
      if(grid_data->gridBoundaries[1] == BOUNDARY && j == 0){
        grid[i][j].velocity[0] = 0;
        grid[i][j].velocity[1] = 0;
      }
      if(grid_data->gridBoundaries[2] == BOUNDARY && i == num_nodesX-1){
        grid[i][j].velocity[0] = 0;
        grid[i][j].velocity[1] = 0;
      }
      if(grid_data->gridBoundaries[3] == BOUNDARY && j == num_nodesY-1){
        //grid[i][j].velocity[0] = 0;
        grid[i][j].velocity[0] = grid[i][j].momentum[0]/grid[i][j].mass;
        grid[i][j].velocity[1] = 0;
      }
      /**
      if(grid[i][j].has_mass == MPM_TRUE){
        if(grid_data->gridBoundaries[0] == BOUNDARY && i == 0){
          if(grid_data->gridBoundaryType[0] == DIRICHLET){
            grid[i][j].velocity[0] = 0;
            grid[i][j].velocity[1] = 0;
          }else if(grid_data->gridBoundaryType[0] == REFLEXIVE){
            grid[i][j].velocity[0] = -grid[i][j].momentum[0]/grid[i][j].mass;
            grid[i][j].velocity[1] = grid[i][j].momentum[1]/grid[i][j].mass;
          }else{
            grid[i][j].velocity[0] = 0;
            grid[i][j].velocity[1] = 0;
          }
        }else if(grid_data->gridBoundaries[1] == BOUNDARY && j == 0){
          if(grid_data->gridBoundaryType[0] == DIRICHLET){
            grid[i][j].velocity[0] = 0;
            grid[i][j].velocity[1] = 0;
          }else if(grid_data->gridBoundaryType[0] == REFLEXIVE){
            grid[i][j].velocity[0] = grid[i][j].momentum[0]/grid[i][j].mass;
            grid[i][j].velocity[1] = -grid[i][j].momentum[1]/grid[i][j].mass;
          }else{
            grid[i][j].velocity[0] = 0;
            grid[i][j].velocity[1] = 0;
          }
        }else if(grid_data->gridBoundaries[2] == BOUNDARY && i == num_nodesX-1){
          if(grid_data->gridBoundaryType[0] == DIRICHLET){
            grid[i][j].velocity[0] = 0;
            grid[i][j].velocity[1] = 0;
          }else if(grid_data->gridBoundaryType[0] == REFLEXIVE){
            grid[i][j].velocity[0] = -grid[i][j].momentum[0]/grid[i][j].mass;
            grid[i][j].velocity[1] = grid[i][j].momentum[1]/grid[i][j].mass;
          }else{
            grid[i][j].velocity[0] = 0;
            grid[i][j].velocity[1] = 0;
          }
        }else if(grid_data->gridBoundaries[3] == BOUNDARY && j == num_nodesY-1){
          if(grid_data->gridBoundaryType[0] == DIRICHLET){
            grid[i][j].velocity[0] = 0;
            grid[i][j].velocity[1] = 0;
          }else if(grid_data->gridBoundaryType[0] == REFLEXIVE){
            grid[i][j].velocity[0] = grid[i][j].momentum[0]/grid[i][j].mass;
            grid[i][j].velocity[1] = -grid[i][j].momentum[1]/grid[i][j].mass;
          }else{
            grid[i][j].velocity[0] = 0;
            grid[i][j].velocity[1] = 0;
          }
        }else{
          grid[i][j].velocity[0] = grid[i][j].momentum[0]/grid[i][j].mass;
          grid[i][j].velocity[1] = grid[i][j].momentum[1]/grid[i][j].mass;
        }
      }
      **/
    }
  }
}

void updateStressAndStrain(GridData *grid_data, Node** grid, Particle* particles, double dt){
  int i,j,index;
  int particle_count = particles[0].particle_count;
  int haloX0 = grid_data->halo_cells[0];
  int haloY0 = grid_data->halo_cells[1];
  int num_nodesX = grid_data->num_nodes[0];
  int num_nodesY = grid_data->num_nodes[1];
  double temp1[2][2];
  double temp2[2][2];
  double D[2][2];
  double mu, lam;

  for(index = 0; index < particle_count; index++){
    i = particles[index].nodes[0].i;
    j = particles[index].nodes[0].j;
    mu = particles[index].E/(1 + particles[index].poisson);
    lam = mu*(particles[index].poisson/(1 - 2 * particles[index].poisson));

    //Compute the velocity gradient
    particles[index].vel_grad[0][0]  = 0;
    particles[index].vel_grad[0][0] += grid[i][j].velocity[0]     * particles[index].dbasisX[0];
    particles[index].vel_grad[0][0] += grid[i+1][j].velocity[0]   * particles[index].dbasisX[1];
    particles[index].vel_grad[0][0] += grid[i+1][j+1].velocity[0] * particles[index].dbasisX[2];
    particles[index].vel_grad[0][0] += grid[i][j+1].velocity[0]   * particles[index].dbasisX[3];

    particles[index].vel_grad[0][1]  = 0;
    particles[index].vel_grad[0][1] += grid[i][j].velocity[0]     * particles[index].dbasisY[0];
    particles[index].vel_grad[0][1] += grid[i+1][j].velocity[0]   * particles[index].dbasisY[1];
    particles[index].vel_grad[0][1] += grid[i+1][j+1].velocity[0] * particles[index].dbasisY[2];
    particles[index].vel_grad[0][1] += grid[i][j+1].velocity[0]   * particles[index].dbasisY[3];

    particles[index].vel_grad[1][0]  = 0;
    particles[index].vel_grad[1][0] += grid[i][j].velocity[1]     * particles[index].dbasisX[0];
    particles[index].vel_grad[1][0] += grid[i+1][j].velocity[1]   * particles[index].dbasisX[1];
    particles[index].vel_grad[1][0] += grid[i+1][j+1].velocity[1] * particles[index].dbasisX[2];
    particles[index].vel_grad[1][0] += grid[i][j+1].velocity[1]   * particles[index].dbasisX[3];

    particles[index].vel_grad[1][1]  = 0;
    particles[index].vel_grad[1][1] += grid[i][j].velocity[1]     * particles[index].dbasisY[0];
    particles[index].vel_grad[1][1] += grid[i+1][j].velocity[1]   * particles[index].dbasisY[1];
    particles[index].vel_grad[1][1] += grid[i+1][j+1].velocity[1] * particles[index].dbasisY[2];
    particles[index].vel_grad[1][1] += grid[i][j+1].velocity[1]   * particles[index].dbasisY[3];

    //Update the deformation gradient
    temp1[0][0] = 1 + particles[index].vel_grad[0][0]*dt;
    temp1[0][1] = 0 + particles[index].vel_grad[0][1]*dt;
    temp1[1][0] = 0 + particles[index].vel_grad[1][0]*dt;
    temp1[1][1] = 1 + particles[index].vel_grad[1][1]*dt;

    temp2[0][0] = temp1[0][0] * particles[index].F[0][0] + temp1[0][1] * particles[index].F[1][0];
    temp2[0][1] = temp1[0][0] * particles[index].F[0][1] + temp1[0][1] * particles[index].F[1][1];
    temp2[1][0] = temp1[1][0] * particles[index].F[0][0] + temp1[1][1] * particles[index].F[1][0];
    temp2[1][1] = temp1[1][0] * particles[index].F[0][1] + temp1[1][1] * particles[index].F[1][1];

    particles[index].F[0][0] = temp2[0][0];
    particles[index].F[0][1] = temp2[0][1];
    particles[index].F[1][0] = temp2[1][0];
    particles[index].F[1][1] = temp2[1][1];

    D[0][0] =  particles[index].vel_grad[0][0];
    D[0][1] =  .5*(particles[index].vel_grad[0][1] + particles[index].vel_grad[1][0]);
    D[1][0] =  .5*(particles[index].vel_grad[1][0] + particles[index].vel_grad[0][1]);
    D[1][1] =  particles[index].vel_grad[1][1];

    particles[index].strain[0][0] += D[0][0] * dt;
    particles[index].strain[0][1] += D[0][1] * dt;
    particles[index].strain[1][0] += D[1][0] * dt;
    particles[index].strain[1][1] += D[1][1] * dt;

    particles[index].stress[0][0] = (mu * particles[index].strain[0][0]
                                  + lam * (particles[index].strain[0][0] + particles[index].strain[1][1]));
    particles[index].stress[0][1] = (mu * particles[index].strain[0][1]);
    particles[index].stress[1][0] = (mu * particles[index].strain[1][0]);
    particles[index].stress[1][1] = (mu * particles[index].strain[1][1]
                                  + lam * (particles[index].strain[0][0] + particles[index].strain[1][1]));

    /**
    particles[index].stress[0][0] = (mu * D[0][0] + lam * (D[0][0] + D[1][1])) * dt;
    particles[index].stress[0][1] = (mu * D[0][1]) * dt;
    particles[index].stress[1][0] = (mu * D[1][0]) * dt;
    particles[index].stress[1][1] = (mu * D[1][1] + lam * (D[0][0] + D[1][1])) * dt;
    **/
    //printf("dt: %f, vel_grad [0]: %f, [1]: %f, [2]: %f, [3]: %f\n", dt, particles[index].vel_grad[0][0], particles[index].vel_grad[0][1], particles[index].vel_grad[1][0], particles[index].vel_grad[1][1]);
  }


}

void computeForces(GridData *grid_data, Node **grid, Particle *particles)
{
  int i,j,index;
  int particle_count = particles[0].particle_count;
  int num_nodesX = grid_data->num_nodes[0];
  int num_nodesY = grid_data->num_nodes[1];

  for(index = 0; index < particle_count; index++){
    i = particles[index].nodes[0].i;
    j = particles[index].nodes[0].j;
    grid[i][j].int_force[0] -= particles[index].stress[0][0] * particles[index].dbasisX[0]
                             + particles[index].stress[0][1] * particles[index].dbasisY[0];
    grid[i][j].int_force[1] -= particles[index].stress[1][0] * particles[index].dbasisX[0]
                             + particles[index].stress[1][1] * particles[index].dbasisY[0];

    grid[i+1][j].int_force[0] -= particles[index].stress[0][0] * particles[index].dbasisX[1]
                               + particles[index].stress[0][1] * particles[index].dbasisY[1];
    grid[i+1][j].int_force[1] -= particles[index].stress[1][0] * particles[index].dbasisX[1]
                               + particles[index].stress[1][1] * particles[index].dbasisY[1];

    grid[i+1][j+1].int_force[0] -= particles[index].stress[0][0] * particles[index].dbasisX[2]
                                 + particles[index].stress[0][1] * particles[index].dbasisY[2];
    grid[i+1][j+1].int_force[1] -= particles[index].stress[1][0] * particles[index].dbasisX[2]
                                 + particles[index].stress[1][1] * particles[index].dbasisY[2];

    grid[i][j+1].int_force[0] -= particles[index].stress[0][0] * particles[index].dbasisX[3]
                               + particles[index].stress[0][1] * particles[index].dbasisY[3];
    grid[i][j+1].int_force[1] -= particles[index].stress[1][0] * particles[index].dbasisX[3]
                               + particles[index].stress[1][1] * particles[index].dbasisY[3];

    grid[i][j].body_force[0] += (particles[index].mass * grid_data->gravity[0]) * particles[index].basis[0];
    grid[i][j].body_force[1] += (particles[index].mass * grid_data->gravity[1]) * particles[index].basis[0];

    grid[i+1][j].body_force[0] += (particles[index].mass * grid_data->gravity[0]) * particles[index].basis[1];
    grid[i+1][j].body_force[1] += (particles[index].mass * grid_data->gravity[1]) * particles[index].basis[1];

    grid[i+1][j+1].body_force[0] += (particles[index].mass * grid_data->gravity[0]) * particles[index].basis[2];
    grid[i+1][j+1].body_force[1] += (particles[index].mass * grid_data->gravity[1]) * particles[index].basis[2];

    grid[i][j+1].body_force[0] += (particles[index].mass * grid_data->gravity[0]) * particles[index].basis[3];
    grid[i][j+1].body_force[1] += (particles[index].mass * grid_data->gravity[1]) * particles[index].basis[3];
  }
}

void computeAcceleration(GridData *grid_data, Node **grid)
{
  int i,j;
  int num_nodesX = grid_data->num_nodes[0];
  int num_nodesY = grid_data->num_nodes[1];
  for(i = 0; i < num_nodesX; i++){
    for(j = 0; j < num_nodesY; j++){
      if(grid[i][j].has_mass == MPM_TRUE){
        grid[i][j].acceleration[0] = (grid[i][j].int_force[0] + grid[i][j].body_force[0])/grid[i][j].mass;
        grid[i][j].acceleration[1] = (grid[i][j].int_force[1] + grid[i][j].body_force[1])/grid[i][j].mass;
      }
      if(grid_data->gridBoundaries[0] == BOUNDARY && i == 0){
        grid[i][j].acceleration[0] = 0;
        grid[i][j].acceleration[1] = 0;
      }
      if(grid_data->gridBoundaries[1] == BOUNDARY && j == 0){
        grid[i][j].acceleration[0] = 0;
        grid[i][j].acceleration[1] = 0;
      }
      if(grid_data->gridBoundaries[2] == BOUNDARY && i == num_nodesX-1){
        grid[i][j].acceleration[0] = 0;
        grid[i][j].acceleration[1] = 0;
      }
      if(grid_data->gridBoundaries[3] == BOUNDARY && j == num_nodesY-1){
        //grid[i][j].acceleration[0] = 0;
        grid[i][j].acceleration[0] = (grid[i][j].int_force[0] + grid[i][j].body_force[0])/grid[i][j].mass;
        grid[i][j].acceleration[1] = 0;
      }
    }
  }
}

void updateNodalValues(GridData *grid_data, Node **grid, double dt)
{
  int i,j,index;
  int num_nodesX = grid_data->num_nodes[0];
  int num_nodesY = grid_data->num_nodes[1];
  for(i = 0; i < num_nodesX; i++){
    for(j = 0; j < num_nodesY; j++){
      grid[i][j].velocity[0] += grid[i][j].acceleration[0] * dt;
      grid[i][j].velocity[1] += grid[i][j].acceleration[1] * dt;

      if(grid_data->gridBoundaries[0] == BOUNDARY && i == 0){
        grid[i][j].velocity[0] = 0;
        grid[i][j].velocity[1] = 0;
      }
      if(grid_data->gridBoundaries[1] == BOUNDARY && j == 0){
        grid[i][j].velocity[0] = 0;
        grid[i][j].velocity[1] = 0;
      }
      if(grid_data->gridBoundaries[2] == BOUNDARY && i == num_nodesX-1){
        grid[i][j].velocity[0] = 0;
        grid[i][j].velocity[1] = 0;
      }
      if(grid_data->gridBoundaries[3] == BOUNDARY && j == num_nodesY-1){
        //grid[i][j].velocity[0] += 0;
        grid[i][j].velocity[0] += grid[i][j].acceleration[0] * dt;
        grid[i][j].velocity[1] = 0;
      }
    }
  }
}

void updateParticlePosition(GridData *grid_data, Node **grid, Particle *particles, double dt)
{
  int i,j,index;
  int particle_count = particles[0].particle_count;
  int num_nodesX = grid_data->num_nodes[0];
  int num_nodesY = grid_data->num_nodes[1];

  for(index = 0; index < particle_count; index++){
    i = particles[index].nodes[0].i;
    j = particles[index].nodes[0].j;

    particles[index].velocity[0] += (grid[i][j].acceleration[0]     * particles[index].basis[0]
                                   + grid[i+1][j].acceleration[0]   * particles[index].basis[1]
                                   + grid[i+1][j+1].acceleration[0] * particles[index].basis[2]
                                   + grid[i][j+1].acceleration[0]   * particles[index].basis[3]) * dt;

    particles[index].velocity[1] += (grid[i][j].acceleration[1]     * particles[index].basis[0]
                                   + grid[i+1][j].acceleration[1]   * particles[index].basis[1]
                                   + grid[i+1][j+1].acceleration[1] * particles[index].basis[2]
                                   + grid[i][j+1].acceleration[1]   * particles[index].basis[3]) * dt;

    particles[index].position[0] += (grid[i][j].velocity[0]     * particles[index].basis[0]
                                   + grid[i+1][j].velocity[0]   * particles[index].basis[1]
                                   + grid[i+1][j+1].velocity[0] * particles[index].basis[2]
                                   + grid[i][j+1].velocity[0]   * particles[index].basis[3]) * dt;

    particles[index].position[1] += (grid[i][j].velocity[1]     * particles[index].basis[0]
                                   + grid[i+1][j].velocity[1]   * particles[index].basis[1]
                                   + grid[i+1][j+1].velocity[1] * particles[index].basis[2]
                                   + grid[i][j+1].velocity[1]   * particles[index].basis[3]) * dt;

  }
}

void pUpdateParticlePosition(GridData *grid_data, Node **grid, Particle *particles,
                             Particle *post_particles, Particle **shalo_parts, double dt)
{
  int i,j,index,pindex;
  int rank;
  int particle_count = particles[0].particle_count;
  int pre_pc = particle_count;
  int post_pc = 0;
  int num_nodesX = grid_data->num_nodes[0];
  int num_nodesY = grid_data->num_nodes[1];
  double px,py;
  MpmBool particle_removed;
  Particle* pholder;


  for(i = 0; i < 8; i++){
    rank = grid_data->rank[i];
    if(rank >= 0){
      shalo_parts[i][0].particle_count = 0;
    }
  } 

  post_pc = 0;
  for(index = 0; index < particle_count; index++){
    particle_removed = MPM_FALSE;
    i = particles[index].nodes[0].i;
    j = particles[index].nodes[0].j;

    particles[index].velocity[0] += (grid[i][j].acceleration[0]     * particles[index].basis[0]
                                   + grid[i+1][j].acceleration[0]   * particles[index].basis[1]
                                   + grid[i+1][j+1].acceleration[0] * particles[index].basis[2]
                                   + grid[i][j+1].acceleration[0]   * particles[index].basis[3]) * dt;

    particles[index].velocity[1] += (grid[i][j].acceleration[1]     * particles[index].basis[0]
                                   + grid[i+1][j].acceleration[1]   * particles[index].basis[1]
                                   + grid[i+1][j+1].acceleration[1] * particles[index].basis[2]
                                   + grid[i][j+1].acceleration[1]   * particles[index].basis[3]) * dt;

    particles[index].position[0] += (grid[i][j].velocity[0]     * particles[index].basis[0]
                                   + grid[i+1][j].velocity[0]   * particles[index].basis[1]
                                   + grid[i+1][j+1].velocity[0] * particles[index].basis[2]
                                   + grid[i][j+1].velocity[0]   * particles[index].basis[3]) * dt;

    particles[index].position[1] += (grid[i][j].velocity[1]     * particles[index].basis[0]
                                   + grid[i+1][j].velocity[1]   * particles[index].basis[1]
                                   + grid[i+1][j+1].velocity[1] * particles[index].basis[2]
                                   + grid[i][j+1].velocity[1]   * particles[index].basis[3]) * dt;
    
    
    px = particles[index].position[0];
    py = particles[index].position[1];

    if(px < grid_data->box[0] && grid_data->gridBoundaries[0] == HALO){
      if(py < grid_data->box[1] && grid_data->gridBoundaries[1] == HALO){
        pindex = shalo_parts[1][0].particle_count;
        shalo_parts[1][pindex] = particles[index];
        shalo_parts[1][0].particle_count = pindex + 1;
        particle_removed = MPM_TRUE;
      }else if(py > grid_data->box[3] && grid_data->gridBoundaries[3] == HALO){
        pindex = shalo_parts[7][0].particle_count;
        shalo_parts[7][pindex] = particles[index];
        shalo_parts[7][0].particle_count = pindex + 1;
        particle_removed = MPM_TRUE;
      }else{
        pindex = shalo_parts[0][0].particle_count;
        shalo_parts[0][pindex] = particles[index];
        shalo_parts[0][0].particle_count = pindex + 1;
        particle_removed = MPM_TRUE;
      }
    }else if(px > grid_data->box[2] && grid_data->gridBoundaries[2] == HALO){
      if(py < grid_data->box[1] && grid_data->gridBoundaries[1] == HALO){
        pindex = shalo_parts[3][0].particle_count;
        shalo_parts[3][pindex] = particles[index];
        shalo_parts[3][0].particle_count = pindex + 1;
        particle_removed = MPM_TRUE;
      }else if(py > grid_data->box[3] && grid_data->gridBoundaries[3] == HALO){
        pindex = shalo_parts[5][0].particle_count;
        shalo_parts[5][pindex] = particles[index];
        shalo_parts[5][0].particle_count = pindex + 1;
        particle_removed = MPM_TRUE;
      }else{
        pindex = shalo_parts[4][0].particle_count;
        shalo_parts[4][pindex] = particles[index];
        shalo_parts[4][0].particle_count = pindex + 1;
        particle_removed = MPM_TRUE;
      }
    }

    if(py < grid_data->box[1] && grid_data->gridBoundaries[1] == HALO){
      if(px >= grid_data->box[0] && px <= grid_data->box[2]){
        pindex = shalo_parts[2][0].particle_count;
        shalo_parts[2][pindex] = particles[index];
        shalo_parts[2][0].particle_count = pindex + 1;
        particle_removed = MPM_TRUE;
      }
    }else if(py > grid_data->box[3] && grid_data->gridBoundaries[3] == HALO){
      if(px >= grid_data->box[0] && px <= grid_data->box[2]){
        pindex = shalo_parts[6][0].particle_count;
        shalo_parts[6][pindex] = particles[index];
        shalo_parts[6][0].particle_count = pindex + 1;
        particle_removed = MPM_TRUE;
      }
    }
    
    if(particle_removed == MPM_FALSE){
      post_particles[post_pc] = particles[index];
      post_pc++;
    }else{
      printf("tparticle moved: %f, %f\n", px, py);
      printf("b0: %f, b1: %f, b2: %f, b3: %f\n", grid_data->box[0], grid_data->box[1], grid_data->box[2], grid_data->box[3]);
    }
  }
  post_particles[0].particle_count = post_pc;

  pholder = particles;
  particles = post_particles;
  post_particles = pholder;
}
