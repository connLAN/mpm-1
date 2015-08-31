/**
 *
 * mpmio.c
 *
 **/

#include <mpmio.h>

void writeParticlesVtkFile(GridData *grid_data, Node **grid, Particle *particles, const char* filename)
{
	FILE * outfile;
    outfile = fopen(filename, "w");
    writeHeader(outfile, "particle data", ASCII);
    writeParticleData(outfile, grid_data, grid, particles);
    fclose(outfile);
}

void writeGridVtkFile(GridData *grid_data, Node **grid, const char* filename){
	FILE * outfile;
    outfile = fopen(filename, "w");
    writeHeader(outfile, "grid data", ASCII);
    writeGridData(outfile, grid_data, grid);
    fclose(outfile);
}

void writeHeader(FILE * outfile, const char * title, int datatype )
{
    fprintf(outfile, "# vtk DataFile Version 2.0\n");
    fprintf(outfile, "%s\n", title);
    if(datatype == ASCII)
        fprintf(outfile, "ASCII\n");
    else
        fprintf(outfile, "BINARY\n");
}

void writeParticleData(FILE * outfile, GridData *grid_data, Node **grid, Particle *particles)
{
	int i;
    int particle_count = particles[0].particle_count;
   
    fprintf(outfile, "\n");
    fprintf(outfile, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(outfile, "POINTS %d float\n", particle_count);
    for(i = 0; i < particle_count; i++){
        fprintf(outfile, "%20.12e %20.12e %20.12e\n", particles[i].position[0], particles[i].position[1], 0.0);
    }

    fprintf(outfile, "\n");
    fprintf(outfile, "CELLS %d %d\n", particle_count, 2*particle_count);
    for(i = 0; i < particle_count; i++){
        fprintf(outfile, "1 %d\n", i);
    }

    fprintf(outfile, "\n");
    fprintf(outfile, "CELL_TYPES %d\n", particle_count);
    for(i = 0; i < particle_count; i++){
        fprintf(outfile, "1\n");
    }

    fprintf(outfile, "\n");
    fprintf(outfile, "POINT_DATA %d\n", particle_count);
    fprintf(outfile, "SCALARS mass float\n");
    fprintf(outfile, "LOOKUP_TABLE default\n");
    for(i = 0; i < particle_count; i++){
        fprintf(outfile, "%20.12e\n", particles[i].mass);
    }

    fprintf(outfile, "VECTORS velocity float\n");
    for(i = 0; i < particle_count; i++){
        fprintf(outfile, "%20.12e %20.12e %20.12e\n", particles[i].velocity[0], particles[i].velocity[1], 0.0);
    }
}

void writeGridData(FILE * outfile, GridData *grid_data, Node **grid)
{
	int i,j,index;
    int grid_countI = grid_data->num_nodes[0];
    int grid_countJ = grid_data->num_nodes[1];
    int node_count = grid_countI * grid_countJ;
   
    fprintf(outfile, "\n");
    fprintf(outfile, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(outfile, "POINTS %d float\n", node_count);
    for(i = 0; i < grid_countI; i++){
      for(j = 0; j < grid_countJ; j++){
        fprintf(outfile, "%20.12e %20.12e %20.12e\n", grid[i][j].position[0], grid[i][j].position[1], 0.0);
      }
    }

    index = 0;
    fprintf(outfile, "\n");
    fprintf(outfile, "CELLS %d %d\n", node_count, 2*node_count);
    for(i = 0; i < grid_countI; i++){
      for(j = 0; j < grid_countJ; j++){
      	fprintf(outfile, "1 %d\n", index);
        index++;
      }
    }

    fprintf(outfile, "\n");
    fprintf(outfile, "CELL_TYPES %d\n", node_count);
    for(i = 0; i < grid_countI; i++){
      for(j = 0; j < grid_countJ; j++){
      	fprintf(outfile, "1\n");
      }
    }

    fprintf(outfile, "\n");
    fprintf(outfile, "POINT_DATA %d\n", node_count);
    fprintf(outfile, "SCALARS mass float\n");
    fprintf(outfile, "LOOKUP_TABLE default\n");
    for(i = 0; i < grid_countI; i++){
      for(j = 0; j < grid_countJ; j++){
      	fprintf(outfile, "%20.12e\n", grid[i][j].mass);
      }
    }

    fprintf(outfile, "VECTORS velocity float\n");
    for(i = 0; i < grid_countI; i++){
      for(j = 0; j < grid_countJ; j++){
      	fprintf(outfile, "%20.12e %20.12e %20.12e\n", grid[i][j].velocity[0], grid[i][j].velocity[1], 0.0);
      }
    }

    fprintf(outfile, "VECTORS acceleration float\n");
    for(i = 0; i < grid_countI; i++){
      for(j = 0; j < grid_countJ; j++){
      	fprintf(outfile, "%20.12e %20.12e %20.12e\n", grid[i][j].acceleration[0], grid[i][j].acceleration[1], 0.0);
      }
    }
}
