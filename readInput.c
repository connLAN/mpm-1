/**
 *
 * readInput.c
 *
 */

#include<stdio.h>
#include<stdlib.h>
#include<mpm_structs.h>

void readFile(char* filename, InputData* data){

  FILE* input;
  input = fopen(filename, "r");
  
  fclose(input);
}

int main(int argc, char* argc){

  return 0;
}
