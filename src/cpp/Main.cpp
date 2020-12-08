
#include "iS3D.h"

int main(int argc, char *argv[])
{
  //create an instance of IS3D class
  IS3D particlization;

  //run iS3D
  //if argument == 1, freeeout surface is read from file
  //otherwise freezeout surface is read from memory
  particlization.run_particlization(1);


  // reading freezeout info from file or memory depends on whether any command arguments are given
  // if(argc == 1)
  // {
  // 	particlization.run_particlization(0);
  // }
  // else
  // {
  // 	particlization.run_particlization(1);
  // }
}
