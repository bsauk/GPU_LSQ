#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include "lsq.h"
#include "sub.h"

#include "CycleTimer.h"

void compare_results(int first, int max_size, int nbest, int lopt_dim1, double** ressGold, double** ressGPU, int** loptGold, int** loptGPU) {
  bool ressTest = true;
  bool loptTest = true;
  for(int i=0; i<max_size; i++) {
    for(int j=0; j<nbest; j++) {
      if(ressGold[i][j] != ressGPU[i][j]) {
	std::cout << "INCORRECT OUTPUT: ressGold[" << i << "][" << j << "] = " << ressGold[i][j] << "; ressGPU[" << i << "][" << j << "] = " << ressGPU[i][j] << std::endl;
	ressTest = false;
      } 
    }
  }
  if(ressTest) std::cout << "ress is correct!" << std::endl;
  
  for(int i=0; i<nbest; i++) {
    for(int j=0; j<lopt_dim1; j++) {
      if(loptGold[i][j] != loptGPU[i][j]) {
	std::cout << "INCORRECT OUTPUT: loptGold[" << i << "][" << j << "] = " << loptGold[i][j] << "; loptGPU[" << i << "][" << j << "] = " << loptGPU[i][j] << std::endl;
	loptTest = false;
      }
    }
  }
  if(loptTest) std::cout << "lopt is correct!" << std::endl;
}
