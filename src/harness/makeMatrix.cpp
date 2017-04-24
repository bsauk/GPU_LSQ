#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>


int main(int argc, char* argv[]) {
  int m = atoi(argv[1]);
  int n = atoi(argv[2]);
  int seeds[5] = {79, 6153, 1641, 9551, 20770};
  int in = atoi(argv[3]);
  double div = RAND_MAX/1000;
  srand((unsigned)seeds[in]);
  for(int row=0; row<m; row++) {
    for(int col=0; col<n; col++) {
      std::cout << rand()/div << "   ";
    }
    std::cout << 1.00 << "   " << rand()/div << "   " << std::endl;
  }
  return 0;
}
