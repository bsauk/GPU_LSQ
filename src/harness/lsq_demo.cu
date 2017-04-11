#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>

// Will need to include some type of header file with all of my other functions defined in it
// Need to accept an input file
int main(int argc, char** argv) {
  // Error handling for # of inputs
  if(argc == 1) {
    std::cout << "Please provide .dat file and number of rows and cols!" << std::endl;
    return 0;
  }
  if(argc == 2) {
    std::cout << "Please provide number of rows and cols!" << std::endl;
    return 0;
  }
  if(argc == 3){
    std::cout << "Please provide number of cols!" << std::endl;
    return 0;
  }
  if(argc > 4) {
    std:: cout << "Too many arguments please only provide file name, rows, cols" << std::endl;
  }
  
  int rows, cols;
  rows = strtol(argv[2], NULL, 10);
  cols = strtol(argv[3], NULL, 10);
  double A[rows][cols];

  std::string filename(argv[1]);
  std::ifstream file;
  file.open(filename);
  file >> rows >> cols;
  for(int i=0; i<rows; i++) {
    for(int j=0; j<cols; j++) {
      file >> A[i][j];
    }
  }
 
}
