#include <stdlib.h>
#include <chrono>
#include <iostream>
#include <math.h>
#include "SPP.h"
using namespace std;

int main(){
  SPP bio(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, false);
  bio.add_cell(0, 0, 0, 0);
  bio.add_cell(2.1, 0, 0, 3.1415);
  for (int i=0; i<100; i++){
    printf("los centros de 1 son: %f y %f\n", bio.cell_centers[0], bio.cell_centers[1]);
    printf("los centros de 2 son: %f y %f\n", bio.cell_centers[3], bio.cell_centers[4]);
    printf("Esto debería ser 0: %f, %f\n", bio.cell_centers[2], bio.cell_centers[5]);
    bio.step();
  }
  return 0;
}