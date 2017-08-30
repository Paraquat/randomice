#include "randomice.h"

int main(int argc, char* argv[]){

  Cell cell;

  cell.read_cell("MgSiO3.cell");
  cell.frac2cart_all();
  std::cout << cell;

  return 0;
}
