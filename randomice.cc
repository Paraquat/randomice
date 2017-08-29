#include "randomice.h"

int main(int argc, char* argv[]){

  Cell cell;

  cell.read_cell("MgSiO3.cell");
  std::cout << cell;

  return 0;
}
