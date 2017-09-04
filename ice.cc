#include "ice.h"

// Default constructor
Ice::Ice()
{
}

// Destructor
Ice::~Ice()
{
}

Ice::Ice(Cell& cell)
{
  lat = cell.lat;
  atoms = cell.atoms;
  natoms = cell.natoms;
  frac = cell.frac; 
}

// Compute all hydrogen positions for a given oxygen lattice
void Ice::get_h_pos(void)
{
  Eigen::Vector3d oo, r_h;
  std::string species = "H";
  Atom a;

  if (frac == true) frac2cart_all();
  get_nn(oo_max);
  for (int i=0; i<natoms; i++){
    for (int j=0; j<atoms[i].nn.size(); j++){
      oo = mic_cart(atoms[i], *atoms[i].nn[j]);
      // oo = (*atoms[i].nn[j]).r - atoms[i].r;
      r_h = atoms[i].r + oh_default*oo.normalized();
      a = Atom(species, natoms+1, r_h, false);
      add_atom(a);
    }
  }
}