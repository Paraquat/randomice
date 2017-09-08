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
      r_h = atoms[i].r + oh_default*oo.normalized();
      a = Atom(species, natoms+1, r_h, false);
      add_atom(a);
    }
  }
}

void Ice::add_water(Water& w)
{
  waters.push_back(w);
  nwater++;
}

void Ice::add_hbond(Hbond& h)
{
  hbonds.push_back(h);
  nhbond++;
}

void Ice::get_waters(void)
{
  nwater = 0;
  get_nn(oh_max);
  for (int i=0; i<natoms; i++){
    if (atoms[i].name == "O"){
      assert(atoms[i].nn.size() == 4);
      Atom::atom_ptr o(new Atom(atoms[i]));
      Atom::atom_ptr h1(new Atom(*atoms[i].nn[0]));
      Atom::atom_ptr h2(new Atom(*atoms[i].nn[1]));
      Atom::atom_ptr h3(new Atom(*atoms[i].nn[2]));
      Atom::atom_ptr h4(new Atom(*atoms[i].nn[3]));

      Water w(nwater, o, h1, h2, h3, h4);
      add_water(w);
      // std::cout << waters.back() << std::endl;
    }
  }
}

void Ice::get_water_nn(double cutoff)
{
  Eigen::Vector3d d;

  assert(frac == false);
  for (int i=0; i<nwater; i++){
    for (int j=0; j<nwater; j++){
      if (i != j){
        d = mic_cart(*waters[i].O, *waters[j].O);
        if (d.norm() < cutoff){
          Water::water_ptr ptr(new Water(waters[j]));
          waters[i].add_nn(ptr);
        }
      }
    }
  }
}

void Ice::get_hbonds(void)
{
  Atom::atom_ptr h1, h2;
  nhbond = 0;
  get_nn(oo_max);
  for (int i=0; i<nwater; i++){
    for (int j=i+1; j<nwater; j++){
      if (isPointOnLine(*waters[i].O, *waters[j].O, *waters[i].H1)){
        h1.reset(new Atom(*waters[i].H1));
      }
      else if (isPointOnLine(*waters[i].O, *waters[j].O, *waters[i].H2)){
        h1.reset(new Atom(*waters[i].H2));
      }
      else if (isPointOnLine(*waters[i].O, *waters[j].O, *waters[i].H3)){
        h1.reset(new Atom(*waters[i].H3));
      }
      else if(isPointOnLine(*waters[i].O, *waters[j].O, *waters[i].H4)){
        h1.reset(new Atom(*waters[i].H4));
      }
      // else std::cout << "Failed to find h1" << std::endl;

      if(isPointOnLine(*waters[i].O, *waters[j].O, *waters[j].H1)){
        h2.reset(new Atom(*waters[j].H1));
      }
      else if(isPointOnLine(*waters[i].O, *waters[j].O, *waters[j].H2)){
        h2.reset(new Atom(*waters[j].H2));
      }
      else if(isPointOnLine(*waters[i].O, *waters[j].O, *waters[j].H3)){
        h2.reset(new Atom(*waters[j].H3));
      }
      else if(isPointOnLine(*waters[i].O, *waters[j].O, *waters[j].H4)){
        h2.reset(new Atom(*waters[j].H4));
      }
      // else std::cout << "Failed to find h2" << std::endl;
      Hbond hb(waters[i].O, waters[j].O, h1, h2);

      add_hbond(hb);
    }
  }
}

void Ice::print_ice(void)
{
  std::cout << nwater << " water molecules" << std::endl << std::endl;
  for (int i=0; i<nwater; i++){
    std::cout << *waters[i].O << std::endl;
    std::cout << *waters[i].H1 << std::endl;
    std::cout << *waters[i].H2 << std::endl;
    std::cout << *waters[i].H3 << std::endl;
    std::cout << *waters[i].H4 << std::endl;
    std::cout << std::endl;
  }
}