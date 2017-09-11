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
  bool occ = false;
  Atom a;

  if (frac == true) frac2cart_all();
  get_nn(oo_max);
  for (int i=0; i<natoms; i++){
    for (int j=0; j<atoms[i].nn.size(); j++){
      oo = mic_cart(atoms[i], *atoms[i].nn[j]);
      r_h = atoms[i].r + oh_default*oo.normalized();
      a = Atom(species, natoms+1, r_h, occ);
      std::cout << a.occupied << std::endl;
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
  nhbond = 0;
  get_nn(oo_max);
  for (int i=0; i<nwater; i++){
    Water::water_ptr w1(new Water(waters[i]));
    assert(waters[i].nn.size() == 4);
    for (int j=0; j<4; j++){
      Water::water_ptr w2(new Water(*waters[i].nn[j]));
      Atom::atom_ptr h1, h2;

      if (isPointOnLine(*w1->O, *w2->O, *w1->H1)){
        h1 = w1->H1;
      }
      else if (isPointOnLine(*w1->O, *w2->O, *w1->H2)){
        h1 = w1->H2;
      }
      else if (isPointOnLine(*w1->O, *w2->O, *w1->H3)){
        h1 = w1->H3;
      }
      else if (isPointOnLine(*w1->O, *w2->O, *w1->H4)){
        h1 = w1->H4;
      }
      else std::cout << "Failed to find h1 for water " << w1->label 
                     << std::endl;

      if (isPointOnLine(*w1->O, *w2->O, *w2->H1)){
        h2 = w2->H1;
      }
      else if (isPointOnLine(*w1->O, *w2->O, *w2->H2)){
        h2 = w2->H2;
      }
      else if (isPointOnLine(*w1->O, *w2->O, *w2->H3)){
        h2 = w2->H3;
      }
      else if (isPointOnLine(*w1->O, *w2->O, *w2->H4)){
        h2 = w2->H4;
      }
      else std::cout << "Failed to find h2 for water " << w2->label 
                     << std::endl;

      Hbond hb((*w1).O, (*w2).O , h1, h2);

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