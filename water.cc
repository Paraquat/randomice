#include "water.h"

Water::Water()
{
}

Water::~Water()
{
}

Water::Water(const Water& w)
{
  O = w.O;
  H1 = w.H1;
  H2 = w.H2;
  H3 = w.H3;
  H4 = w.H4;
  label = w.label;
}

Water& Water::operator= (const Water& w)
{
  if (this == &w) return *this;
  O = w.O;
  H1 = w.H1;
  H2 = w.H2;
  H3 = w.H3;
  H4 = w.H4;
  label = w.label;
  return *this;
}

bool Water::operator== (const Water& w)
{
  if (label == w.label) return true;
  else return false;
}

Water::Water(int l, Atom::atom_ptr o, Atom::atom_ptr h1, Atom::atom_ptr h2, Atom::atom_ptr h3, Atom::atom_ptr h4)
{
  label = l;
  O = o;
  H1 = h1;
  H2 = h2;
  H3 = h3;
  H4 = h4;
}

std::ostream& operator<< (std::ostream& os, Water& w)
{
  os  << *w.O << std::endl
      << *w.H1 << std::endl
      << *w.H2 << std::endl
      << *w.H3 << std::endl
      << *w.H4 << std::endl;
  return os;
}

std::ofstream& operator<< (std::ofstream& ofs, Water& w)
{
  ofs  << *w.O << std::endl;
  if (w.H1->occupied) ofs << *w.H1 << std::endl;
  if (w.H2->occupied) ofs << *w.H2 << std::endl;
  if (w.H3->occupied) ofs << *w.H3 << std::endl;
  if (w.H4->occupied) ofs << *w.H4 << std::endl;
  return ofs;
}

void Water::add_nn(water_ptr n)
{
  nn.push_back(n);
}

void Water::add_hbond(hbond_ptr hb)
{
  hbonds.push_back(hb);
}

void Water::check_defect(void)
{
  nH = 0;
  if (H1 -> occupied) nH++;
  if (H2 -> occupied) nH++;
  if (H3 -> occupied) nH++;
  if (H4 -> occupied) nH++;
  if (nH == 0){
    ionic = "O2-";
  } else if (nH == 1){
    ionic = "OH-";
  } else if (nH == 2){
    ionic = "None";
  } else if (nH == 3){
    ionic = "H3O+";
  } else if (nH == 4){
    ionic = "H4O2+";
  }
}

int Water::coord(void)
{
  int coordination = 0;
  if (H1->occupied) coordination++;
  if (H2->occupied) coordination++;
  if (H3->occupied) coordination++;
  if (H4->occupied) coordination++;
  return coordination;
}