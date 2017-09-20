#include "hbond.h"

Hbond::Hbond()
{
}

Hbond::~Hbond()
{
}

Hbond::Hbond(const Hbond& h)
{
  O1 = h.O1;
  O2 = h.O2;
  H1 = h.H1;
  H2 = h.H2;
  W1 = h.W1;
  W2 = h.W2;
  w = h.w;
}

Hbond& Hbond::operator= (const Hbond& h)
{
  if (this == &h) return *this;
  O1 = h.O1;
  O2 = h.O2;
  H1 = h.H1;
  H2 = h.H2;
  W1 = h.W1;
  W2 = h.W2;
  w = h.w;
  return *this;
}

Hbond::Hbond(Atom::atom_ptr o1, Atom::atom_ptr o2, Atom::atom_ptr h1, Atom::atom_ptr h2, water_ptr w1, water_ptr w2)
{
  O1 = o1;
  O2 = o2;
  H1 = h1;
  H2 = h2;
  W1 = w1;
  W2 = w2;
  w.insert(O1->label);
  w.insert(O2->label);
}

std::ostream& operator<< (std::ostream& os, Hbond& hb)
{
  os  << *hb.O1 << std::endl
      << *hb.O2 << std::endl
      << *hb.H1 << std::endl
      << *hb.H2 << std::endl;
  // if (hb.H1->occupied) std::cout << *hb.H1 << std::endl;
  // if (hb.H2->occupied) std::cout << *hb.H2 << std::endl;
  return os;
}

// Switch the occupied and unoccupied hydrogen positions
void Hbond::swap_h(void)
{
  H1->occupied = !H1->occupied;
  H2->occupied = !H2->occupied;
}

water_ptr Hbond::get_target(void)
{
  if (H1->occupied) return W2;
  else return W1;
}

void Hbond::check_defect(void)
{
  if (H1 -> occupied && H2 -> occupied){
    bjerrum = "D";
  } else if (!(H1 -> occupied) && !(H2 -> occupied)){
    bjerrum = "L";
  } else {
    bjerrum = "None";
  }
}