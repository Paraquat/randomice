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
}

Hbond& Hbond::operator= (const Hbond& h)
{
  if (this == &h) return *this;
  O1 = h.O1;
  O2 = h.O2;
  H1 = h.H1;
  H2 = h.H2;
  return *this;
}

Hbond::Hbond(Atom::atom_ptr o1, Atom::atom_ptr o2, Atom::atom_ptr h1, Atom::atom_ptr h2)
{
  O1 = o1;
  O2 = o2;
  H1 = h1;
  H2 = h2;
}

std::ostream& operator<< (std::ostream& os, Hbond& hb)
{
  os  << *hb.O1 << std::endl
      << *hb.O2 << std::endl
      << *hb.H1 << std::endl
      << *hb.H2 << std::endl;
  return os;
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