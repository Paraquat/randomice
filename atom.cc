#include "atom.h"

Atom::Atom()
{
}

Atom::~Atom()
{
}

Atom::Atom(std::string symbol, int i, double x, double y, double z)
{
  name = symbol;
  label = i;
  r = Eigen::Vector3d(x, y, z);
  occupied = true;
}

Atom::Atom(std::string symbol, int i, Eigen::Vector3d v)
{
  name = symbol;
  label = i;
  r = v;
  occupied = true;
}

// Allow addition of an unoccupied atomic site
Atom::Atom(std::string symbol, int i, Eigen::Vector3d v, bool occ)
{
  name = symbol;
  label = i;
  r = v;
  occupied = occ;
}

Atom::Atom(const Atom& a)
{
  name = a.name;
  label = a.label;
  r = a.r;
  occupied = true;
}

bool Atom::operator== (const Atom& a)
{
  if (label == a.label) return true;
  else return false;
}

Atom& Atom::operator= (const Atom& a)
{
  if (this == &a) return *this;
  name = a.name;
  label = a.label;
  r = a.r;
  return *this;
}

std::ostream& operator<< (std::ostream& os, Atom& a)
{
  std::ostringstream oss;

  oss << a.name;
  os  << std::fixed << std::setprecision(8) \
      << std::left << std::setw(4) << oss.str() \
      << std::right << std::setw(20) << a.r[0] \
      << std::right << std::setw(20) << a.r[1] \
      << std::right << std::setw(20) << a.r[2];
  return os;
}

std::ofstream& operator<< (std::ofstream& ofs, Atom& a)
{
  std::ostringstream oss;

  oss << a.name;
  ofs << std::fixed << std::setprecision(8) \
      << std::left << std::setw(4) << oss.str()
      << std::right << std::setw(20) << a.r[0] \
      << std::right << std::setw(20) << a.r[1] \
      << std::right << std::setw(20) << a.r[2];
  return ofs;
}

void Atom::add_nn(atom_ptr n)
{
  nn.push_back(n);
}

void Atom::print_nn(void)
{
  std::cout << label << ": ";
  for (int i=0; i<nn.size(); i++){
    std::cout << (*nn[i]).label << " ";
  }
  std::cout << std::endl;
}