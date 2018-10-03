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
  remove = false;
  nneighbour = 0;
  comment = "";
}

Atom::Atom(std::string symbol, int i, Eigen::Vector3d v)
{
  name = symbol;
  label = i;
  r = v;
  occupied = true;
  remove = false;
  nneighbour = 0;
  nneighbour = 0;
  comment = "";
}

// Allow addition of an unoccupied atomic site
Atom::Atom(std::string symbol, int i, Eigen::Vector3d v, bool occ)
{
  name = symbol;
  label = i;
  r = v;
  occupied = occ;
  remove = false;
  nneighbour = 0;
  nneighbour = 0;
  comment = "";
}

Atom::Atom(const Atom& a)
{
  // std::cout << "copy constructor invoked" << std::endl;
  name = a.name;
  label = a.label;
  r = a.r;
  occupied = a.occupied;
  remove = false;
  nneighbour = 0;
  nneighbour = 0;
  comment = "";
}

bool Atom::operator== (const Atom& a)
{
  if (label == a.label) return true;
  else return false;
}

Atom& Atom::operator= (const Atom& a)
{
  // std::cout << "operator= invoked" << std::endl;
  if (this == &a) return *this;
  name = a.name;
  label = a.label;
  r = a.r;
  occupied = a.occupied;
  nneighbour = 0;
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
      << std::right << std::setw(20) << a.r[2] \
      << std::right << std::setw(8) << a.occupied;
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

std::string Atom::write_cq(int spec_int, std::string cx, 
                     std::string cy, std::string cz)
{
  std::ostringstream oss;

  oss  << std::fixed << std::setprecision(8) \
       << std::right << std::setw(20) << r[0]*ang2bohr \
       << std::right << std::setw(20) << r[1]*ang2bohr \
       << std::right << std::setw(20) << r[2]*ang2bohr \
       << std::right << std::setw(4) << spec_int \
       << std::right << std::setw(2) << cx \
       << std::right << std::setw(2) << cy \
       << std::right << std::setw(2) << cz \
       << std::right << comment;

  return oss.str();
}

void Atom::move(Eigen::Vector3d trans)
{
  r += trans;
}

void Atom::occupy(bool occ)
{
  occupied = occ;
}

void Atom::add_nn(atom_ptr n)
{
  nn.push_back(n);
  nneighbour++;
}

void Atom::print_nn(void)
{
  std::cout << label << ": ";
  for (int i=0; i<nn.size(); i++){
    std::cout << (*nn[i]).label << " ";
  }
  std::cout << std::endl;
}