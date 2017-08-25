#include "cell.h"

Cell::Cell()
{
  natoms = 0
}

virtual Cell::~Cell()
{
}

Cell::Cell(double ai, double bi, double ci, double alp, double bet, \
           double gam, std::vector<Atom> r)
{
  a = ai;
  b = bi;
  c = ci;
  alpha = alp;
  beta = bet;
  gamma = gam;
  atoms = r;
  natoms = r.size();
}

Cell::Cell(Eigen::Matrix3f l, std::vector<Atom> r)
{
  lat = l;
  atoms = r;
}

Cell::Cell(const Cell& cell)
{
  a = cell.a;
  b = cell.b;
  c = cell.c;
  alpha = cell.alpha;
  beta = cell.beta;
  gamma = cell.gamma;
  lat = cell.lat;
  atoms = cell.atoms;
  natoms = cell.natoms;
}
