#include "atom.h"

class Cell {
  private:
  public:
    std::string name;
    int natoms;
    double a, b, c, alpha, beta, gamma;
    Eigen::Matrix3f lat, lat_inv;
    std::vector<Atom> atoms;

    Cell();
    virtual ~Cell();
    Cell(double, double, double, double, double, double, std::vector<Atom>);
    Cell(Eigen::Matrix3f, std::vector<Atom>);
    Cell(const Cell&);
};
