#include "atom.h"
#include <algorithm>
#include <regex>
#include <iostream>

class Cell {
  private:
  public:
    std::string name;
    int natoms;
    Eigen::Matrix3f lat = Eigen::Matrix3f::Zero();
    Eigen::Matrix3f lat_inv = Eigen::Matrix3f::Zero();
    std::vector<Atom> atoms;

    Cell();
    virtual ~Cell();
    Cell(Eigen::Matrix3f, std::vector<Atom>);
    Cell(const Cell&);

    friend std::ostream& operator<< (std::ostream&, Cell&);
    // friend std::ofstream& operator<< (std::ofstream&, Cell&);

    void add_atom(Atom);
    void read_cell(std::string);
    // void write_cell(std::string);
    // void recell(void);
    // void uncell(void);

    // double frac2cart_disp(Eigen::Vector3d);
    // Eigen::Vector3d frac2cart(Eigen::Vector3d);
    // Eigen::Vector3d cart2frac(Eigen::Vector3d);
    // void frac2cart_cell(void);
    // void cart2frac_cell(void);
};
