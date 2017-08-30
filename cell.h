#include "atom.h"
#include <algorithm>
#include <regex>
#include <iostream>
#include <Eigen/LU>     // required for matrix.inverse()

class Cell {
  private:
  public:
    std::string name;
    int natoms;
    bool frac;
    Eigen::Matrix3d lat, lat_inv;
    std::vector<Atom> atoms;
    Eigen::MatrixXd dt;

    Cell();
    virtual ~Cell();
    Cell(Eigen::Matrix3d, std::vector<Atom>);
    Cell(const Cell&);

    friend std::ostream& operator<< (std::ostream&, Cell&);
    // friend std::ofstream& operator<< (std::ofstream&, Cell&);

    void add_atom(Atom);
    void read_cell(std::string);
    void write_cell(std::string);
    // void recell(void);
    // void uncell(void);

    // double frac2cart_disp(Eigen::Vector3d);
    void wrap(void);
    void frac2cart(Atom&);
    void frac2cart_all(void);
    void cart2frac(Atom&);
    void cart2frac_all(void);
    Eigen::Vector3d mic_cart(Atom&, Atom&);
    Eigen::Vector3d mic_frac(Atom&, Atom&);
    void get_dt(void);
};
