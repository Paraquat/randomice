#include "atom.h"
#include "constants.h"
#include <algorithm>
#include <regex>
#include <Eigen/LU>     // required for matrix.inverse()
#include <Eigen/Dense>  // required for dot and cross products

class Cell {
  private:
  protected:
    Eigen::MatrixXd dt;
    Eigen::Matrix3d lat_inv;

    void wrap(void);
    void frac2cart(Atom&);
    void frac2cart_all(void);
    void cart2frac(Atom&);
    void cart2frac_all(void);
    Eigen::Vector3d mic_cart(Atom&, Atom&);
    Eigen::Vector3d mic_frac(Atom&, Atom&);
    Atom& get_atom(int);
    void get_dt(void);
    void get_nn(double);
    bool isPointOnLine(Atom&, Atom&, Atom&);
  public:
    std::string name;
    int natoms;
    Eigen::Matrix3d lat;
    std::deque<Atom> atoms;
    bool frac;

    Cell();
    virtual ~Cell();
    Cell(Eigen::Matrix3d);
    Cell(Eigen::Matrix3d, std::deque<Atom>);
    Cell(const Cell&);

    Cell& operator= (const Cell&);

    friend std::ostream& operator<< (std::ostream&, Cell&);

    void add_atom(Atom&);
    void read_cell(std::string);
    virtual void write_cell(std::string);
    void write_xyz(std::string, std::string);
    Cell super(int, int, int);
    void shift(double, double, double);
};