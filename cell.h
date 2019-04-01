#include "atom.h"
#include <algorithm>
#include <regex>
#include <Eigen/LU>     // required for matrix.inverse()
#include <Eigen/Dense>  // required for dot and cross products

class Cell {
  private:
  protected:
    Eigen::MatrixXd dt;
    std::deque<Atom> ghosts;
    int nghosts;

    void wrap(void);
    void frac2cart(Atom&);
    void cart2frac(Atom&);
    Eigen::Vector3d mic_cart(Atom&, Atom&);
    Eigen::Vector3d mic_frac(Atom&, Atom&);
    Atom& get_atom(int);
    void get_ghosts(void);
    void get_dt(void);
    void get_nn(double);
    bool isPointOnLine(Atom&, Atom&, Atom&);
  public:
    std::string name;
    int natoms;
    Eigen::Matrix3d lat, lat_inv;
    std::deque<Atom> atoms;
    bool frac;
    Eigen::Vector3d scdim;

    Cell();
    virtual ~Cell();
    Cell(Eigen::Matrix3d);
    Cell(Eigen::Matrix3d, std::deque<Atom>);
    Cell(const Cell&);

    Cell& operator= (const Cell&);

    friend std::ostream& operator<< (std::ostream&, Cell&);

    void add_atom(Atom&);
    void read_cell(std::string);
    void cart2frac_all(void);
    void frac2cart_all(void);
    virtual void write_cell(std::string);
    virtual void write_cq(std::string);
    virtual void write_vasp(std::string);
    virtual void write_xsf(std::string, int, bool);
    void write_xyz(std::string, std::string);
    Cell super(int, int, int);
    void shift(double, double, double);
};
