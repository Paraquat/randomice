#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <Eigen/Core>
#include <deque>
#include <boost/shared_ptr.hpp>

#ifndef ATOM_H
#define ATOM_H

class Atom {
  private:
  public:
    std::string name;
    Eigen::Vector3d r;
    int label, nneighbour;
    bool occupied;
    typedef boost::shared_ptr<Atom> atom_ptr;
    std::deque<atom_ptr> nn;

    Atom();
    virtual ~Atom();
    Atom(std::string, int, double, double, double);
    Atom(std::string, int, Eigen::Vector3d);
    Atom(std::string, int, Eigen::Vector3d, bool);
    Atom(const Atom&);
    Atom& operator= (const Atom&);
    friend std::ostream& operator<< (std::ostream&, Atom&);
    friend std::ofstream& operator<< (std::ofstream&, Atom&);
    bool operator== (const Atom&);

    void add_nn(atom_ptr);
    void print_nn(void);
};

#endif  // ATOM_H