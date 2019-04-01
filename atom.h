#include "constants.h"
#include "global.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <Eigen/Core>
#include <deque>
#include <memory>

#ifndef ATOM_H
#define ATOM_H

class Atom {
  private:
  public:
    std::string name;
    std::string comment;
    Eigen::Vector3d r;
    int label, nneighbour, species;
    bool occupied;  // whether the position is occupied
    bool remove;    // whether to remove the atom when writing to file
    bool fixed;     // whether this atom can be made occupied
    typedef std::shared_ptr<Atom> atom_ptr;
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
    std::string write_cq(int, std::string, std::string, std::string);
    std::string write_cell_highlight(std::string, std::string);
    std::string write_cell();
    bool operator== (const Atom&);

    void move(Eigen::Vector3d);
    void occupy(bool);
    void toggle_occupied(void);
    void add_nn(atom_ptr);
    void print_nn(void);
};

#endif  // ATOM_H