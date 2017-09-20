#include "atom.h"

#ifndef WATER_H
#define WATER_H

class Hbond;
typedef std::shared_ptr<Hbond> hbond_ptr;

class Water {
  private:
    int nH;
  protected:
  public:
    Atom::atom_ptr O, H1, H2, H3, H4;
    std::string ionic;
    int label;
    typedef std::shared_ptr<Water> water_ptr;
    std::deque<water_ptr> nn;
    std::deque<hbond_ptr> hbonds;

    Water();
    virtual ~Water();
    Water(const Water&);
    Water& operator= (const Water&);
    bool operator== (const Water&);
    Water(int, Atom::atom_ptr, Atom::atom_ptr, Atom::atom_ptr, Atom::atom_ptr, Atom::atom_ptr);

    friend std::ostream& operator<< (std::ostream&, Water&);
    friend std::ofstream& operator<< (std::ofstream&, Water&);
    void add_nn(water_ptr);
    void add_hbond(hbond_ptr);
    void check_defect(void);
    int coord(void);
};

#endif  // WATER_H