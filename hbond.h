#include "atom.h"

#ifndef HBOND_H
#define HBOND_H

class Water;
typedef std::shared_ptr<Water> water_ptr;

class Hbond {
  private:
  protected:
  public:
    Atom::atom_ptr O1, O2, H1, H2;
    std::string bjerrum;
    typedef std::shared_ptr<Hbond> hbond_ptr;

    Hbond();
    virtual ~Hbond();
    Hbond(const Hbond&);

    friend std::ostream& operator<< (std::ostream&, Hbond&);
    Hbond& operator= (const Hbond&);
    Hbond(Atom::atom_ptr, Atom::atom_ptr, Atom::atom_ptr, Atom::atom_ptr);

    void check_defect(void);
};

#endif  // HBOND_H