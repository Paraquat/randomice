#include "atom.h"

#ifndef WATER_H
#define WATER_H

class Hbond;

class Water {
  private:
  protected:
  public:
    int O, H1, H2, H3, H4, bilayer;
    std::string ionic;
    std::deque<int> hbonds;
    bool surface1, surface2; // is the molecule on either surface?
    bool remove;             // remove this molecule to make a step?
    bool dOH;                // does the molecule have a dangling H?
    int step;                // label for step

    Water();
    virtual ~Water();
    Water(const Water&);
    Water& operator= (const Water&);
    bool operator== (const Water&);
    Water(int, int, int, int, int);

    void add_hbond(int);
};

#endif  // WATER_H