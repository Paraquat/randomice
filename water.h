#include "atom.h"

#ifndef WATER_H
#define WATER_H

class Hbond;

class Water {
  private:
  protected:
  public:
    int O, H1, H2, H3, H4, target;
    std::string ionic;
    std::deque<int> hbonds;
    bool surface;

    Water();
    virtual ~Water();
    Water(const Water&);
    Water& operator= (const Water&);
    bool operator== (const Water&);
    Water(int, int, int, int, int);

    void add_hbond(int);
};

#endif  // WATER_H