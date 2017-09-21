#include "atom.h"

#ifndef WATER_H
#define WATER_H

class Hbond;

class Water {
  private:
  protected:
  public:
    int O, H1, H2, H3, H4;
    std::string ionic;
    std::deque<int> nn;
    std::deque<int> hbonds;

    Water();
    virtual ~Water();
    Water(const Water&);
    Water& operator= (const Water&);
    bool operator== (const Water&);
    Water(int, int, int, int, int);

    void add_nn(int);
    void add_hbond(int);
};

#endif  // WATER_H