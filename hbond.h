#include "atom.h"
#include <set>

#ifndef HBOND_H
#define HBOND_H

class Water;

class Hbond {
  private:
  protected:
  public:
    int O1, O2, H1, H2, W1, W2;
    std::string bjerrum;
    std::set<int> w;

    Hbond();
    virtual ~Hbond();
    Hbond(const Hbond&);

    Hbond& operator= (const Hbond&);
    bool operator== (const Hbond&);
    Hbond(int, int, int, int, int, int);
};

#endif  // HBOND_H