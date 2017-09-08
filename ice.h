#include "cell.h"
#include "water.h"
#include "hbond.h"
#include <gsl/gsl_rng.h>

class Ice: public Cell {
  private:
  public:
    std::deque<Water> waters;
    std::deque<Hbond> hbonds;
    int nwater, nhbond;

    Ice();
    virtual ~Ice();
    Ice(Cell&);

    void get_h_pos(void);
    void add_water(Water&);
    void add_hbond(Hbond&);
    void get_waters(void);
    void get_water_nn(double);
    void get_hbonds(void);
    void print_ice(void);
};