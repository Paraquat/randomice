#include "cell.h"
#include "water.h"
#include "hbond.h"
#include <gsl/gsl_rng.h>

class Ice: public Cell {
  private:
    gsl_rng *r;
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
    void get_hbonds(void);
    void print_ice(void);

    void init_rng(void);
    int rng_int(int);
    double rng_uniform(void);
    void populate_h_random(void);
    int water_coord(int);
    int o_two_coordinated(void);
    void swap_h(int);
    void buch_mc_correct(void);

    void write_cell(std::string);
    void print_water(int);
    void print_hbond(int);
};