#include "cell.h"
#include <gsl/gsl_rng.h>

class Ice: public Cell {
  private:
  public:
    Ice();
    virtual ~Ice();
    Ice(Cell&);

    void get_h_pos(void);

};