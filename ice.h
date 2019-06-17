#include "cell.h"
#include "water.h"
#include "hbond.h"
#include <vector>
#include <boost/progress.hpp>
#include <boost/format.hpp>

struct Node {
    int water;
    int hbond;

    Node(){};
    virtual ~Node(){}; 
    Node(int w, int h){
        water = w;
        hbond = h;
    }
    Node(const Node& n){
        water = n.water;
        hbond = n.hbond;
    }
};

class Ice: public Cell {
  private:
    gsl_rng *rp;
    bool ghost_method;
    double oh_bond_length;
  public:
    std::deque<Water> waters;
    std::deque<Hbond> hbonds;
    std::deque<int> s1list, s2list;
    int nwater, nhbond, nbilayer, noccupied;
    int xlayers, ylayers, zlayers;

    Ice();
    virtual ~Ice();
    Ice(const Ice&);
    Ice(Cell&);

    void set_oh_length(double);
    void read_h_pos(Cell&);
    Ice super(int, int, int);
    void super_self(int, int, int);
    void get_h_pos(void);
    void add_water(Water&);
    void add_hbond(Hbond&);
    void get_waters(void);
    void get_hbonds(void);
    void init_bulk_random(double);
    void init_bulk_ordered(double, Cell, int, int, int);
    void print_ice(void);

    void init_rng(void);
    int rng_int(int);
    double rng_uniform(void);
    void populate_h_random(void);
    void fix_water(int);
    int water_coord(int);
    int hbond_occ(int);
    int hbond_single_occupied(void); 
    int o_two_coordinated(void);
    bool swap_h(int);
    bool rotate_water(int);
    void buch_mc_correct(void);
    void add_defects(int, int);
    void protonate_water(int);
    void deprotonate_water(int);
    int get_random_h(int);
    int hb_target(int);

    std::vector<bool> save_config(void);
    void revert_config(std::vector<bool>);
    std::deque<Node> get_loop(int);
    void rick_move(int);
    void rick_randomise(int);

    int check_ionic_defects(void);
    int check_bjerrum_defects(void);
    void check_defects(void);
    Eigen::Vector3d water_dipole(int);
    double c1_dipole(void);

    void build_slab(double, int);
    std::deque<int> find_dOH(int);
    double order_parameter(double);
    void build_ordered_slab(double, int, double, int);
    void build_step(std::string, int, double, double, bool, std::string);

    void write_xsf(std::string, int, bool);
    void write_cell(std::string);
    void write_chunk_cell(std::string, std::deque<int>);
    void write_highlight_cell(std::string, std::deque<int>);
    void print_water(int);
    void print_hbond(int);
    void print_network(void);
};
