#include "ice.h"

// Default constructor
Ice::Ice()
{
}

// Destructor
Ice::~Ice()
{
}

Ice::Ice(Cell& cell)
{
  lat = cell.lat;
  atoms = cell.atoms;
  natoms = cell.natoms;
  frac = cell.frac; 
}

// Compute all hydrogen positions for a given oxygen lattice
void Ice::get_h_pos(void)
{
  std::cout << "Computing hydrogen positions from oxygens" << std::endl;
  Eigen::Vector3d oo, r_h;
  std::string species = "H";
  bool occ = false;
  std::deque<Atom> hlist;

  if (frac == true) frac2cart_all();

  get_dt();
  for (int i=0; i<natoms; i++){
    for (int j=0; j<natoms; j++){
      if (i != j){
      if (dt(i,j) < oo_max){
        oo = mic_cart(atoms[i], atoms[j]);
        r_h = atoms[i].r + oh_default*oo.normalized();
        Atom h = Atom(species, natoms+1, r_h, occ);
        hlist.push_back(h);
      }
      }
    }
  }

  for (int i=0; i<hlist.size(); i++){
    add_atom(hlist[i]);
  }

}

void Ice::add_water(Water& w)
{
  waters.push_back(w);
  nwater++;
}

void Ice::add_hbond(Hbond& h)
{
  hbonds.push_back(h);
  nhbond++;
}

void Ice::get_waters(void)
{
  get_dt();
  for (int i=0; i<natoms; i++){
    if (atoms[i].name == "O"){
      int io = i;
      int ih1, ih2, ih3, ih4;
      int hcount = 0;
      for (int j=0; j<natoms; j++){
        if (dt(i,j) < oh_max){
          if (atoms[j].name == "H"){
            hcount++;
            if (hcount == 1) ih1 = j;
            else if (hcount == 2) ih2 = j;
            else if (hcount == 3) ih3 = j;
            else if (hcount == 4) ih4 = j;
          }
        } 
      }
      if (hcount != 4){
        throw std::runtime_error("Did not find 4 hydrogens for this water");
      }
      Water w(io, ih1, ih2, ih3, ih4);
      add_water(w);
    }
  }
}

void Ice::get_hbonds(void)
{
  std::cout << "Finding hydrogen bonds" << std::endl;
  bool done;

  nhbond = 0;
  get_dt();

  for (int i=0; i<nwater; i++){
    int io1 = waters[i].O;
    for (int j=i+1; j<nwater; j++){
      int io2 = waters[j].O;
      if (dt(io1,io2) < oo_max){
        int ih1, ih2;

        std::set<int> ws;
        ws.insert(i);
        ws.insert(j);

        done = false;
        for (int k=0; k<hbonds.size(); k++){
          if (ws == hbonds[k].w){
            done = true;
            break;
          }
        }
        if (done) continue;

        int nh = 0;
        if (isPointOnLine(atoms[io1], atoms[io2], atoms[waters[i].H1])){
          nh++; 
          ih1 = waters[i].H1;
        }
        if (isPointOnLine(atoms[io1], atoms[io2], atoms[waters[i].H2])){
          nh++; 
          ih1 = waters[i].H2;
        }
        if (isPointOnLine(atoms[io1], atoms[io2], atoms[waters[i].H3])){
          nh++; 
          ih1 = waters[i].H3;
        }
        if (isPointOnLine(atoms[io1], atoms[io2], atoms[waters[i].H4])){
          nh++; 
          ih1 = waters[i].H4;
        }
        if (isPointOnLine(atoms[io1], atoms[io2], atoms[waters[j].H1])){
          nh++; 
          ih2 = waters[j].H1;
        }
        if (isPointOnLine(atoms[io1], atoms[io2], atoms[waters[j].H2])){
          nh++; 
          ih2 = waters[j].H2;
        }
        if (isPointOnLine(atoms[io1], atoms[io2], atoms[waters[j].H3])){
          nh++; 
          ih2 = waters[j].H3;
        }
        if (isPointOnLine(atoms[io1], atoms[io2], atoms[waters[j].H4])){
          nh++; 
          ih2 = waters[j].H4;
        }
        assert(nh == 2);
        Hbond hb(io1, io2, ih1, ih2, i, j);
        add_hbond(hb);

      }
    }
  }

  for (int i=0; i<nhbond; i++){
    int io1 = hbonds[i].O1;
    int io2 = hbonds[i].O2;
    for (int j=0; j<nwater; j++){
      if (waters[j].O == io1) waters[j].add_hbond(i);
      else if (waters[j].O == io2) waters[j].add_hbond(i);
    }
  }
  for (int i=0; i<nwater; i++) assert(waters[i].hbonds.size() == 4);
}

void Ice::print_ice(void)
{
  std::cout << nwater << " water molecules" << std::endl << std::endl;
  for (int i=0; i<nwater; i++){
    std::cout << get_atom(waters[i].O) << std::endl;
    std::cout << get_atom(waters[i].H1) << std::endl;
    std::cout << get_atom(waters[i].H2) << std::endl;
    std::cout << get_atom(waters[i].H3) << std::endl;
    std::cout << get_atom(waters[i].H4) << std::endl;
    std::cout << std::endl;
  }
}

void Ice::init_rng(void)
{
  r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_env_setup();
  gsl_rng_set(r, time(0));
}

int Ice::rng_int(int a)
{
  return gsl_rng_get(r)%a;
}

double Ice::rng_uniform(void)
{
  return gsl_rng_uniform(r);
}

// Populate each hydrogen bond with one hydrogen, pick one of two randoms sites
void Ice::populate_h_random(void)
{
  std::cout << "Randomly assigning hydrogens, one per hydrogen bond" << std::endl;
  int rn;
  init_rng();
  for (int i=0; i<nhbond; i++){
    rn = rng_int(2);
    if (rn == 0) {
      atoms[hbonds[i].H1].occupy(true);
      atoms[hbonds[i].H2].occupy(false);
    } else{
      atoms[hbonds[i].H1].occupy(false);
      atoms[hbonds[i].H2].occupy(true);
    }
  }
}

// Check the H coordination number of an oxygen
int Ice::water_coord(int w)
{
  int nH = 0;
  if (get_atom(waters[w].H1).occupied) nH++;
  if (get_atom(waters[w].H2).occupied) nH++;
  if (get_atom(waters[w].H3).occupied) nH++;
  if (get_atom(waters[w].H4).occupied) nH++;
  return nH;
}

// count the number of two-coordinated oxygens
int Ice::o_two_coordinated(void)
{
  int ntwocoord = 0;
  for (int i=0; i<nwater; i++){
    if (water_coord(i) == 2){
      ntwocoord++;
    }
  }
  return ntwocoord;
}

// Swap the H position on a hydrogen bond
void Ice::swap_h(int hb)
{
  atoms[hbonds[hb].H1].occupied = !atoms[hbonds[hb].H1].occupied;
  atoms[hbonds[hb].H2].occupied = !atoms[hbonds[hb].H2].occupied;
}

// Monte Carlo algorithm to correct ice rules (Buch et al JCP B 102:8641, 1998)
void Ice::buch_mc_correct(void)
{
  std::cout << "Enforcing ice rules via Buch algorithm..." << std::endl;
  int c1a, c2a, cdiffa, c1b, c2b, cdiffb;

  init_rng();
  while (o_two_coordinated() != nwater){
    // Pick a random h-bond
    int hb_ind = rng_int(nhbond);
    // before the h swap
    c1a = water_coord(hbonds[hb_ind].W1);
    c2a = water_coord(hbonds[hb_ind].W2);
    cdiffa = std::abs(c1a - c2a);
    // after the h swap
    swap_h(hb_ind);
    c1b = water_coord(hbonds[hb_ind].W1);
    c2b = water_coord(hbonds[hb_ind].W2);
    cdiffb = std::abs(c1b - c2b);
    //No change in the difference
    if (cdiffb - cdiffa == 0){
      // move with probability 1/2
      double rn = rng_uniform();
      if (rn < 0.5) swap_h(hb_ind);
    // Decrease in the coordination difference after move
    } else if (cdiffb - cdiffa < 0){
      // accept move

    // Increase in the coordination difference after move
    } else if (cdiffb - cdiffa > 0){
      // reject move
      swap_h(hb_ind);
    }
  }
}

// Return the direction of a hydrogen bond. If it is contains a Bjerrum
// defect, return -1.
int Ice::hb_target(int hb)
{
  int target;
  if (get_atom(hbonds[hb].H1).occupied && get_atom(hbonds[hb].H2).occupied){
    target = -1;
  } 
  if (get_atom(hbonds[hb].H1).occupied) target = hbonds[hb].W2;
  else if (get_atom(hbonds[hb].H2).occupied) target = hbonds[hb].W1;
  else target = -1;
  if (target == -1){
    throw std::runtime_error("Bjerrum defect");
  } else return target;
}

// Save the configuration in case of reset after Monte Carlo move
std::vector<bool> Ice::save_config(void)
{
  std::vector<bool> conf;
  for (int i=0; i<natoms; i++) conf.push_back(atoms[i].occupied);
  return conf;
}

// Find a closed loop of hydrogen bonds (in direction of donation). May cross
// cell boundary and end in a image
std::vector<Ice::i2ple> Ice::get_loop(void)
{
  int max_steps = 1000;
  std::vector<i2ple> loop;
  int w, hb, target;
  bool end;

  // Pick a random starting point and a random donor hbond
  init_rng();
  w = rng_int(nwater);
  while (true){
    hb = rng_int(4);
    target = hb_target(waters[w].hbonds[hb]);
    if (target != w) break;
  }
  i2ple step(w, hb);
  loop.push_back(step);

  for (int i=0; i<max_steps; i++){
    end = false;
    w = hb_target(std::get<1>(loop.back()));
    // Pick a hbond at random
    hb = rng_int(4);
    target = hb_target(waters[w].hbonds[hb]);
    if (target == w) continue; // if the bond is a donor to this water
    else {
      for (int j=0; j<loop.size(); j++){
        if (target == std::get<0>(loop[j])){
          end = true;
          break;
        }
      }
    }

    if (end) break;
    else {
      i2ple step(target, hb);
      loop.push_back(step);
    std::cout << w << ": " << hbonds[waters[w].hbonds[hb]].W1 << " " << hbonds[waters[w].hbonds[hb]].W2 << " -> " << target << std::endl;
    std::cout << std::get<0>(step) << " " << std::get<1>(step) << std::endl;
    }
  }
  return loop;
}

// Randomise ice lattice using rick algorithm
void Ice::rick_algo(void)
{
  std::vector<i2ple> loop;

  loop = get_loop();
  // for (int i=0; i<loop.size(); i++){
  //   std::cout << std::get<0>(loop[i]) << " " << std::get<1>(loop[i]) << std::endl;
  // }
  // std::cout << std::endl;
}

// Write object to a .cell file
void Ice::write_cell(std::string fname)
{
  std::ofstream outfile(fname.c_str());

  if (outfile.fail()) {
    std::cout << "Error opening file " << fname << std::endl;
  }

  outfile << "%BLOCK lattice_cart" << std::endl;
  for (int i=0; i<=2; i++) {
    outfile << std::fixed << std::setprecision(8) \
            << lat.row(i) << std::endl;
  }
  outfile << "%ENDBLOCK lattice_cart" << std::endl << std::endl;
  if (frac) outfile << "%BLOCK positions_frac" << std::endl;
  else outfile << "%BLOCK positions_abs" << std::endl;

  for (int i=0; i<nwater; i++) {
    outfile << get_atom(waters[i].O) << std::endl;
    if (atoms[waters[i].H1].occupied){
      outfile << get_atom(waters[i].H1) << std::endl;
    }
    if (atoms[waters[i].H2].occupied){
      outfile << get_atom(waters[i].H2) << std::endl;
    }
    if (atoms[waters[i].H3].occupied){
      outfile << get_atom(waters[i].H3) << std::endl;
    }
    if (atoms[waters[i].H4].occupied){
      outfile << get_atom(waters[i].H4) << std::endl;
    }
  }

  if (frac) outfile << "%ENDBLOCK positions_frac" << std::endl << std::endl;
  else outfile << "%ENDBLOCK positions_abs" << std::endl << std::endl;
}

void Ice::print_water(int w)
{
  std::cout << get_atom(waters[w].O) << std::endl
            << get_atom(waters[w].H1) << std::endl
            << get_atom(waters[w].H2) << std::endl
            << get_atom(waters[w].H3) << std::endl
            << get_atom(waters[w].H4) << std::endl << std::endl;
}

void Ice::print_hbond(int hb)
{
  std::cout << get_atom(hbonds[hb].O1) << std::endl
            << get_atom(hbonds[hb].O2) << std::endl
            << get_atom(hbonds[hb].H1) << std::endl
            << get_atom(hbonds[hb].H2) << std::endl << std::endl;
}

void Ice::print_network(void)
{
  for (int i=0; i<nwater; i++){
    std::cout << waters[i].O << ": ";
    for (int j=0; j<waters[i].hbonds.size(); j++){
      std::cout << std::setw(6) << hb_target(waters[i].hbonds[j]) << " ";
    }
    std::cout << std::endl;
  }
}