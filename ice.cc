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

// Check how many hydrogens there are on a hydrogen bond
int Ice::hbond_occ(int hb)
{
  int nH = 0;
  if (get_atom(hbonds[hb].H1).occupied) nH++;
  if (get_atom(hbonds[hb].H2).occupied) nH++;
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

// Revert the configuration to saved
void Ice::revert_config(std::vector<bool> conf)
{
  for (int i=0; i<natoms; i++) atoms[i].occupied = conf[i];
}

// Find a closed loop of hydrogen bonds (in direction of donation). May cross
// cell boundary and end in a image
std::deque<Node> Ice::get_loop(void)
{
  std::deque<Node> loop;
  int current, hb, target;
  bool end;

  init_rng();
  current = rng_int(nwater);
  while (true){
    hb = rng_int(4);
    target = hb_target(waters[current].hbonds[hb]);
    if (target != current) break;
  }
  loop.push_back(Node(current, waters[current].hbonds[hb]));

  while (true){
    end = false;
    while (true){
      hb = rng_int(4);
      target = hb_target(waters[current].hbonds[hb]);
      if (target != current) break;
    }
    for (int j=0; j<loop.size(); j++){
      if (target == loop[j].water){
        end = true;
        break;
      }
    }
    // std::cout << current << ": " << hbonds[waters[current].hbonds[hb]].W1 << " " << hbonds[waters[current].hbonds[hb]].W2 << " -> " << target << std::endl;
    loop.push_back(Node(current, waters[current].hbonds[hb]));
    current = target;
    if (end) break;
  }
  // Trim the vector to get just the loop part
  while (true){
    if (loop[0].water == target) break;
    else loop.pop_front();
  }
  return loop;
}

// Perform a move of the Rick algorithm (Rick/Haymet JCP  118, 9291 (2003))
void Ice::rick_move(void)
{
  std::deque<Node> loop;
  loop = get_loop();

  // write_highlight_cell("test1.cell", loop2);
  for (int i=0; i<loop.size(); i++){
    swap_h(loop[i].hbond);
  }
  // write_highlight_cell("test2.cell", loop2);
}

// Randomise water orientations using Rick algorithm with C1 dipole contraint
void Ice::rick_randomise(int max_loops)
{
  std::cout << "Randomising water orientations via Rick algorithm..." 
            << std::endl;
  int ndefect;
  double cell_dipole, cell_dipole_old, rn, mcp;
  std::vector<bool> conf;

  init_rng();

  cell_dipole_old = c1_dipole();
  conf = save_config();
  for (int i=0; i<max_loops; i++){
    while (true){
      while (true){
        rick_move(); 
        ndefect = check_ionic_defects();
        if (ndefect == 0) break;
        else revert_config(conf);
      }
      // Monte Carlo criterion
      cell_dipole = c1_dipole();
      if (cell_dipole > cell_dipole_old){
        mcp = exp(cell_dipole_old - cell_dipole);
        rn = rng_uniform();
        if (rn > mcp) revert_config(conf);
      }
      else break;
    }
    conf = save_config();
    std::cout << "Cell dipole = " << cell_dipole_old << std::endl;
    cell_dipole_old = cell_dipole;
  }
}

// Count ionic defects
int Ice::check_ionic_defects(void)
{
  int coord;
  int nionic = 0;
  for (int i=0; i<nwater; i++){
    coord = water_coord(i);
    switch(coord) {
      case 0:
        waters[i].ionic = "O2-";
        nionic++;
        break;
      case 1:
        waters[i].ionic = "OH-";
        nionic++;
        break;
      case 2:
        waters[i].ionic = "None";
        break;
      case 3:
        waters[i].ionic = "H3O+";
        nionic++;
        break;
      case 4:
        waters[i].ionic = "H4O2+";
        nionic++;
        break;
    }
  }
  return nionic;
}

// Count Bjerrum defects
int Ice::check_bjerrum_defects(void)
{
  int occ;
  int nbjerrum = 0;
  for (int i=0; i<nhbond; i++){
    occ = hbond_occ(i);
    switch(occ) {
      case 0:
        hbonds[i].bjerrum = "L";
        nbjerrum++;
        break;
      case 1:
        hbonds[i].bjerrum = "None";
        break;
      case 2:
        hbonds[i].bjerrum = "D";
        nbjerrum++;
        break;
    }
  }
  return nbjerrum;
}

// Get the dipole vector of water w (multiply orientation vector by constant
// 3.09 -- Batista et al JCP 109, 4546 (1998))
Eigen::Vector3d Ice::water_dipole(int w)
{
  Eigen::Vector3d oh1, oh2;

  int nh = 0;
  int io = waters[w].O;
  int ih1, ih2;
  assert(water_coord(w) == 2);
  if (get_atom(waters[w].H1).occupied) {
    if (nh == 0) ih1 = waters[w].H1;
    if (nh == 1) ih2 = waters[w].H1;
    nh++;
  }
  if (get_atom(waters[w].H1).occupied) {
    if (nh == 0) ih1 = waters[w].H2;
    if (nh == 1) ih2 = waters[w].H2;
    nh++;
  }
  if (get_atom(waters[w].H1).occupied) {
    if (nh == 0) ih1 = waters[w].H3;
    if (nh == 1) ih2 = waters[w].H3;
    nh++;
  }
  if (get_atom(waters[w].H1).occupied) {
    if (nh == 0) ih1 = waters[w].H4;
    if (nh == 1) ih2 = waters[w].H4;
    nh++;
  }

  oh1 = mic_cart(get_atom(io), get_atom(ih1));
  oh2 = mic_cart(get_atom(io), get_atom(ih2));
  return (oh1 + oh2).normalized()*ice_h2o_dipole_mag;
}

// Check the dipole moment of the unit cell using the C1 constraint
// (Hayward/Reimers JCP 106, 1518 (1997))
double Ice::c1_dipole(void)
{
  Eigen::Vector3d cell_dipole;
  cell_dipole << 0.0, 0.0, 0.0;

  for (int i=0; i<nwater; i++){
    cell_dipole += water_dipole(i);
  }
  return cell_dipole.norm();
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

// Write object to a .cell file
void Ice::write_chunk_cell(std::string fname, std::deque<int> wlist)
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

  for (int i=0; i<wlist.size(); i++) {
    outfile << get_atom(waters[wlist[i]].O) << std::endl;
    if (atoms[waters[wlist[i]].H1].occupied){
      outfile << get_atom(waters[wlist[i]].H1) << " " << wlist[i] << std::endl;
    }
    if (atoms[waters[wlist[i]].H2].occupied){
      outfile << get_atom(waters[wlist[i]].H2) << " " << wlist[i] << std::endl;
    }
    if (atoms[waters[wlist[i]].H3].occupied){
      outfile << get_atom(waters[wlist[i]].H3) << " " << wlist[i] << std::endl;
    }
    if (atoms[waters[wlist[i]].H4].occupied){
      outfile << get_atom(waters[wlist[i]].H4) << " " << wlist[i] << std::endl;
    }
  }

  if (frac) outfile << "%ENDBLOCK positions_frac" << std::endl << std::endl;
  else outfile << "%ENDBLOCK positions_abs" << std::endl << std::endl;
}

// Write object to a .cell file
void Ice::write_highlight_cell(std::string fname, std::deque<int> wlist)
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
    Atom o = get_atom(waters[i].O);
    bool replace = false;
    for (int j=0; j<wlist.size(); j++){
      if (wlist[j] == waters[i].O) {
        replace = true;
      }
    }
    if (replace){
      Atom temp("N", o.label, o.r, o.occupied);
      outfile << temp << std::endl;
    } else{
      outfile << get_atom(waters[i].O) << std::endl;
    }
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