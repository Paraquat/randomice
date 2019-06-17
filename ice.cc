#include "ice.h"

// Default constructor
Ice::Ice()
{
}

// Destructor
Ice::~Ice()
{
}

Ice::Ice(const Ice& ice)
{
  lat = ice.lat;
  lat_inv = lat.inverse();
  atoms = ice.atoms;
  natoms = ice.natoms;
  frac = ice.frac;
  ghost_method = false;
  nwater = 0;
}

Ice::Ice(Cell& cell)
{
  lat = cell.lat;
  lat_inv = lat.inverse();
  atoms = cell.atoms;
  natoms = cell.natoms;
  frac = cell.frac; 
  scdim = cell.scdim;
  ghost_method = false;
  nwater = 0;
}

// set the OH bond length
void Ice::set_oh_length(double oh)
{
  oh_bond_length = oh;
}

// Given a cell object containing hydrogens, construct and ice object
void Ice::read_h_pos(Cell& cell)
{
  ghost_method = false;
  lat = cell.lat;
  lat_inv = lat.inverse();
  frac = cell.frac; 
  assert(frac == true);

  // Add the oxygen atoms
  for (int i=0; i<cell.natoms; i++){
    if (cell.atoms[i].name == "O"){
      add_atom(cell.atoms[i]);
    }
  }
  if (natoms < 30) ghost_method = true;
  get_h_pos();
  get_waters();
  cart2frac_all();
  wrap();

  Eigen::Vector3d s;
  double d;
  double thresh = 0.2;
  int nH = 0;
  for (int i=0; i<cell.natoms; i++){
    if (cell.atoms[i].name == "H"){
      for (int j=0; j<natoms; j++){
        if (atoms[j].name == "H"){
          s = mic_frac(cell.atoms[i], atoms[j]);
          d = s.norm();
          if (d <= thresh){
            atoms[j].occupied = true;
            nH++;
          }
        }
      }
    }
  }
}

// Generate a axbxc supercell from an *ordered* unit cell
Ice Ice::super(int a, int b, int c)
{
  Eigen::Vector3d t, rt;
  Eigen::Matrix3d lat_super;

  assert(frac == true);
  lat_super.row(0) = lat.row(0)*static_cast<double>(a);
  lat_super.row(1) = lat.row(1)*static_cast<double>(b);
  lat_super.row(2) = lat.row(2)*static_cast<double>(c);
  Cell sc(lat_super);
  sc.lat_inv = sc.lat.inverse();

  scdim << static_cast<double>(a), static_cast<double>(b), \
           static_cast<double>(c);
  int label = 1;
  for (int i=0; i<a; i++){
    for (int j=0; j<b; j++){
      for (int k=0; k<c; k++){
        t << static_cast<double>(i), static_cast<double>(j), \
             static_cast<double>(k);

        for (int n=0; n<natoms; n++){
          rt = atoms[n].r + t;
          for (int m=0; m<=2; m++){
            rt(m) = rt(m)/scdim(m);
          }
          Atom a(atoms[n].name, label, rt);
          a.occupied = atoms[n].occupied;
          sc.add_atom(a);
          label += 1;
        }
      }
    }
  }
  Ice ice_sc(sc);
  ice_sc.get_waters();
  std::string fname = "test.cell";
  ice_sc.write_cell(fname);
  return ice_sc; 
}

// Convert self to a supercell
void Ice::super_self(int a, int b, int c)
{
  Eigen::Vector3d t, rt;
  Eigen::Matrix3d lat_super;

  std::cout << "Making supercell" << std::endl;

  Ice ice_unit(*this);
  ice_unit.get_waters();
  // ice_unit.print_ice();

  assert(frac == true);
  lat_super.row(0) = lat.row(0)*static_cast<double>(a);
  lat_super.row(1) = lat.row(1)*static_cast<double>(b);
  lat_super.row(2) = lat.row(2)*static_cast<double>(c);
  lat = lat_super;
  lat_inv = lat.inverse();


  int natoms_old = natoms;
  atoms.clear();
  natoms = 0;
  waters.clear();
  nwater = 0;

  scdim << static_cast<double>(a), static_cast<double>(b), \
           static_cast<double>(c);
  int label = 1;
  for (int i=0; i<a; i++){
    for (int j=0; j<b; j++){
      for (int k=0; k<c; k++){
        t << static_cast<double>(i), static_cast<double>(j), \
             static_cast<double>(k);

        for (int n=0; n<natoms_old; n++){
          rt = ice_unit.atoms[n].r + t;
          for (int m=0; m<=2; m++){
            rt(m) = rt(m)/scdim(m);
          }
          Atom a(ice_unit.atoms[n].name, label, rt, 
                 ice_unit.atoms[n].occupied);
          if (a.occupied) noccupied++;
          add_atom(a);
        }
      }
    }
  }
  ghost_method = false;
  get_waters();
}

// Compute all hydrogen positions for a given oxygen lattice
void Ice::get_h_pos(void)
{
  std::cout << "Computing hydrogen positions from oxygens" << std::endl;
  Eigen::Vector3d oo, r_h;
  std::string species = "H";
  bool occ = false;
  std::deque<Atom> hlist;
  int nH = 0;

  if (frac == true) frac2cart_all();
  // if (frac == false) cart2frac_all();

  if (ghost_method){
    get_ghosts(); 
    for (int i=0; i<natoms; i++){
      for (int j=0; j<nghosts; j++){
        if (atoms[i].label != ghosts[j].label){
          oo = ghosts[j].r - atoms[i].r;
          if (oo.norm() < oo_max){
            r_h = atoms[i].r + oh_bond_length*oo.normalized();
            Atom h = Atom(species, natoms+1, r_h, occ);
            hlist.push_back(h);
            nH++;
          }
        }
      }
    }
  } else {
    get_dt();
    for (int i=0; i<natoms; i++){
      for (int j=0; j<natoms; j++){
        if (i != j){
          if (dt(i,j) < oo_max){
            oo = mic_cart(atoms[i], atoms[j]);
            r_h = atoms[i].r + oh_bond_length*oo.normalized();
            Atom h = Atom(species, natoms+1, r_h, occ);
            hlist.push_back(h);
            nH++;
          }
        }
      }
    }
  }

  std::cout << "Found " << nH << " hydrogen positions" << std::endl;
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
  std::cout << "Finding water molecules" << std::endl;

  if (!ghost_method) get_dt();
  double d;
  for (int i=0; i<natoms; i++){
    if (atoms[i].name == "O"){
      int io = i;
      int ih1, ih2, ih3, ih4;
      int hcount = 0;
      for (int j=0; j<natoms; j++){
        if (ghost_method){
          d = (atoms[i].r - atoms[j].r).norm();
        } else{
          d = dt(i,j);
        }
        if (d < oh_max){
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
        std::cout << hcount << " hydrogens" << std::endl;
        throw std::runtime_error("Did not find 4 hydrogens for this water");
      }
      Water w(io, ih1, ih2, ih3, ih4);
      add_water(w);
    }
  }
  std::cout << "Found " << nwater << " water molecules" <<std::endl;
}

void Ice::get_hbonds(void)
{
  std::cout << "Finding hydrogen bonds" << std::endl;
  bool done;

  nhbond = 0;

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
  if (flag_debug){
    std::cout << "Found " << nhbond << " hydrogen bonds" << std::endl;
  }
}

void Ice::init_bulk_random(double oh)
{
  set_oh_length(oh);
  get_h_pos();
  get_waters();
  get_hbonds();
  populate_h_random();
  buch_mc_correct();
}

void Ice::init_bulk_ordered(double oh, Cell cell, int a, int b, int c)
{
  set_oh_length(oh);
  read_h_pos(cell);
  super_self(a, b, c);
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
  rp = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_env_setup();
  gsl_rng_set(rp, time(0));
}

int Ice::rng_int(int a)
{
  return gsl_rng_get(rp)%a;
}

double Ice::rng_uniform(void)
{
  return gsl_rng_uniform(rp);
}

// Populate each hydrogen bond with one hydrogen, pick one of two randoms sites
void Ice::populate_h_random(void)
{
  noccupied = natoms;
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
    noccupied--;
  }
}

// Fix a water molecule (don't let its configuration change)
void Ice::fix_water(int w)
{
  waters[w].fix();
  get_atom(waters[w].H1).fixed = true;
  get_atom(waters[w].H2).fixed = true;
  get_atom(waters[w].H3).fixed = true;
  get_atom(waters[w].H4).fixed = true;
  get_atom(waters[w].O).fixed = true;
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

// check that all hydrogen bonds are singly occupied
int Ice::hbond_single_occupied(void)
{
  int n_single_occ = 0;
  for (int i=0; i<nhbond; i++){
    if (hbond_occ(i) == 1) n_single_occ++;
  }
  return n_single_occ;
}

// count the number of two-coordinated oxygens
int Ice::o_two_coordinated(void)
{
  int n_two_coord = 0;
  for (int i=0; i<nwater; i++){
    if (water_coord(i) == 2){
      n_two_coord++;
    }
  }
  return n_two_coord;
}

// Swap the H position on a hydrogen bond
bool Ice::swap_h(int hb)
{
  int H1, H2;
  H1 = hbonds[hb].H1;
  H2 = hbonds[hb].H2;
  if (get_atom(H1).fixed || get_atom(H2).fixed){ 
    return false;
  }
  else {
    // get_atom(H1).toggle_occupied();
    // get_atom(H2).toggle_occupied();
    get_atom(H1).occupied = !get_atom(H1).occupied;
    get_atom(H2).occupied = !get_atom(H2).occupied;
    return true;
  }
}

bool Ice::rotate_water(int w)
{
  if (waters[w].fixed) return false;
  init_rng();
  int H1, H2;
  // pick two positions at random, one occupied and one unoccupied 
  H1 = get_random_h(w);
  bool h1_occ = get_atom(H1).occupied;
  bool h2_occ;
  bool done = false;
  while (!done){
    H2 = get_random_h(w);
    h2_occ = get_atom(H2).occupied;
    if (H1 == H2) continue;
    if (h1_occ){
      if (h2_occ){
        continue;
      } else {
        break;
      }
    }
    else {
      if (h2_occ){
        break;
      } else {
        continue;
      }
    }
  }
  get_atom(H1).toggle_occupied();
  get_atom(H2).toggle_occupied();
  return true;
}

// Monte Carlo algorithm to correct ice rules (Buch et al JCP B 102:8641, 1998)
void Ice::buch_mc_correct(void)
{
  std::cout << "Enforcing ice rules via Buch algorithm..." << std::endl;
  int c1a, c2a, cdiffa, c1b, c2b, cdiffb, iter;
  bool swapped;
  std::string trajfilename = "trajectory.xsf";

  init_rng();
  iter = 0;
  while (o_two_coordinated() != nwater){
    // Pick a random h-bond
    int hb_ind = rng_int(nhbond);
    // before the h swap
    c1a = water_coord(hbonds[hb_ind].W1);
    c2a = water_coord(hbonds[hb_ind].W2);
    cdiffa = std::abs(c1a - c2a);
    // after the h swap
    swapped = swap_h(hb_ind);
    c1b = water_coord(hbonds[hb_ind].W1);
    c2b = water_coord(hbonds[hb_ind].W2);
    cdiffb = std::abs(c1b - c2b);
    //No change in the difference
    if (cdiffb - cdiffa == 0){
      // move with probability 1/2
      double rn = rng_uniform();
      // if (rn < 0.5) swapped = swap_h(hb_ind);
      if (rn < 0.5) {
        swapped = swap_h(hb_ind);
      } else{
        if (flag_debug){
          write_xsf(trajfilename, iter, true);
          iter++;
        }
      }
    // Decrease in the coordination difference after move
    } else if (cdiffb - cdiffa < 0){
      // accept move
      if (flag_debug){
        write_xsf(trajfilename, iter, true);
        iter++;
      }
    // Increase in the coordination difference after move
    } else if (cdiffb - cdiffa > 0){
      // reject move
      swapped = swap_h(hb_ind);
    }
  }
  std::cout << "...done" << std::endl;
}

// Randomise the configuration, adding specified number of defects
void Ice::add_defects(int nBjerrum, int nIonic)
{
  int c1a, c2a, cdiffa, c1b, c2b, cdiffb;
  bool swapped;
  std::string trajfilename = "trajectory.xsf";
  if (nBjerrum > 1 || nIonic > 1){
    throw std::runtime_error("Addition of >1 defect pair not implemented yet!");
  }
  init_rng();
  std::cout << "Randomising configuration and adding defects:" << std::endl;
  if (nBjerrum > 0){
    std::cout << nBjerrum << " Bjerrum defect pairs" << std::endl;
    if (nIonic > 0){
      throw std::runtime_error("Mixture of defect types not implemented yet!");
    }
  }
  if (nIonic > 0){
    std::cout << nIonic << " Ionic defect pairs" << std::endl;
    if (nBjerrum > 0){
      throw std::runtime_error("Mixture of defect types not implemented yet!");
    }
  }

  double max_dist;
  int defect1, defect2, H1, H2;

  // Pick a random oxygen atom as the first defect
  defect1 = rng_int(nwater);
  int o1_ind = waters[defect1].O;
  max_dist = 0.0;
  // Find an oxygen sufficiently far away as the second
  std::string oxy = "O";
  for (int at=0; at<natoms; at++){
    if (atoms[at].name == oxy){
      if (dt(o1_ind,at) > max_dist){
        max_dist = dt(defect1, at);
        defect2 = at;
      }
    }
  }
  for (int i=0; i<nwater; i++){
    if (defect2 == waters[i].O){
      defect2 = i;
      break;
    }
  }

  std::string comment = "defect";
  int hb_ind, w_ind;
  int iter = 0;
  int maxiters = 100000;
  double rn;

  if (nBjerrum > 0){
    // Pick a random h-bond attached to the molecule, add Bjerrum D defect
    int rand_hb = rng_int(4);
    int hb1, hb2, O1, O2;
    hb1 = waters[defect1].hbonds[rand_hb];
    O1 = waters[defect1].O;
    H1 = hbonds[hb1].H1;
    H2 = hbonds[hb1].H2;
    if (flag_debug) std::cout << "Defect 1 indices H1: " << H1 \
      << " H2: " << H2 << std::endl;
    get_atom(H1).occupied = true;
    get_atom(H2).occupied = true;
    get_atom(H1).fixed = true;
    get_atom(H2).fixed = true;
    get_atom(O1).comment = comment;
    // Pick a random h-bond attached to the molecule, add Bjerrum L defect
    do {
      rand_hb = rng_int(4);
      hb2 = waters[defect2].hbonds[rand_hb];
    } while (hb2 == hb1);
    O2 = waters[defect2].O;
    H1 = hbonds[hb2].H1;
    H2 = hbonds[hb2].H2;
    if (flag_debug) std::cout << "Defect 2 indices H1: " << H1 \
      << " H2: " << H2 << std::endl;
    get_atom(H1).occupied = false;
    get_atom(H2).occupied = false;
    get_atom(H1).fixed = true;
    get_atom(H2).fixed = true;
    get_atom(O2).comment = comment;

    // use Buch algorithm with new constraints to enforce ice rules elsewhere
    std::cout << "Randomising configuration via Buch algorithm..."<< std::endl;
    while (o_two_coordinated() != nwater){
      // Pick a random h-bond
      hb_ind = rng_int(nhbond);
      // before the h swap
      c1a = water_coord(hbonds[hb_ind].W1);
      c2a = water_coord(hbonds[hb_ind].W2);
      cdiffa = std::abs(c1a - c2a);
      // after the h swap
      swapped = swap_h(hb_ind);
      if (!swapped) continue;
      c1b = water_coord(hbonds[hb_ind].W1);
      c2b = water_coord(hbonds[hb_ind].W2);
      cdiffb = std::abs(c1b - c2b);
      //No change in the difference
      if (cdiffb - cdiffa == 0){
        // move with probability 1/2
        rn = rng_uniform();
        if (rn < 0.5) {
          swapped = swap_h(hb_ind);
        } else {
          write_xsf(trajfilename, iter, true);
          iter++;
        }
      // Decrease in the coordination difference after move
      } else if (cdiffb - cdiffa < 0){
        // accept move
        write_xsf(trajfilename, iter, true);
        iter++;
      // Increase in the coordination difference after move 
      } else if (cdiffb - cdiffa > 0){
        // reject move
        swapped = swap_h(hb_ind);
      }
      if (iter >= maxiters) {
        std::cout << "Maximum iterations (" << maxiters  << ") reached." \
          << std::endl;
        break;
      }
    }
  } else if (nIonic > 0){
    bool rotated;
    int occ, n1occa, n1occb;
    Water w_tmp; 
    // Add a random H to defect1 to make H3O+ ion
    protonate_water(defect1);
    fix_water(defect1);
    // Remove a random H from defect2 to make OH- ion
    deprotonate_water(defect2);
    fix_water(defect2);
    while (hbond_single_occupied() != nhbond){
      // Pick a random water molecule
      w_ind = rng_int(nwater);
      w_tmp = waters[w_ind];
      // before the rotation
      n1occa = 0;
      for (int i=0; i<4; i++){
        occ = hbond_occ(waters[defect1].hbonds[i]);
        n1occa++;
      }
      // after the rotation
      rotated = rotate_water(w_ind);
      if (!rotated) continue; // don't swap if a molecule is fixed
      n1occb = 0;
      for (int i=0; i<4; i++){
        occ = hbond_occ(waters[defect2].hbonds[i]);
        n1occb++;
      }
      // Number of singly occupied Hbonds is the same
      if (n1occa == n1occb){
        // move with probability 1/2
        rn = rng_uniform();
        if (rn < 0.5) {
          waters[w_ind] = w_tmp; // revert the configuration
          iter++;
        } else {
          iter++;
          write_xsf(trajfilename, iter, true);
        }
      // Number of singly occupied H bonds increases
      } else if (n1occa < n1occb){
        // accept move
        iter++;
        write_xsf(trajfilename, iter, true);
      // Increase in the coordination difference after move
      } else if (n1occa > n1occb){
        // reject move
        // if (rn > 0.01) {
          waters[w_ind] = w_tmp; // revert the configuration
          iter++;
        // } else {
        //   iter++;
        //   write_xsf(trajfilename, iter, true);
        // }
      }
      if (iter >= maxiters) {
        std::cout << "Maximum iterations (" << maxiters  << ") reached." \
          << std::endl;
        break;
      }
    }
  }
  std::cout << "...done" << std::endl;
}

// Add a random hydrogen to water w
void Ice::protonate_water(int w)
{
  bool done = false;
  int new_h;
  while (!done){
    new_h = get_random_h(w);
    if (get_atom(new_h).occupied){
      continue;
    } else {
      get_atom(new_h).occupied = true;
      done = true;
    }
  }
}

// Remove a random hydrongen from water w
void Ice::deprotonate_water(int w)
{
  bool done = false;
  int new_h;
  while (!done){
    new_h = get_random_h(w);
    if (get_atom(new_h).occupied){
      get_atom(new_h).occupied = false;
      done = true;
    } else {
      continue;
    }
  }
}

// Pick a random H bond from water w
int Ice::get_random_h(int w)
{
  int rn = rng_int(4);
  int rand_h;
  switch (rn){
    case 0:
      rand_h = waters[w].H1;
      break;
    case 1:
      rand_h = waters[w].H2;
      break;
    case 2:
      rand_h = waters[w].H3;
      break;
    case 3:
      rand_h = waters[w].H4;
      break;
  }
  return rand_h;
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
std::deque<Node> Ice::get_loop(int w_start)
{
  std::deque<Node> loop;
  int current, hb, target;
  bool end;

  init_rng();
  current = w_start;
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
void Ice::rick_move(int w_start)
{
  std::deque<Node> loop;
  bool swapped;
  assert(w_start < nwater);
  do loop = get_loop(w_start);
  while (loop.size() < 3);

  for (int i=0; i<loop.size(); i++){
    swapped = swap_h(loop[i].hbond);
  }
}

// Randomise water orientations using Rick algorithm with C1 dipole contraint
void Ice::rick_randomise(int max_loops)
{
  std::cout << "Randomising water orientations via Rick algorithm..." 
            << std::endl;
  int ndefect, nmove, naccepted, loop_start;
  double cell_dipole, cell_dipole_old, rn, mcp;
  std::vector<bool> conf;
  unsigned int loop = 0;
  boost::progress_display show_progress(loop);

  init_rng();

  cell_dipole_old = c1_dipole();
  conf = save_config();
  nmove = 0;
  naccepted = 0;
  for (loop=0; loop<max_loops; loop++){
    while (true){
      while (true){
        loop_start = rng_int(nwater);   // Pick a randomw water to start loop
        rick_move(loop_start); 
        nmove++;
        ndefect = check_ionic_defects();
        if (ndefect == 0) break;
        else revert_config(conf);
      }
      // Monte Carlo criterion
      cell_dipole = c1_dipole();
//      if (cell_dipole > cell_dipole_old){
      mcp = exp(cell_dipole_old - cell_dipole);
      rn = rng_uniform();
      if (rn > mcp) revert_config(conf);
//      }
      else break;
    }
    conf = save_config();
    if (flag_debug){
      std::cout << "Iteration " << loop << ": Cell dipole = " 
                << cell_dipole_old << std::endl;
    }
    else ++show_progress;
    if (cell_dipole < cell_dipole_thresh) break;
    cell_dipole_old = cell_dipole;
  }
  std::cout << "Total number of moves = " << nmove << std::endl;
  std::cout << "Finished Rick algorithm. Final dipole = " << cell_dipole 
            << std::endl;
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

void Ice::check_defects(void)
{
  std::cout << "Checking ice rules" << std::endl;
  int nionic = check_ionic_defects();
  int nbjerrum = check_bjerrum_defects();
  int n_two_coord = o_two_coordinated();
  int n_good_hbonds = hbond_single_occupied();

  std::cout << n_two_coord << " / " << nwater \
            << " two-coordinated water molecules" << std::endl;
  std::cout << n_good_hbonds << " / " << nhbond \
            << " singly-occupied hydrogen bonds" << std::endl;

  if (nionic != 0){
    std::cout << "WARNING: " << nionic << " ionic defects" << std::endl;
  }
  if (nbjerrum != 0){
    std::cout << "WARNING: " << nbjerrum << " Bjerrum defects" << std::endl;
  }
}

// Get the dipole vector of water w (multiply orientation vector by constant
// 3.09 -- Batista et al JCP 109, 4546 (1998))
Eigen::Vector3d Ice::water_dipole(int w)
{
  Eigen::Vector3d oh1, oh2;

  int nh = 0;
  int io = waters.at(w).O;
  int ih1, ih2;
  assert(water_coord(w) == 2);
  if (get_atom(waters.at(w).H1).occupied) {
    if (nh == 0) ih1 = waters.at(w).H1;
    if (nh == 1) ih2 = waters.at(w).H1;
    nh++;
  }
  if (get_atom(waters.at(w).H2).occupied) {
    if (nh == 0) ih1 = waters.at(w).H2;
    if (nh == 1) ih2 = waters.at(w).H2;
    nh++;
  }
  if (get_atom(waters.at(w).H3).occupied) {
    if (nh == 0) ih1 = waters.at(w).H3;
    if (nh == 1) ih2 = waters.at(w).H3;
    nh++;
  }
  if (get_atom(waters.at(w).H4).occupied) {
    if (nh == 0) ih1 = waters.at(w).H4;
    if (nh == 1) ih2 = waters.at(w).H4;
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

// Construct a slab in the defined direction (perpendicular to bilayers)
// Put bilayer 1 at the bottom of the slab
void Ice::build_slab(double dhkl, int direction)
{
  std::cout << "Constructing a slab in direction " << direction 
            << " with dhkl " << dhkl << std::endl;
  wrap();
  if (frac) frac2cart_all();
  // Assign each water molecule to a bilayer
  int bilayer;
  double r;
  nbilayer = static_cast<int>(round(lat(direction,direction)/dhkl));

  for (int i=0; i<waters.size(); i++){
    r = get_atom(waters[i].O).r(direction);
    bilayer = static_cast<int>(ceil(r/dhkl));
    waters[i].bilayer = bilayer;
  }

  int ibilayer1 = 0;
  for (int i=0; i<nwater; i++){
    if (waters[i].bilayer == 1){
      ibilayer1 = i;
      break;
    }
  }

  // Shift all atoms such that bilayer 1 is at the bottom of the cell
  // double rbilayer1 = get_atom(waters[ibilayer1].O).r(direction);
  // double lshift = -rbilayer1 + dhkl/2.0;
  // switch(direction){
  //   case 0:
  //     shift(lshift, 0.0, 0.0);
  //   case 1:
  //     shift(0.0, lshift, 0.0);
  //   case 2:
  //     shift(0.0, 0.0, lshift);
  // }

  // label the surface water molecules
  double ormin = lat(direction, direction);
  double ormax = 0.0;
  double thresh = 0.1;
  for (int i=0; i<nwater; i++){
    r = get_atom(waters[i].O).r(direction);
    if (r < ormin) ormin = r;
    if (r > ormax) ormax = r;
  }
  for (int i=0; i<nwater; i++){
    r = get_atom(waters[i].O).r(direction);
    if (r < (ormin+thresh)){
      waters[i].surface1 = true; 
      waters[i].step = 0;
      s1list.push_back(i);
    }
    else waters[i].surface1 = false;
    if (r > (ormax-thresh)){
      waters[i].surface2 = true;
      s2list.push_back(i);
    }
    else waters[i].surface2 = false;
  }
}

// Find and label dangling H
std::deque<int> Ice::find_dOH(int direction)
{
  // std::cout << "Finding dangling OH bonds..." << std::endl;
  // reference vector to compare direction
  Eigen::Vector3d ref;
  switch (direction){
    case 0:
      ref << 1.0, 0.0, 0.0;
    case 1:
      ref << 0.0, 1.0, 0.0;
    case 2:
      ref << 0.0, 0.0, 1.0;
  }

  int w;
  std::deque<int> dOHlist;
  Eigen::Vector3d oh1, oh2, oh3, oh4;
  for (int i=0; i<s1list.size(); i++){
    w = s1list[i];
    waters[w].dOH = false;
    oh1 = mic_cart(get_atom(waters[w].O), get_atom(waters[w].H1));
    oh2 = mic_cart(get_atom(waters[w].O), get_atom(waters[w].H2));
    oh3 = mic_cart(get_atom(waters[w].O), get_atom(waters[w].H3));
    oh4 = mic_cart(get_atom(waters[w].O), get_atom(waters[w].H4));
    if (get_atom(waters[w].H1).occupied){
      if (oh1.cross(ref).norm() < smallish){
        waters[w].dOH = true;
        dOHlist.push_back(w);
      }
    }
    if (get_atom(waters[w].H2).occupied){
      if (oh2.cross(ref).norm() < smallish){
        waters[w].dOH = true;
        dOHlist.push_back(w);
      }
    }
    if (get_atom(waters[w].H3).occupied){
      if (oh3.cross(ref).norm() < smallish){
        waters[w].dOH = true;
        dOHlist.push_back(w);
      }
    }
    if (get_atom(waters[w].H4).occupied){
      if (oh4.cross(ref).norm() < smallish){
        waters[w].dOH = true;
        dOHlist.push_back(w);
      }
    }
  }  
  for (int i=0; i<s2list.size(); i++){
    w = s2list[i];
    waters[w].dOH = false;
    oh1 = mic_cart(get_atom(waters[w].O), get_atom(waters[w].H1));
    oh2 = mic_cart(get_atom(waters[w].O), get_atom(waters[w].H2));
    oh3 = mic_cart(get_atom(waters[w].O), get_atom(waters[w].H3));
    oh4 = mic_cart(get_atom(waters[w].O), get_atom(waters[w].H4));
    if (get_atom(waters[w].H1).occupied){
      if (oh1.cross(ref).norm() < small){
        waters[w].dOH = true;
        dOHlist.push_back(w);
      }
    }
    if (get_atom(waters[w].H2).occupied){
      if (oh2.cross(ref).norm() < small){
        waters[w].dOH = true;
        dOHlist.push_back(w);
      }
    }
    if (get_atom(waters[w].H3).occupied){
      if (oh3.cross(ref).norm() < small){
        waters[w].dOH = true;
        dOHlist.push_back(w);
      }
    }
    if (get_atom(waters[w].H4).occupied){
      if (oh4.cross(ref).norm() < small){
        waters[w].dOH = true;
        dOHlist.push_back(w);
      }
    }
  }  
  return dOHlist;
}

// Compute surface order parameter
double Ice::order_parameter(double surface_nn_cutoff)
{
  if (flag_debug){
    std::cout << "Computing order parameter" << std::endl;
  }
  int s_nn_total = 0; // running total of dangling OH surface neighbours 
  int ndOH = 0;
  int w1, w2, o1, o2;

  for (int i=0; i<s1list.size(); i++){
    w1 = s1list[i];
    if (!waters[w1].dOH) continue;
    else ndOH++;
    o1 = waters[w1].O;
    for(int j=0; j<s1list.size(); j++){
      w2 = s1list[j];
      o2 = waters[w2].O;
      if (i != j){
        if (dt(o1, o2) < surface_nn_cutoff){
          if (waters[w2].dOH) s_nn_total++;
        }
      } 
    } 
  }

  for (int i=0; i<s2list.size(); i++){
    w1 = s2list[i];
    if (!waters[w1].dOH) continue;
    else ndOH++;
    o1 = waters[w1].O;
    for(int j=0; j<s2list.size(); j++){
      w2 = s2list[j];
      o2 = waters[w2].O;
      if (i != j){
        if (dt(o1, o2) < surface_nn_cutoff){
          if (waters[w2].dOH) s_nn_total++;
        }
      } 
    } 
  }
  // std::cout << s_nn_total << " " << ndOH << std::endl;

  return static_cast<double>(s_nn_total)/static_cast<double>(ndOH);
}

void Ice::build_ordered_slab(double dhkl, int direction, double target_cOH, int max_loops)
{
  build_slab(dhkl, direction);
  std::cout << "Constructing slab via Rick algorithm" << std::endl;
  std::cout << "Target c_OH = " << target_cOH << std::endl;

  std::deque<int> dOHlist;
  int ndefect, nmove, naccepted, loop_start;
  double cell_dipole, cell_dipole_old, rn, mcp, cOH, cOH_diff, cOH_diff_old;
  std::vector<bool> conf;

  init_rng();

  cell_dipole_old = c1_dipole();
  dOHlist = find_dOH(direction);
  cOH = order_parameter(surface_nn_cut);
  cOH_diff_old = std::abs(cOH - target_cOH);
  conf = save_config();
  nmove = 0;
  naccepted = 0;
  for (int i=0; i<max_loops; i++){
    while (true){
      while (true){
        // Start the loop from the surface to force cOH to change more 
        // frequently
        rn = rng_uniform();
        if (rn < 0.5) loop_start = s1list[rng_int(s1list.size())];
        else loop_start = s2list[rng_int(s2list.size())];
        rick_move(loop_start); 
        nmove++;
        ndefect = check_ionic_defects();
        if (ndefect == 0) break;
        else revert_config(conf);
      }
      // Monte Carlo criterion
      cell_dipole = c1_dipole();
      dOHlist = find_dOH(direction);
      cOH = order_parameter(surface_nn_cut);
      cOH_diff = std::abs(cOH - target_cOH);

      mcp = exp(cell_dipole_old - cell_dipole);
      rn = rng_uniform();
      if (rn > mcp) revert_config(conf);
      else {
        mcp = exp((cOH_diff_old - cOH_diff)*50.0);
        rn = rng_uniform();
        // if (rn > 0.5) revert_config();
        if (rn > mcp) revert_config(conf);
        else break;
      }
    }
    conf = save_config();
    if (flag_debug) {
      std::cout << "Iteration " << i << std::endl;
      std::cout << "  Cell dipole = " << cell_dipole_old << std::endl;
      std::cout << "  cOH         = " << cOH << std::endl;
      cell_dipole_old = cell_dipole;
      cOH_diff_old = cOH_diff;
      if (cell_dipole < cell_dipole_thresh){
        if (cOH_diff <= cOH_thresh) break;
      }
      if (i == max_loops-1){
        std::cout << "Exceeded maximum iterations (" << max_loops << ")" 
                  << std::endl;
      }
    }
  }
  std::cout << "Finished Rick algorithm. " << std::endl;
  std::cout << "Total number of moves = " << nmove << std::endl;
  std::cout << "  Final dipole = " << cell_dipole << std::endl;
  std::cout << "  Final cOH     = " << cOH << std::endl;
}

// Construct a step: compute which atoms to remove, and move them to 
// the bottom of the file
void Ice::build_step(std::string direction, int layers, 
                     double step_width, double vacuum_gap,
                     bool oneside, std::string fname)
{
  if (flag_debug){
    std::cout << "Building step..." << std::endl;
  }
  int dir = 0;
  std::string t = "T";

  if (frac) frac2cart_all();

  if (direction == "a") dir = 0;
  else if (direction == "b") dir = 1;
  else if (direction == "c") dir = 2;
  else std::cout << "Invalid direction" << std::endl;

  int iO, iH1, iH2, iH3, iH4, layer;
  double valley_width = (lat(dir,dir) - step_width)/2.0;
  double bound_incr = 0.5;
  double bound1 = valley_width;
  double bound2 = step_width + valley_width;
  double binwidth = 0.0;
  binwidth = lat(dir,dir)/static_cast<double>(layers);

  bound1 = bound_incr;
  bound2 = lat(dir,dir) - bound_incr;
  boost::format fmt("%03d");
  for (int i=0; i<nwater; i++){
    if (waters[i].bilayer == nbilayer) {
      std::string tag = " # step ";
      iO = waters[i].O;
      iH1 = waters[i].H1;
      iH2 = waters[i].H2;
      iH3 = waters[i].H3;
      iH4 = waters[i].H4;
      layer = ceil(atoms[iO].r(dir)/binwidth);
      fmt%layer;
      tag += fmt.str();
      // tag += std::to_string(layer);
      atoms[iO].comment = tag;
      atoms[iH1].comment = tag;
      atoms[iH2].comment = tag;
      atoms[iH3].comment = tag;
      atoms[iH4].comment = tag;
    } 
    if (!oneside){
      if (waters[i].bilayer == 1) {
        std::string tag = " # step ";
        iO = waters[i].O;
        iH1 = waters[i].H1;
        iH2 = waters[i].H2;
        iH3 = waters[i].H3;
        iH4 = waters[i].H4;
        layer = ceil(atoms[iO].r(dir)/binwidth);
        tag += std::to_string(layer);
        atoms[iO].comment = tag;
        atoms[iH1].comment = tag;
        atoms[iH2].comment = tag;
        atoms[iH3].comment = tag;
        atoms[iH4].comment = tag;
      }
    }
  }

  std::string file_bulk = fname + "_bulk.in";
  std::string file_slab = fname + "_slab.in";
  std::string file_step = fname + "_step.in";
  std::string atom_str;
  int species_label;

  std::ofstream ofs;

  ofs.open(file_bulk);
  if (ofs.fail()) {
    std::cout << "Error opening file " << file_bulk << std::endl;
  }
  for (int i=0; i<=2; i++) {
    ofs << std::fixed << std::setprecision(8) \
        << lat.row(i)*ang2bohr << std::endl;
  }
  ofs << noccupied << std::endl;
  for (int i=0; i<natoms; i++){
    if (atoms[i].occupied){
      if (!atoms[i].remove){
        if (atoms[i].name == "H") species_label = 1;
        else if (atoms[i].name == "O") species_label = 2;
        atom_str = atoms[i].write_cq(species_label, t, t, t);
        ofs << atom_str << std::endl;
      }
    }
  }
  for (int i=0; i<natoms; i++){
    if (atoms[i].occupied){
      if (atoms[i].remove){
        if (atoms[i].name == "H") species_label = 1;
        else if (atoms[i].name == "O") species_label = 2;
        atom_str = atoms[i].write_cq(species_label, t, t, t);
        ofs << atom_str << std::endl;
      }
    }
  }
  ofs.close();

  ofs.open(file_slab);
  if (ofs.fail()) {
    std::cout << "Error opening file " << file_slab << std::endl;
  }
  Eigen::Vector3d lat_c = lat.row(2);
  lat_c(2) += vacuum_gap;
  lat_c *= ang2bohr;
  ofs << std::fixed << std::setw(14) << std::setprecision(8) \
      << std::left << lat.row(0)*ang2bohr << std::endl \
      << std::left << lat.row(1)*ang2bohr << std::endl \
      << std::left << std::setw(14) << lat_c(0) \
      << std::left << std::setw(14) << lat_c(1) \
      << std::left << std::setw(14) << lat_c(2) << std::endl;
  ofs << noccupied << std::endl;
  for (int i=0; i<natoms; i++){
    if (atoms[i].occupied){
      if (!atoms[i].remove){
        if (atoms[i].name == "H") species_label = 1;
        else if (atoms[i].name == "O") species_label = 2;
        atom_str = atoms[i].write_cq(species_label, t, t, t);
        // atom_str = atoms[i].write_cell();
        ofs << atom_str << std::endl;
      }
    }
  }
  for (int i=0; i<natoms; i++){
    if (atoms[i].occupied){
      if (atoms[i].remove){
        if (atoms[i].name == "H") species_label = 1;
        // else if (atoms[i].name == "O") species_label = 2;
        atom_str = atoms[i].write_cq(species_label, t, t, t);
        atom_str = atoms[i].write_cell();
        ofs << atom_str << std::endl;
      }
    }
  }
  ofs.close();
}


// Write object to a .xsf file
void Ice::write_xsf(std::string fname, int iter, bool append)
{
  std::vector<std::string> species;
  std::ofstream outfile;
  if (iter==0){
    outfile.open(fname.c_str());
  } else{
    outfile.open(fname.c_str(), std::ios_base::app);
  }

  if (flag_debug && iter==0){
    std::cout << "Writing ice configuration to " << fname \
      << " (xsf)" << std::endl;
  }

  if (iter==0){
    outfile << "ANIMSTEP " << std::endl;
  }
  outfile << "CRYSTAL" << std::endl;
  outfile << "PRIMVEC   " << iter << std::endl;
  for (int i=0; i<=2; i++) {
    outfile << std::fixed << std::right << std::setprecision(8)
            << lat.row(i) << std::endl;
  }
  outfile << "PRIMCOORD " << iter << std::endl;
  outfile << nwater*3 << " 1" << std::endl;
  for (int i=0; i<natoms; i++) {
    if (atoms[i].occupied){
      outfile << get_atom(i) << std::endl;
    }
  }
}

// Write object to a .cell file
void Ice::write_cell(std::string fname)
{
  if (flag_debug){
    std::cout << "Writing ice configuration to " << fname << std::endl;
  }
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

  std::string comment = "defect";
  std::string newname = "N";
  for (int i=0; i<nwater; i++) {
    // outfile << get_atom(waters[i].O) << std::endl;
    outfile << get_atom(waters[i].O).write_cell_highlight(newname, comment) << std::endl;
    if (atoms[waters[i].H1].occupied){
      // outfile << get_atom(waters[i].H1).write_cell_highlight(newname) << std::endl;
      outfile << get_atom(waters[i].H1) << std::endl;
    }
    if (atoms[waters[i].H2].occupied){
      // outfile << get_atom(waters[i].H2).write_cell_highlight(newname) << std::endl;
      outfile << get_atom(waters[i].H2) << std::endl;
    }
    if (atoms[waters[i].H3].occupied){
      // outfile << get_atom(waters[i].H3).write_cell_highlight(newname) << std::endl;
      outfile << get_atom(waters[i].H3) << std::endl;
    }
    if (atoms[waters[i].H4].occupied){
      // outfile << get_atom(waters[i].H4).write_cell_highlight(newname) << std::endl;
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
