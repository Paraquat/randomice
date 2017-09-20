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
  // get_nn(oo_max);
  // for (int i=0; i<natoms; i++){
  //   for (int j=0; j<atoms[i].nn.size(); j++){
  //     oo = mic_cart(atoms[i], *atoms[i].nn[j]);
  //     r_h = atoms[i].r + oh_default*oo.normalized();
  //     Atom a = Atom(species, natoms+1, r_h, occ);
  //     add_atom(a);
  //   }
  // }

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
  std::cout << "Finding water molecules" << std::endl;
  nwater = 0;
  get_nn(oh_max);
  for (int i=0; i<natoms; i++){
    if (atoms[i].name == "O"){
      assert(atoms[i].nn.size() == 4);
      Atom::atom_ptr o = std::make_shared<Atom>(atoms[i]);
      Atom::atom_ptr h1 = std::make_shared<Atom>(*atoms[i].nn[0]);
      Atom::atom_ptr h2 = std::make_shared<Atom>(*atoms[i].nn[1]);
      Atom::atom_ptr h3 = std::make_shared<Atom>(*atoms[i].nn[2]);
      Atom::atom_ptr h4 = std::make_shared<Atom>(*atoms[i].nn[3]);
      // Atom::atom_ptr o(new Atom(atoms[i]));
      // Atom::atom_ptr h1(new Atom(*atoms[i].nn[0]));
      // Atom::atom_ptr h2(new Atom(*atoms[i].nn[1]));
      // Atom::atom_ptr h3(new Atom(*atoms[i].nn[2]));
      // Atom::atom_ptr h4(new Atom(*atoms[i].nn[3]));

      Water w(o->label, o, h1, h2, h3, h4);
      add_water(w);
      // std::cout << waters.back() << std::endl;
    }
  }

  // get_dt();
  // for (int i=0; i<natoms; i++){
  //   if (atoms[i].name == "O"){
  //     Atom::atom_ptr o = std::make_shared<Atom>(atoms[i]);
  //     Atom::atom_ptr h1, h2, h3, h4;
  //     int hcount = 0;
  //     for (int j=0; j<natoms; j++){
  //       if (dt(i,j) < oh_max){
  //         if (atoms[j].name == "H"){
  //           hcount++;
  //           if (hcount == 1){
  //             h1 = std::make_shared<Atom>(atoms[j]);
  //           } else if (hcount == 2){
  //             h2 = std::make_shared<Atom>(atoms[j]);
  //           } else if (hcount == 3){
  //             h3 = std::make_shared<Atom>(atoms[j]);
  //           } else if (hcount == 4){
  //             h4 = std::make_shared<Atom>(atoms[j]);
  //           }
  //         }
  //       } 
  //     }
  //     if (hcount != 4){
  //       throw std::runtime_error("Did not find 4 hydrogens for this water");
  //     }
  //     Water w(o->label, o, h1, h2, h3, h4);
  //     add_water(w);
  //   }
  // }
}

void Ice::get_water_nn(double cutoff)
{
  Eigen::Vector3d d;

  assert(frac == false);
  for (int i=0; i<nwater; i++){
    for (int j=0; j<nwater; j++){
      if (i != j){
        d = mic_cart(*waters[i].O, *waters[j].O);
        if (d.norm() < cutoff){
          Water::water_ptr ptr(new Water(waters[j]));
          waters[i].add_nn(ptr);
        }
      }
    }
  }
}

void Ice::get_hbonds(void)
{
  std::cout << "Finding hydrogen bonds" << std::endl;
  // bool done, adjacent;
  bool done;

  nhbond = 0;
  get_dt();

  // for (int i=0; i<natoms; i++){
  //   if (atoms[i].name == "O"){
  //     Atom::atom_ptr o1 = std::make_shared<Atom>(atoms[i]);
  //     for (int j=i+1; j<natoms; j++){
  //       if ((atoms[j].name == "O") && (dt(i,j) < oo_max)){
  //         Atom::atom_ptr o2 = std::make_shared<Atom>(atoms[j]);

  //         std::set<int> ws;
  //         ws.insert(atoms[i].label);
  //         ws.insert(atoms[j].label);

  //         done = false;
  //         for (int k=0; k<hbonds.size(); k++){
  //           if (ws == hbonds[k].w){
  //             done = true;
  //             break;
  //           }
  //         }

  //         if (done) continue;
  //         // Atom::atom_ptr h1, h2;
  //         int nh = 0; 
  //         for (int k=0; k<natoms; k++){
  //           if (dt(i,k) < oh_max) adjacent = true;
  //           else if (dt(j,k) < oh_max) adjacent = true;
  //           else adjacent = false;
  //           if (adjacent == false) continue;

  //           if (atoms[k].name == "H"){
  //             if (isPointOnLine(atoms[i], atoms[j], atoms[k])){
  //               nh++;
  //               if (nh == 1){
  //                 h1_ind = k;
  //                 // h1.reset(new Atom(atoms[k]));
  //                 // h1 = std::make_shared<Atom>(atoms[k]);
  //               } else if (nh == 2){
  //                 h2_ind = k;
  //                 // h2.reset(new Atom(atoms[k]));
  //                 // h2 = std::make_shared<Atom>(atoms[k]);
  //               }
  //             }
  //           }
  //         } 
  //         assert(nh == 2);
  //         Atom::atom_ptr h1(new Atom(atoms[h1_ind]));
  //         Atom::atom_ptr h2(new Atom(atoms[h2_ind]));
  //         // Atom::atom_ptr h2 = std::make_shared<Atom>(atoms[h2_ind]);
  //         Hbond hb(o1, o2, h1, h2);
  //         add_hbond(hb);
  //       }
  //     }
  //   }
  // }

  get_nn(oo_max);
  for (int i=0; i<nwater; i++){
    // water_ptr w1(new Water(waters[i]));
    water_ptr w1 = std::make_shared<Water>(waters[i]);
    assert(waters[i].nn.size() == 4);
    for (int j=0; j<4; j++){
      // water_ptr w2(new Water(*waters[i].nn[j]));
      water_ptr w2 = std::make_shared<Water>(*waters[i].nn[j]);

      // check whether this hbond is already in the list
      std::set<int> ws;
      ws.insert(w1->label);
      ws.insert(w2->label);

      done = false;
      for (int k=0; k<hbonds.size(); k++){
        if (ws == hbonds[k].w){
          done = true;
          break;
        }
      }
      if (done) continue;

      Atom::atom_ptr h1, h2;
      if (isPointOnLine(*w1->O, *w2->O, *w1->H1)){
        h1 = w1->H1;
      }
      else if (isPointOnLine(*w1->O, *w2->O, *w1->H2)){
        h1 = w1->H2;
      }
      else if (isPointOnLine(*w1->O, *w2->O, *w1->H3)){
        h1 = w1->H3;
      }
      else if (isPointOnLine(*w1->O, *w2->O, *w1->H4)){
        h1 = w1->H4;
      }
      else std::cout << "Failed to find h1 for water " << w1->label 
                     << std::endl;

      if (isPointOnLine(*w1->O, *w2->O, *w2->H1)){
        h2 = w2->H1;
      }
      else if (isPointOnLine(*w1->O, *w2->O, *w2->H2)){
        h2 = w2->H2;
      }
      else if (isPointOnLine(*w1->O, *w2->O, *w2->H3)){
        h2 = w2->H3;
      }
      else if (isPointOnLine(*w1->O, *w2->O, *w2->H4)){
        h2 = w2->H4;
      }
      else std::cout << "Failed to find h2 for water " << w2->label 
                     << std::endl;

      Hbond hb((*w1).O, (*w2).O , h1, h2, w1, w2);

      add_hbond(hb);
      hbond_ptr hbp = std::make_shared<Hbond>(hbonds.back());
      w1->add_hbond(hbp);
      w2->add_hbond(hbp);
    }
  }
}

void Ice::print_ice(void)
{
  std::cout << nwater << " water molecules" << std::endl << std::endl;
  for (int i=0; i<nwater; i++){
    std::cout << *waters[i].O << std::endl;
    std::cout << *waters[i].H1 << std::endl;
    std::cout << *waters[i].H2 << std::endl;
    std::cout << *waters[i].H3 << std::endl;
    std::cout << *waters[i].H4 << std::endl;
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
    if (rn == 0) hbonds[i].H1->occupy(true);
    else hbonds[i].H2->occupy(true);
  }
}

// count the number of two-coordinated oxygens
int Ice::o_two_coordinated(void)
{
  int ntwocoord = 0;
  for (int i=0; i<nwater; i++){
    if (waters[i].coord() == 2){
      ntwocoord++;
    }
  }
  return ntwocoord;
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
    c1a = hbonds[hb_ind].W1->coord();
    c2a = hbonds[hb_ind].W2->coord();
    cdiffa = std::abs(c1a - c2a);
    // after the h swap
    hbonds[hb_ind].swap_h();
    c1b = hbonds[hb_ind].W1->coord();
    c2b = hbonds[hb_ind].W2->coord();
    cdiffb = std::abs(c1b - c2b);
    //No change in the difference
    if (cdiffb - cdiffa == 0){
      // move with probability 1/2
      double rn = rng_uniform();
      if (rn < 0.5) hbonds[hb_ind].swap_h();
    // Decrease in the coordination difference after move
    } else if (cdiffb - cdiffa < 0){
      // accept move

    // Increase in the coordination difference after move
    } else if (cdiffb - cdiffa > 0){
      // reject move
      hbonds[hb_ind].swap_h();
    }
  }
}


// Write object to a .cell file
void Ice::write_cell(std::string fname)
{
  std::ofstream     outfile(fname.c_str());

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
    if (atoms[i].occupied){
      outfile << waters[i];
    }
  }

  if (frac) outfile << "%ENDBLOCK positions_frac" << std::endl << std::endl;
  else outfile << "%ENDBLOCK positions_abs" << std::endl << std::endl;
}