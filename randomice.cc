#include "randomice.h"

namespace po = boost::program_options;

bool flag_debug;

void handler(int sig) {
  void *array[20];
  size_t size;

  // get void*s for all entries on the stack
  size = backtrace(array, 20);

  // print out all frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  exit(1);
}

int main(int argc, char* argv[]){

  po::options_description desc("Allowed options");
  std::string infname, phase;
  std::vector<int> scdim;
  int scx, scy, scz, nBjerrum, nIonic, slices;
  int maxiter = 1000;
  bool cq_out = false;
  bool ordered = false;
  bool defect = false;
  bool oneside = false;
  bool castep_out = false;
  bool vasp_out = false;
  bool bulk = true;
  bool slab = false;
  bool step = false;
  bool frac = false;
  bool randomise = false;
  double oh = oh_default;

  signal(SIGSEGV, handler);

  desc.add_options()
    ("help,h", "Print help information")
    ("infile,i", po::value< std::string >(), "Input structure file")
    ("phase,p", po::value< std::string >(), "Input ice phase (Ih, Ic, etc")
    ("debug,d", po::bool_switch(&flag_debug), "Run in debug mode")
    ("maxiter,m", po::value< int >(),
     "Maximum number of Rick algorithm iterations")
    ("cq", po::bool_switch(&cq_out), "Write to Conquest file")
    ("ordered,o", po::bool_switch(&ordered), "Generate an ordered supercell")
    ("defect", po::bool_switch(&defect), "Generate a disordered supercell containing a defect pair")
    ("oneside", po::bool_switch(&oneside), "Build step on only one side of slab")
    ("cell", po::bool_switch(&castep_out), "Write to CASTEP cell file")
    ("vasp", po::bool_switch(&vasp_out), "Write to VASP POSCAR file")
    ("bulk", po::bool_switch(&step), "Construct bulk ice")
    ("slab", po::bool_switch(&step), "Construct a slab")
    ("step", po::bool_switch(&step), "Construct a step")
    ("frac", po::bool_switch(&frac), "Output in fractional coordinates")
    ("step_direction", po::value< std::string >(), "Direction of step (a,b,c)")
    ("OH", po::value< double >(), "OH bond length")
    ("step_width", po::value< double >(), "Width of step (a0)")
    ("vacuum_gap", po::value< double >(), "Size of vacuum gap (a0)")
    ("supercell,s", po::value<std::vector<int> >() -> multitoken(), 
     "Supercell dimensions (a x b x c)")
    ("randomise,r", po::bool_switch(&randomise), "Randomise an existing cell")
    ("nBjerrum", po::value< int >(), "Number of Bjerrum defect pairs")
    ("nIonic", po::value< int >(), "Number of ionic defect pairs")
    ("slices", po::value< int >(), "Number of steps of varying d to generate")
    ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
  if (vm.count("help")){
    std::cout << desc << "\n";
    return 1;
  }
  if (vm.count("infile")){
    infname = vm["infile"].as<std::string>();
    if (flag_debug) std::cout << "Input file is " <<  infname << std::endl;
  } else {
    std::cout << "No input structure file specified" << std::endl;
  }
  if (vm.count("phase")){
    phase = vm["phase"].as<std::string>();
  } else {
    phase = "Ih";
  }
  if (!vm["supercell"].empty() &&
      (scdim = vm["supercell"].as<std::vector<int> >()).size() == 3){
    scx = scdim[0];
    scy = scdim[1];
    scz = scdim[2];
  } else {
    scx = 1;
    scy = 1;
    scz = 1;
  }
  if (vm.count("maxiter")){
    maxiter = vm["maxiter"].as<int>();
  }
  if (vm.count("OH")){
    oh = vm["OH"].as<double>();
  }
  if (vm.count("nBjerrum")){
    nBjerrum = vm["nBjerrum"].as<int>();
  } else {
    nBjerrum = 0;
  }
  if (vm.count("nIonic")){
    nIonic = vm["nIonic"].as<int>();
  } else {
    nIonic = 0;
  }
  if (vm.count("slices")){
    slices = vm["slices"].as<int>();
  } else {
    slices = 0;
  }

  if (bulk) std::cout << "Generating bulk ice" << std::endl;
  if (slab) std::cout << "Generating ice slab" << std::endl;
  if (step) std::cout << "Generating ice step" << std::endl;
  if (nBjerrum > 0) std::cout << "Adding " << nBjerrum \
    << " Bjerrum defect pair" << std::endl;
  if (nIonic > 0) std::cout << "Adding " << nIonic \
    << " Ionic defect pair" << std::endl;

  double step_width = 0.0;
  double vacuum_gap = 0.0;
  std::string step_direction = "a";
  if (step){
    if (vm["step_width"].empty()){
      step_width = 0.0;
    } else {
      step_width = vm["step_width"].as<double>();
    }
    vacuum_gap = vm["vacuum_gap"].as<double>();
    step_direction = vm["step_direction"].as<std::string>();
  }

  Cell cell;

  // Generate the bulk ice supercell
  cell.read_cell(infname);
  double dhkl = cell.lat(2,2);
  if (phase == "Ih"){
    dhkl = cell.lat(2,2)/2.0;
  }
  else if (phase == "Ic"){
    dhkl = cell.lat(2,2)/3.0;
  }
  Ice* ice;

  if (ordered){
    std::cout << "Reading initial hydrogen positions from file" << std::endl;
    ice = new Ice();
    ice -> init_bulk_ordered(oh, cell, scx, scy, scz);
  } else if (defect){
    std::cout << "Generating bulk ice with defect pair(s)" << std::endl;
    Cell sc = cell.super(scx, scy, scz);
    ice = new Ice(sc);
    ice -> init_bulk_random(oh);
    ice -> add_defects(nBjerrum, nIonic);
  } else{
    std::cout << "Generating proton disordered bulk ice Ih" << std::endl;
    Cell sc = cell.super(scx, scy, scz);
    ice = new Ice(sc);
    ice -> init_bulk_random(oh);
  }
  ice -> check_defects();
  if (frac) ice -> cart2frac_all();
  if (cq_out) ice -> write_cq("iceIh.coord");
  else if (vasp_out) ice -> write_vasp("iceIh.vasp");
  else if (castep_out) ice -> write_cell("iceIh.cell");

  // if (ordered){
  //   std::cout << "Reading initial hydrogen positions from file" << std::endl;
  //   Ice unit;
  //   unit.set_oh_length(oh);
  //   unit.read_h_pos(cell);
  //   ice.set_oh_length(oh);
  //   Ice* ice_sc = new Ice:
  //   ice = unit.super(scx, scy, scz);

  // } else if (defect){
  //   std::cout << "Generating bulk ice with defect pair(s)" << std::endl;
  //   Cell sc;
  //   sc = cell.super(scx, scy, scz);
  //   // Ice ice(sc);
  //   Ice* ice_defect = new Ice(sc);
  //   ice.set_oh_length(oh);
  //   ice.get_h_pos();
  //   ice.get_waters();
  //   ice.get_hbonds();

  //   std::cout << ice.nwater << " waters" << std::endl;
  //   std::cout << ice.nhbond << " hydrogen bonds" << std::endl;

  //   ice.populate_h_random();
  //   ice.add_defects(nBjerrum, nIonic);

  // } else{
  //   std::cout << "Generating proton disordered bulk ice Ih" << std::endl;
  //   Cell sc;
  //   sc = cell.super(scx, scy, scz);
  //   Ice ice(sc);
  //   ice.set_oh_length(oh);
  //   ice.get_h_pos();
  //   ice.get_waters();
  //   ice.get_hbonds();

  //   ice.populate_h_random();
  //   ice.buch_mc_correct();
  //   // ice.write_cell("init.cell");
  //   int nionic = ice.check_ionic_defects();
  //   int nbjerrum = ice.check_bjerrum_defects();
  //   if (nionic != 0){
  //     std::cout << "WARNING: " << nionic << " ionic defects" << std::endl;
  //   }
  //   if (nbjerrum != 0){
  //     std::cout << "WARNING: " << nbjerrum << " Bjerrum defects" << std::endl;
  //   }
  //   if (cq_out) ice.write_cq("iceIh.coord");
  //   else if (vasp_out) ice.write_vasp("iceIh.vasp");
  //   else if (castep_out) ice.write_cell("iceIh.cell");
  // }

  if (slab){ // Construct the slab and/or step
    std::deque<int> slist;
    for (int i=0; i<ice -> s1list.size(); i++) slist.push_back(ice -> s1list[i]);
    for (int i=0; i<ice -> s2list.size(); i++) slist.push_back(ice -> s2list[i]);

    ice -> write_highlight_cell("surface.cell", slist);
    std::deque<int> dOHlist = ice -> find_dOH(2);
    ice -> write_highlight_cell("dOH.cell", dOHlist);
    ice -> order_parameter(surface_nn_cut);

    std::cout << "Building surface-ordered slab" << std::endl;
    ice -> build_ordered_slab(dhkl, 2, 2.0, maxiter);
    ice -> check_defects();
  } else if (step){
    std::cout << "Building slab with step" << std::endl;
    ice -> build_slab(dhkl, 2);
    std::string fname = "ice";
    ice -> build_step(step_direction, slices, step_width, vacuum_gap, 
                      oneside, fname);
  }

  return 0;
}
