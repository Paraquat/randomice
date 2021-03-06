#include "randomice.h"

namespace po = boost::program_options;

int main(int argc, char* argv[]){

  po::options_description desc("Allowed options");
  std::string infname, phase;
  std::vector<int> scdim;
  int scx, scy, scz;
  int maxiter = 1000;
  bool cq_out = false;
  bool ordered = false;
  bool oneside = false;
  bool castep_out = false;
  bool vasp_out = false;
  bool flag_debug = false;
  bool step = false;

  desc.add_options()
    ("help,h", "Print help information")
    ("infile,i", po::value< std::string >(), "Input structure file")
    ("phase,p", po::value< std::string >(), "Input ice phase (Ih, Ic, etc")
    ("debug,d", po::bool_switch(&flag_debug), "Run in debug mode")
    ("maxiter,m", po::value< int >(),
     "Maximum number of Rick algorithm iterations")
    ("cq", po::bool_switch(&cq_out), "Write to Conquest file")
    ("ordered,o", po::bool_switch(&ordered), "Generate an ordered supercell")
    ("oneside", po::bool_switch(&oneside), "Build step on only one side of slab")
    ("cell", po::bool_switch(&castep_out), "Write to CASTEP cell file")
    ("vasp", po::bool_switch(&vasp_out), "Write to VASP POSCAR file")
    ("step", po::bool_switch(&step), "Construct a step")
    ("step_direction", po::value< std::string >(), "Direction of step (a,b,c)")
    ("step_width", po::value< double >(), "Width of step (a0)")
    ("vacuum_gap", po::value< double >(), "Size of vacuum gap (a0)")
    ("supercell,s", po::value<std::vector<int> >() -> multitoken(), 
     "Supercell dimensions (a x b x c)")
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
    std::cout << "Input file is " <<  infname << std::endl;
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

  bool build_step = false;
  double step_width = 0.0;
  double vacuum_gap = 0.0;
  std::string step_direction = "a";
  if (step){
    build_step = true;
    if (vm["step_width"].empty()){
      step_width = 0.0;
    } else {
      step_width = vm["step_width"].as<double>();
    }
    vacuum_gap = vm["vacuum_gap"].as<double>();
    step_direction = vm["step_direction"].as<std::string>();
  }

  Cell cell;

  cell.read_cell(infname);
  double dhkl = cell.lat(2,2)/2.0;
  Ice ice;
  if (ordered){
    Ice unit;
    unit.read_h_pos(cell);
    ice = unit.super(scx, scy, scz);
    ice.build_slab(dhkl, 2);

  } else{
    Cell sc;
    sc = cell.super(scx, scy, scz);
    Ice ice(sc);
    if (flag_debug) ice.flag_debug = true;
    else ice.flag_debug = false;
    ice.get_h_pos();
    ice.get_waters();
    ice.get_hbonds();

    std::cout << ice.nwater << " waters" << std::endl;
    std::cout << ice.nhbond << " hydrogen bonds" << std::endl;
    ice.populate_h_random();
    ice.buch_mc_correct();
    // ice.write_cell("init.cell");
    int nionic = ice.check_ionic_defects();
    int nbjerrum = ice.check_bjerrum_defects();
    if (nionic != 0){
      std::cout << "WARNING: " << nionic << " ionic defects" << std::endl;
    }
    if (nbjerrum != 0){
      std::cout << "WARNING: " << nbjerrum << " Bjerrum defects" << std::endl;
    }
    ice.build_ordered_slab(dhkl, 2, 2.0, maxiter);
    nionic = ice.check_ionic_defects();
    nbjerrum = ice.check_bjerrum_defects();
    if (nionic != 0){
      std::cout << "WARNING: " << nionic << " ionic defects" << std::endl;
    }
    if (nbjerrum != 0){
      std::cout << "WARNING: " << nbjerrum << " Bjerrum defects" << std::endl;
    }
  }

  if (build_step) {
    if (phase == "Ih"){
      ice.xlayers = 2*scx;
      ice.ylayers = 2*scy;    
    }
    std::string fname = "ice";
    ice.build_step(step_direction, step_width, vacuum_gap, oneside, fname);
  } else{
    std::deque<int> slist;
    for (int i=0; i<ice.s1list.size(); i++) slist.push_back(ice.s1list[i]);
    for (int i=0; i<ice.s2list.size(); i++) slist.push_back(ice.s2list[i]);

    ice.write_highlight_cell("surface.cell", slist);
    std::deque<int> dOHlist = ice.find_dOH(2);
    ice.write_highlight_cell("dOH.cell", dOHlist);
    ice.order_parameter(surface_nn_cut);

    if (cq_out) ice.write_cq("iceIh.coord");
    else if (vasp_out) ice.write_vasp("iceIh.vasp");
    else ice.write_cell("iceIh.cell");
  }

  return 0;
}
