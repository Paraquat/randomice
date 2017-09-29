#include "randomice.h"

namespace po = boost::program_options;

int main(int argc, char* argv[]){

  po::options_description desc("Allowed options");
  std::string infname;
  std::vector<int> scdim;
  int scx, scy, scz;
  int maxiter = 1000;

  desc.add_options()
    ("help,h", "Print help information")
    ("infile,i", po::value< std::string >(), "Input structure file")
    ("debug,d", po::value< bool >(), "Run in debug mode")
    ("maxiter,m", po::value< int >(),
     "Maximum number of Rick algorithm iterations")
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

  Cell cell;
  Cell sc;

  cell.read_cell(infname);
  sc = cell.super(scx, scy, scz);
  Ice ice(sc);
  ice.get_h_pos();
  ice.get_waters();
  ice.get_hbonds();

  std::cout << ice.nwater << " waters" << std::endl;
  std::cout << ice.nhbond << " hydrogen bonds" << std::endl;
  ice.populate_h_random();
  ice.buch_mc_correct();
  ice.write_cell("init.cell");
  int nionic = ice.check_ionic_defects();
  int nbjerrum = ice.check_bjerrum_defects();
  if (nionic != 0){
    std::cout << "WARNING: " << nionic << " ionic defects" << std::endl;
  }
  if (nbjerrum != 0){
    std::cout << "WARNING: " << nbjerrum << " Bjerrum defects" << std::endl;
  }
  // ice.rick_randomise(100000);
  // std::cout << "Randomisation complete";
  // ice.build_slab(dhkl_default, 2);
  ice.build_ordered_slab(dhkl_default, 2, 2.0, maxiter);
  nionic = ice.check_ionic_defects();
  nbjerrum = ice.check_bjerrum_defects();
  if (nionic != 0){
    std::cout << "WARNING: " << nionic << " ionic defects" << std::endl;
  }
  if (nbjerrum != 0){
    std::cout << "WARNING: " << nbjerrum << " Bjerrum defects" << std::endl;
  }
  std::deque<int> slist;
  for (int i=0; i<ice.s1list.size(); i++) slist.push_back(ice.s1list[i]);
  for (int i=0; i<ice.s2list.size(); i++) slist.push_back(ice.s2list[i]);

  ice.write_highlight_cell("surface.cell", slist);
  std::deque<int> dOHlist = ice.find_dOH(2);
  ice.write_highlight_cell("dOH.cell", dOHlist);
  ice.order_parameter(surface_nn_cut);

  ice.write_cell("iceIh.cell");

  return 0;
}
