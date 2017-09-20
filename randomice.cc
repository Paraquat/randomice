#include "randomice.h"

namespace po = boost::program_options;

int main(int argc, char* argv[]){

  po::options_description desc("Allowed options");
  std::string infname;
  std::vector<int> scdim;
  int scx, scy, scz;

  desc.add_options()
    ("help,h", "Print help information")
    ("infile,i", po::value< std::string >(), "Input structure file")
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

  Cell cell;
  Cell sc;

  cell.read_cell(infname);
  sc = cell.super(scx, scy, scz);
  Ice ice(sc);
  ice.get_h_pos();
  ice.get_waters();
  ice.get_water_nn(oo_max);
  ice.get_hbonds();
  std::cout << ice.nwater << " waters" << std::endl;
  std::cout << ice.nhbond << " hydrogen bonds" << std::endl;
  ice.populate_h_random();
  ice.buch_mc_correct();
  for (int i=0; i<ice.nwater; i++){
    std::cout << ice.waters[i].hbonds.size() << std::endl;
    for (int j=0; j<ice.waters[i].hbonds.size(); j++){
      std::cout << *(ice.waters[i].hbonds[j]);
    }
    std::cout << std::endl;
  }
  // for (int i=0; i<ice.nhbond; i++){
  //   std::cout << ice.hbonds[i] << std::endl;
  // }

  ice.write_cell("test.cell");

  return 0;
}
