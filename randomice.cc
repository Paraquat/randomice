#include "randomice.h"

namespace po = boost::program_options;

int main(int argc, char* argv[]){

  po::options_description desc("Allowed options");
  std::string infname;
  std::vector<int> scdim;
  int scx, scy, scz;

  using namespace backward;
  StackTrace st; st.load_here(32);
  Printer p;
  p.object = true;
  p.color_mode = ColorMode::always;
  p.address = true;
  p.print(st, stderr);

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
  ice.get_hbonds();

  std::cout << ice.nwater << " waters" << std::endl;
  std::cout << ice.nhbond << " hydrogen bonds" << std::endl;
  ice.populate_h_random();
  ice.buch_mc_correct();
  int nionic = ice.check_ionic_defects();
  int nbjerrum = ice.check_bjerrum_defects();
  if (nionic != 0){
    std::cout << "WARNING: " << nionic << " ionic defects" << std::endl;
  }
  if (nbjerrum != 0){
    std::cout << "WARNING: " << nbjerrum << " Bjerrum defects" << std::endl;
  }
  ice.rick_randomise(1000);
  std::cout << "Randomisation complete";
  nionic = ice.check_ionic_defects();
  nbjerrum = ice.check_bjerrum_defects();
  if (nionic != 0){
    std::cout << "WARNING: " << nionic << " ionic defects" << std::endl;
  }
  if (nbjerrum != 0){
    std::cout << "WARNING: " << nbjerrum << " Bjerrum defects" << std::endl;
  }
  ice.write_cell("iceIh.cell");

  return 0;
}
