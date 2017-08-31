#include "randomice.h"

namespace po = boost::program_options;

int main(int argc, char* argv[]){

  po::options_description desc("Allowed options");
  std::string infname;

  desc.add_options()
    ("help,h", "Print help information")
    ("infile,i", po::value< std::string >(), "Input structure file")
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

  Cell cell;
  Cell sc;

  cell.read_cell(infname);
  sc = cell.super(2, 2, 2);
  std::cout << sc;

  return 0;
}
