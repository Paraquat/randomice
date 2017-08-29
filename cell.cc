#include "cell.h"

Cell::Cell()
{
  natoms = 0;
}

Cell::~Cell()
{
}

Cell::Cell(Eigen::Matrix3f l, std::vector<Atom> at)
{
  lat = l;
  atoms = at;
}

Cell::Cell(const Cell& cell)
{
  lat = cell.lat;
  atoms = cell.atoms;
  natoms = cell.natoms;
}

std::ostream& operator<< (std::ostream& os, Cell& cell)
{
  os << "Number of atoms: " << cell.natoms << std::endl;
  os << "Lattice vectors:" << std::endl;
  for (int i=0; i<=2; i++){
    os << std::setprecision(8) << std::setw(12) \
       << cell.lat.row(i) << std::endl;
  }
  os << "Atomic positions:" << std::endl;
  for (int i=0; i<cell.natoms; i++){
    os << cell.atoms[i];
  }
  return os;
}

void Cell::add_atom(Atom a)
{
  atoms.push_back(a);
  natoms++;
}

void Cell::read_cell(std::string fname)
{
  std::ifstream       infile(fname.c_str());;
  std::string         line;
  std::string         species;
  double              x, y, z;
  int                 label;

  if (!infile){
    std::cout << "Error opening file " << fname << std::endl;
    std::exit(1);
  }

  while (std::getline(infile,  line))
  {
    std::istringstream  iss(line);
    std::string linestr = iss.str();
    
    if (linestr.empty()) continue;
    if (linestr.at(0) == '%'){
      std::transform(linestr.begin(), linestr.end(), linestr.begin(), \
                     ::tolower);
    }
    if (linestr == "%block lattice_cart"){
      for (int i=0; i<=2; i++){
        std::getline(infile, line);
        std::istringstream  iss(line);
        iss >> lat(i,0) >> lat(i,1) >> lat(i,2);
      }
    }

    label = 1;
    if (linestr == "%block positions_frac"){
      while (std::getline(infile, line)){
        std::istringstream  iss(line);
        linestr = iss.str();
        std::transform(linestr.begin(), linestr.end(), linestr.begin(), \
                       ::tolower);
        if (linestr == "%endblock positions_frac") break;
        iss >> species >> x >> y >> z;
        add_atom(Atom(species, label, x, y, z));
        label++;
      }
    }
  }
}