#include "cell.h"

// Default construcotr
Cell::Cell() 
{
  natoms = 0;
  frac = true;
}

// Destructor
Cell::~Cell()
{
}

// Initialise the box, no atoms
Cell::Cell(Eigen::Matrix3d l)
{
  lat = l;
  natoms = 0;
  frac = true;
}

// Initialise using matrix of lattice vectors and list of atoms
Cell::Cell(Eigen::Matrix3d l, std::deque<Atom> at)
{
  lat = l;
  atoms = at;
  frac = true;
}

// Copy constructor
Cell::Cell(const Cell& cell)
{
  lat = cell.lat;
  atoms = cell.atoms;
  natoms = cell.natoms;
  frac = cell.frac;
}

// Copy constructor
Cell& Cell::operator= (const Cell& cell)
{
  if (this == &cell) return *this;
  lat = cell.lat;
  atoms = cell.atoms;
  natoms = cell.natoms;
  frac = cell.frac;
  return *this;
}

// Write to ostream
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
    os << cell.atoms[i] << std::endl;
  }
  return os;
}

// Add an atom to the list
void Cell::add_atom(Atom& a)
{
  atoms.push_back(a);
  natoms++;
}

// Read a cell file
void Cell::read_cell(std::string fname)
{
  std::ifstream       infile(fname.c_str());
  std::string         line;
  std::string         species;
  double              x, y, z;
  int                 label;
  bool                readatoms = false;

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
    lat_inv = lat.inverse();

    label = 1;
    if (linestr == "%block positions_frac"){
      frac = true;
      readatoms = true;
    } else if (linestr == "%block positions_abs") {
      frac = false;
      readatoms = true;
    }
    if (readatoms){
      while (std::getline(infile, line)){
        std::istringstream  iss(line);
        linestr = iss.str();
        if (linestr.at(0) == '%'){
          std::transform(linestr.begin(), linestr.end(), linestr.begin(), \
                         ::tolower);
        }
        if (linestr == "%endblock positions_frac") break;
        if (linestr == "%endblock positions_abs") break;
        iss >> species >> x >> y >> z;
        Atom a(species, label, x, y, z);
        add_atom(a);
        label++;
      }
    }
    if (frac == false) {
      cart2frac_all();
    }
  }
}

// Write object ot a .cell file
void Cell::write_cell(std::string fname)
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

  for (int i=0; i<natoms; i++) {
    if (atoms[i].occupied){
      outfile << atoms[i] << std::endl;
    }
  }

  if (frac) outfile << "%ENDBLOCK positions_frac" << std::endl << std::endl;
  else outfile << "%ENDBLOCK positions_abs" << std::endl << std::endl;
}

// Wrap the atoms back into the unit cell if necessary
void Cell::wrap(void)
{
  assert(frac == true);
  for (int i=0; i<natoms; i++) {
    for (int j=0; j<=2; j++) {
      while (1) {
        if (atoms[i].r(j) < 1.0 && atoms[i].r(j) >= 0.0) break;
        if (atoms[i].r(j) > 1.0) atoms[i].r(j) -= 1.0;
        if (atoms[i].r(j) < 0.0) atoms[i].r(j) += 1.0;
      }
    }
  }
}

// Convert one atomic coordinate from fractional to Cartesian
void Cell::frac2cart(Atom& a)
{
  a.r = lat*a.r;
}

// Convert all coordinates from fractional to Cartesian
void Cell::frac2cart_all(void)
{
  wrap();
  for (int i=0; i<natoms; i++) {
    frac2cart(atoms[i]);
  }  
  frac = false;
}

// Convert one atomic coordinate from Cartesian to fractional
void Cell::cart2frac(Atom& a)
{
  a.r = lat_inv*a.r ;
}

// Convert all atomic coordinates from Cartesian to fractional
void Cell::cart2frac_all(void)
{
  for (int i=0; i<natoms; i++) {
    cart2frac(atoms[i]);
  }  
  frac = true;
  wrap();
}

// Compute the vector between two atoms with PBCs using the Minimum Image
// Convention (MIC) in Cartesian basis --- only works with an orthorhombic cell
Eigen::Vector3d Cell::mic_cart(Atom a, Atom b)
{
  Eigen::Vector3d d;

  d = b.r - a.r;
  for (int i=0; i<=2; i++) {
    d(i) = d(i) - round(d(i)/lat(i,i))*lat(i,i);
  }
  return d;
}

// Compute the vector between two atoms with PBCs using the Minimum Image
// Convention in fractional coordinates
Eigen::Vector3d Cell::mic_frac(Atom a, Atom b)
{
  Eigen::Vector3d d;

  d = b.r - a.r;
  for (int i=0; i<=2; i++) {
    if (d(i) > 0.5) d(i) -= 1.0;
    else if (d(i) < -0.5) d(i) += 1.0;
  }
  return lat*d;
}

// Compute the distance table
void  Cell::get_dt(void)
{
  Eigen::Vector3d d;
  dt.setZero(natoms,natoms);
  for (int i=0; i<natoms; i++) {
    for (int j=i+1; j<natoms; j++) {
      if (frac) d = mic_frac(atoms[i], atoms[j]);
      else d = mic_cart(atoms[i], atoms[j]);
      dt(i,j) = d.norm();
      dt(j,i) = dt(i,j);
    }
  }
}

// Generate a axbxc supercell
Cell Cell::super(int a, int b, int c)
{
  Eigen::Vector3d t, rt, scdim;
  Eigen::Matrix3d lat_super;
  std::vector<Atom> atoms_super;

  assert(frac == true);
  lat_super.row(0) = lat.row(0)*static_cast<double>(a);
  lat_super.row(1) = lat.row(1)*static_cast<double>(b);
  lat_super.row(2) = lat.row(2)*static_cast<double>(c);
  Cell sc(lat_super);

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
          sc.add_atom(a);
          label += 1;
        }
      }
    }
  }
  return sc;
}

// Find nearest neighbours within cutoff
void Cell::get_nn(double cutoff)
{
  get_dt();
  for (int i=0; i<natoms; i++){
    atoms[i].nn.clear();
    for (int j=0; j<natoms; j++){
      if (i != j){
        if (dt(i,j) < cutoff){
          Atom::atom_ptr ptr(new Atom(atoms[j]));
          atoms[i].add_nn(ptr);
        }
      }
    }
  } 
}

// Check whether atom c lines on a line connecting atoms a and b
bool Cell::isPointOnLine(Atom& a, Atom& b, Atom& c)
{
  Eigen::Vector3d v1, v2, cp;
  double dp;
  bool on_line = false;

  assert(frac == false);
  v1 = mic_cart(b, a);
  v2 = mic_cart(c, a);
  cp = v1.cross(v2);
  dp = v1.dot(v2);
  // std::cout << "cp: " << cp.norm() << " dp: " << dp << " |v1|^2: " << v1.squaredNorm() << std::endl;
  if (cp.norm() < small){
    if (dp > 0.0){
      if (dp < v1.squaredNorm()){
        on_line = true;
      }
    }
  }
  return on_line;
}