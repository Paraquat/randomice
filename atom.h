#include <iomanip>
#include <fstream>
#include <sstream>
#include <Eigen/Core>

class Atom {
  private:
  public:
    std::string name;
    Eigen::Vector3d r;
    int label;

    Atom();
    virtual ~Atom();
    Atom(std::string, int, double, double, double);
    Atom(std::string, int, Eigen::Vector3d);
    Atom(const Atom&);
    Atom& operator= (const Atom&);
    friend std::ostream& operator<< (std::ostream&, Atom&);
    friend std::ofstream& operator<< (std::ofstream&, Atom&);
    bool operator== (const Atom&);
};
