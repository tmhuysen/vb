/*
 * Hard-code exe cause i'm lazy
 */

#include <libwint.hpp>
#include <iomanip>
#include "SelectiveVB.hpp"


int main() {
    std::string filenamebase = "../../exe/H2_data/H2_";
    std::ofstream outfile("H2_data");
    std::ofstream outfile2("H2_coefs");
    outfile<<"Bohr\tHartree"<<std::endl;
    outfile<<std::setprecision(10);
    outfile2<<std::setprecision(10);
    for(double i = 0.5;1<10.05;i+= 0.02) {
        std::stringstream stream;
        stream << std::fixed << std::setprecision(2) << i;
        std::string s = stream.str();
        std::string filename = filenamebase + s + ".xyz";
        libwint::Molecule H2(filename);
        vb::State cov1{{1,2}, {2,1}};
        vb::State ionHe{{1}, {1}};
        vb::State ionH{{2}, {2}};
        libwint::AOBasis ao_basis(H2, "STO-3G");
        ao_basis.calculateIntegrals();
        std::vector<vb::State> states = {cov1,ionHe,ionH};
        vb::SelectiveVB elective_amnesia(ao_basis.get_T() + ao_basis.get_V(), ao_basis.get_g(), ao_basis.get_S(), states);
        double repulsion = H2.calculateInternuclearRepulsionEnergy();
        outfile<<i<<"\t"<<elective_amnesia.solve()+repulsion<<std::endl;
        outfile2<<i<<" : "<<std::endl;
        outfile2<<elective_amnesia.get_eigenvectors()<<std::endl;

    }

}