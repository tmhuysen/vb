/*
 * Hard-code exe cause i'm lazy
 */

#include <libwint.hpp>
#include <iomanip>
#include <hf.hpp>
#include "SelectiveVB.hpp"


int main() {
    std::string filenamebase = "../../exe/H2_data/H2_";
    std::ofstream outfile("H2_datarhf");
    outfile<<"Bohr\tHartree"<<std::endl;
    outfile<<std::setprecision(10);
    for(double i = 0.5;1<10.05;i+= 0.02) {
        std::stringstream stream;
        stream << std::fixed << std::setprecision(2) << i;
        std::string s = stream.str();
        std::string filename = filenamebase + s + ".xyz";
        libwint::Molecule H2(filename);

        libwint::AOBasis ao_basis(H2, "STO-3G");
        ao_basis.calculateIntegrals();
        hf::rhf::RHF elective_amnesia (H2, ao_basis, 1.0e-12, 100000);
        elective_amnesia.solve( hf::rhf::solver::SCFSolverType::DIIS);
        double repulsion = H2.calculateInternuclearRepulsionEnergy();
        outfile<<i<<"\t"<<elective_amnesia.get_electronic_energy()+repulsion<<std::endl;

    }

}
