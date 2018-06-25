/*
 * Hard-code exe cause i'm lazy
 */

#include <libwint.hpp>
#include <iomanip>
#include "SelectiveVB.hpp"
#include <hf.hpp>

int main() {
    std::string filenamebase = "../../exe/HeH_data/HeH_";
    std::ofstream outfile("HeH_datarhf");
    outfile<<"Bohr\tHartree"<<std::endl;
    outfile<<std::setprecision(10);
    for(double i = 0.5;1<10.05;i+= 0.02) {
        std::stringstream stream;
        stream << std::fixed << std::setprecision(2) << i;
        std::string s = stream.str();
        std::string filename = filenamebase + s + ".xyz";
        libwint::Molecule HeH(filename,1);


        libwint::AOBasis ao_basis(HeH, "STO-3G");
        ao_basis.calculateIntegrals();
        hf::rhf::RHF elective_amnesia (HeH, ao_basis, 1.0e-12, 100000);
        elective_amnesia.solve( hf::rhf::solver::SCFSolverType::PLAIN);
        double repulsion = HeH.calculateInternuclearRepulsionEnergy();
        outfile<<i<<"\t"<<elective_amnesia.get_electronic_energy()+repulsion<<std::endl;
    }

}
