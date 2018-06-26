/*
 * Hard-code exe cause i'm lazy
 */

#include <libwint.hpp>
#include <hf.hpp>
#include "CompleteVB.hpp"


int main() {
    std::string filename = "../../exe/h2.xyz";
    libwint::Molecule h2 (filename);
    libwint::AOBasis ao_basis (h2, "STO-3G");
    double repulsion = h2.calculateInternuclearRepulsionEnergy();
    ao_basis.calculateIntegrals();
    vb::CompleteVB vb(ao_basis,1,1);
    std::cout<<std::endl;
    std::cout<<std::endl;
    std::cout<<vb.solve()+repulsion;




}
