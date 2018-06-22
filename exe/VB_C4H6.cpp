/*
 * Hard-code exe cause i'm lazy
 */

#include <libwint.hpp>
#include "CompleteVB.hpp"


int main() {
    std::string filename = "../exe/h2o.xyz";
    libwint::Molecule LiH (filename);
    libwint::AOBasis ao_basis (LiH, "STO-3G");
    double repulsion = LiH.calculateInternuclearRepulsionEnergy();
    ao_basis.calculateIntegrals();
    vb::CompleteVB fci(ao_basis,5,5);
    std::cout<<fci.solve()+repulsion;
}
