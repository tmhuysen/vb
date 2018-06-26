#define BOOST_TEST_MODULE "CompleteVB_h2o"

#include "CompleteVB.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>


BOOST_AUTO_TEST_CASE ( H2O_complete_vb ) {
    // Psi4 and GAMESS' FCI energy
    double reference_energy = -75.0129803939602;

    libwint::Molecule H2O("../../tests/data/h2o.xyz");
    libwint::AOBasis ao_basis(H2O,"STO-3G");
    ao_basis.calculateIntegrals();
    vb::CompleteVB vb (ao_basis,5,5);

    double internuclear_repulsion = H2O.calculateInternuclearRepulsionEnergy();
    double energy = vb.solve() + internuclear_repulsion;

    BOOST_CHECK(std::abs(energy-reference_energy)<1.0e-6);

}


BOOST_AUTO_TEST_CASE ( H4_ref_visual ) {
    // generate reference H4

    libwint::Molecule H4("../../tests/data/h4.xyz");
    libwint::AOBasis ao_basis(H4,"STO-3G");
    ao_basis.calculateIntegrals();
    vb::CompleteVB vb (ao_basis,2,2);

    double internuclear_repulsion = H4.calculateInternuclearRepulsionEnergy();
    double energy = vb.solve() + internuclear_repulsion;
    std::cout<<std::endl<<energy;

}