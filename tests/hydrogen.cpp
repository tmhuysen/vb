#define BOOST_TEST_MODULE "CompleteVB_h2"

#include "CompleteVB.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>


BOOST_AUTO_TEST_CASE ( h2_complete_vb ) {
    // Cristina's H2 FCI energy/OO-DOCI energy
    double reference_energy = -1.1651486697;

    libwint::Molecule H2("../../tests/data/h2_cristina.xyz");
    libwint::AOBasis ao_basis(H2,"6-31g**");
    ao_basis.calculateIntegrals();
    vb::CompleteVB vb (ao_basis,1,1);

    double internuclear_repulsion = H2.calculateInternuclearRepulsionEnergy();
    double energy = vb.solve() + internuclear_repulsion;

    BOOST_CHECK(std::abs(energy-reference_energy)<1.0e-6);

}

