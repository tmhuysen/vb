#define BOOST_TEST_MODULE "CompleteVB_lih"

#include "CompleteVB.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>


BOOST_AUTO_TEST_CASE ( lih_complete_vb ) {
    // LiH FCI calculation GQCG
    double reference_energy = -7.88253;

    libwint::Molecule LiH("../../tests/data/lihvbref.xyz");
    libwint::AOBasis ao_basis(LiH,"STO-3G");
    ao_basis.calculateIntegrals();
    vb::CompleteVB vb (ao_basis,2,2);
    double internuclear_repulsion = LiH.calculateInternuclearRepulsionEnergy();
    double energy = vb.solve() + internuclear_repulsion;
    BOOST_CHECK(std::abs(energy-reference_energy)<1.0e-6);

}

