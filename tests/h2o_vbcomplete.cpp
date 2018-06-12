#define BOOST_TEST_MODULE "CompleteVB_h2o"

#include "CompleteVB.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>


BOOST_AUTO_TEST_CASE ( H2O_complete_vb ) {
    // LiH FCI calculation GQCG
    double reference_energy = -7.88253;

    libwint::Molecule LiH("../../tests/data/h2o.xyz");
    libwint::AOBasis ao_basis(LiH,"STO-3G");
    ao_basis.calculateIntegrals();
    vb::CompleteVB vb (ao_basis,5,5);
    std::cout<<vb.solve()+LiH.calculateInternuclearRepulsionEnergy();

}

