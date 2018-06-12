#define BOOST_TEST_MODULE "CompleteVB_h2"

#include "CompleteVB.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>


BOOST_AUTO_TEST_CASE ( h2_complete_vb ) {
    // LiH FCI calculation GQCG
    double reference_energy = -7.88253;

    libwint::Molecule H2("../../tests/data/h2.xyz");
    libwint::AOBasis ao_basis(H2,"STO-3G");
    ao_basis.calculateIntegrals();
    vb::CompleteVB vb (ao_basis,1,1);
    std::cout<<vb.solve()+H2.calculateInternuclearRepulsionEnergy();
    std::cout<<std::endl<<std::endl<<vb.hamiltonian<<std::endl<<std::endl;
    std::cout<<vb.overlap<<std::endl<<std::endl;
    std::cout<<vb.overlap_alpha<<std::endl<<std::endl;
    std::cout<<vb.overlap_beta<<std::endl<<std::endl;

}

