#define BOOST_TEST_MODULE "SelectiveVB"

#include "SelectiveVB.hpp"
#include "bits.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>


BOOST_AUTO_TEST_CASE ( Heitler_london ) {
    // Cristina's H2 FCI energy/OO-DOCI energy
    double reference_energy = -1.13726;
    vb::State cov{{1,2},{2,1},{1,1}};
    vb::State ion1{{1,2},{1,2},{1,1}};
    libwint::Molecule H2("../../tests/data/h2.xyz");
    libwint::AOBasis ao_basis(H2,"STO-3G");
    ao_basis.calculateIntegrals();
    std::vector<vb::State> states = {cov,ion1};
    vb::SelectiveVB elective_amnesia(ao_basis.get_T()+ao_basis.get_V(),ao_basis.get_g(),ao_basis.get_S(),states);
    double energy = elective_amnesia.solve()+H2.calculateInternuclearRepulsionEnergy();
    BOOST_CHECK(std::abs(energy-reference_energy)<1.0e-6);

}

