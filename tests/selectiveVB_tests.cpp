#define BOOST_TEST_MODULE "SelectiveVB"

#include "SelectiveVB.hpp"
#include "bits.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>


BOOST_AUTO_TEST_CASE ( Heitler_london ) {

    // Visual test

    vb::State cov{{1,2},{2,1},{1,1}};
    vb::State ion1{{1,2},{1,2},{1,1}};
    libwint::Molecule H2("../../tests/data/h2.xyz");
    libwint::AOBasis ao_basis(H2,"STO-3G");
    ao_basis.calculateIntegrals();
    std::vector<vb::State> states = {cov,ion1};
    vb::SelectiveVB elective_amnesia(ao_basis.get_T()+ao_basis.get_V(),ao_basis.get_g(),ao_basis.get_S(),states);

    std::cout<<elective_amnesia.solve()+H2.calculateInternuclearRepulsionEnergy();
    std::cout<<std::endl;
    std::cout<<elective_amnesia.get_eigenvectors();






}

BOOST_AUTO_TEST_CASE ( Heitler_london_H6 ) {
    // H6 calculation of "all" determinants
    libwint::Molecule H4("../../tests/data/h6.xyz");
    libwint::AOBasis ao_basis(H4,"STO-3G");
    ao_basis.calculateIntegrals();
    std::vector<vb::State> states = {};
    std::vector<size_t> tater{7};
    size_t det=7;
    for(int i =0;i<19;i++){
        det = vb::next_bitset_permutation(det);
        tater.push_back(det);

    }
    for(int i=0;i<20;i++){
        for(int j=0;j<20;j++){
            vb::State dummy{{tater[i]},{tater[j]},{1}};
            states.push_back(dummy);

        }
    }
    std::cout<<std::endl<<states.size()<<std::endl;
    vb::SelectiveVB elective_amnesia(ao_basis.get_T()+ao_basis.get_V(),ao_basis.get_g(),ao_basis.get_S(),states);

    std::cout<<std::endl<<elective_amnesia.solve()+H4.calculateInternuclearRepulsionEnergy()<<std::endl;





}