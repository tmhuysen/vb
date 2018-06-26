 * Hard-code exe cause i'm lazy
 */

#include <libwint.hpp>
#include "CompleteVB.hpp"
#include "SelectiveVB.hpp"
#include "bits.hpp"


int main() {
    libwint::Molecule H4("../../tests/data/h6.xyz");
    libwint::AOBasis ao_basis(H4,"STO-3G");
    ao_basis.calculateIntegrals();
    std::vector<vb::State> states = {};
    std::vector<size_t> tater{7};
    size_t sven=7;
    for(int i =0;i<19;i++){
        sven = vb::next_bitset_permutation(sven);
        tater.push_back(sven);

    }
    for(int i=0;i<20;i++){
        for(int j=0;j<20;j++){
            vb::State dummy{{tater[i]},{tater[j]}};
            states.push_back(dummy);

        }
    }
    std::cout<<std::endl<<states.size()<<std::endl;
    vb::SelectiveVB elective_amnesia(ao_basis.get_T()+ao_basis.get_V(),ao_basis.get_g(),ao_basis.get_S(),states);

    std::cout<<std::endl<<elective_amnesia.solve()+H4.calculateInternuclearRepulsionEnergy()<<std::endl;

}
