#ifndef VB_COMMON_HPP
#define VB_COMMON_HPP

#include <cstddef>
#include <vector>
namespace vb {

struct OneElectronCoupling {
    int sign;
    size_t p;
    size_t q;
    double overlap;
    size_t address_target;
};


struct State {

    std::vector<size_t> alpha;
    std::vector<size_t> beta;
    //std::vector<double> weights;

};

size_t bitstring(std::vector<int> indexes);
size_t bitstring(size_t);

}
#endif //VB_COMMON_HPP
