#include "common.hpp"
namespace vb {
size_t bitstring(std::vector<int> indexes)
{
    size_t result = 0;
    for(int x:indexes){
        result += 1<<x;
    }
    return result;
};

}