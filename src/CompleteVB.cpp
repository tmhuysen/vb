#include "CompleteVB.hpp"

namespace vb{

CompleteVB::CompleteVB(libwint::AOBasis& ao_basis, size_t N_A, size_t N_B):
        ao_basis(ao_basis), N_alpha(N_A), N_beta(N_B),K(ao_basis.calculateNumberOfBasisFunctions()),
        addressing_scheme_alpha (bmqc::AddressingScheme(this->K, this->N_alpha)),
        addressing_scheme_beta (bmqc::AddressingScheme(this->K, this->N_beta)),
        dim_alpha (CompleteVB::calculateDimension(K, N_alpha, 0)),
        dim_beta (CompleteVB::calculateDimension(K, 0, N_beta)),
        dim(dim_alpha*dim_beta) {}

void CompleteVB::calculate_matrices() {

}

size_t CompleteVB::calculateDimension(size_t K, size_t N_alpha, size_t N_beta) {
    // K, N_alpha, N_beta are expected to be small, so static-casting them to unsigned (what boost needs) is permitted.
    auto dim_double_alpha = boost::math::binomial_coefficient<double>(static_cast<unsigned>(K), static_cast<unsigned>(N_alpha));
    auto dim_double_beta = boost::math::binomial_coefficient<double>(static_cast<unsigned>(K), static_cast<unsigned>(N_beta));
    auto dim_double = dim_double_alpha * dim_double_beta;

    // Check if the resulting dimension is appropriate to be stored in size_t
    return boost::numeric::converter<double, size_t>::convert(dim_double);
}


}  // namespace vb