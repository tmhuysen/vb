#include "CompleteVB.hpp"

namespace vb{

CompleteVB::CompleteVB(libwint::AOBasis& ao_basis, size_t N_A, size_t N_B):
        ao_basis(ao_basis), N_alpha(N_A), N_beta(N_B),K(ao_basis.calculateNumberOfBasisFunctions()),
        addressing_scheme_alpha (bmqc::AddressingScheme(this->K, this->N_alpha)),
        addressing_scheme_beta (bmqc::AddressingScheme(this->K, this->N_beta)),
        dim_alpha (CompleteVB::calculateDimension(K, N_alpha, 0)),
        dim_beta (CompleteVB::calculateDimension(K, 0, N_beta)),
        dim(dim_alpha*dim_beta) {
    this->oei = ao_basis.get_T()+ao_basis.get_V();
    this->oi = ao_basis.get_S();
    this->tei = ao_basis.get_g();
}

void CompleteVB::calculate_matrices_alpha() {

    bmqc::SpinString<unsigned long> ssa_one (0, this->addressing_scheme_alpha);  // alpha spin string with address 0
    for (size_t I_alpha = 0; I_alpha < this->dim_alpha; I_alpha++) {
        bmqc::SpinString<unsigned long> ssa_two (0, this->addressing_scheme_alpha);  // alpha spin string with address 0
        for (size_t J_alpha = 0; J_alpha < this->dim_alpha; J_alpha++){
            calculate_element(ssa_one.get_representation(),ssa_two.get_representation(),I_alpha,J_alpha);



        }
    }



}

size_t CompleteVB::calculateDimension(size_t K, size_t N_alpha, size_t N_beta) {
    // K, N_alpha, N_beta are expected to be small, so static-casting them to unsigned (what boost needs) is permitted.
    auto dim_double_alpha = boost::math::binomial_coefficient<double>(static_cast<unsigned>(K), static_cast<unsigned>(N_alpha));
    auto dim_double_beta = boost::math::binomial_coefficient<double>(static_cast<unsigned>(K), static_cast<unsigned>(N_beta));
    auto dim_double = dim_double_alpha * dim_double_beta;

    // Check if the resulting dimension is appropriate to be stored in size_t
    return boost::numeric::converter<double, size_t>::convert(dim_double);
}

void CompleteVB::calculate_element(size_t one_string, size_t two_string, size_t index_one, size_t index_two) {




}


}  // namespace vb