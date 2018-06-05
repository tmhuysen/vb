#ifndef VB_COMPLETEVB_HPP
#define VB_COMPLETEVB_HPP

#include "libwint.hpp"
#include "bmqc.hpp"
#include <boost/numeric/conversion/converter.hpp>
#include <boost/math/special_functions.hpp>


namespace vb {

class CompleteVB {
public:
    // Public Member Variables : Developmental Stage (TO:DO Change)
    libwint::AOBasis ao_basis;
    const size_t dim;
    const size_t dim_alpha;  // the dimension of the alpha CI space
    const size_t dim_beta;  // the dimension of the beta CI space

    const size_t K;  // number of spatial orbitals

    const size_t N_alpha;  // number of alpha electrons
    const size_t N_beta;  // number of beta electrons

    const bmqc::AddressingScheme addressing_scheme_alpha;
    const bmqc::AddressingScheme addressing_scheme_beta;

    Eigen::MatrixXd hamiltonian;
    Eigen::MatrixXd overlap;


    // CONSTRUCTORS
    /**
     *  Constructor based on a given @param ao_basis, a number of alpha electrons @param N_alpha and a number of beta electrons
     *  @param N_beta.
     */
    CompleteVB(libwint::AOBasis &ao_basis, size_t N_A, size_t N_B);


    // STATIC PUBLIC METHODS
    /**
     *  Given a number of spatial orbitals @param K, a number of alpha electrons @param N_A, and a number of beta electrons
     *  @param N_B, @return the dimension of the FCI space.
     */
    static size_t calculateDimension(size_t K, size_t N_alpha, size_t N_beta);


    // PUBLIC METHODS
    /**
     *  Calculates the overlap Matrix (overlap) and the (non-orthogonalized) Hamiltonian (hamiltonian)
     */
    void calculate_matrices();
};

}  // namespace vb
#endif //VB_COMPLETEVB_HPP
