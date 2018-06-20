#ifndef VB_COMPLETEVB_HPP
#define VB_COMPLETEVB_HPP

#include "libwint.hpp"
#include "bmqc.hpp"
#include <boost/numeric/conversion/converter.hpp>
#include <boost/math/special_functions.hpp>


namespace vb {

class CompleteVB {
public:
    // Public Member Variables : Developmental Stage TODO: Change
    libwint::AOBasis ao_basis;
    const size_t dim;
    const size_t dim_alpha;  // the dimension of the alpha CI space
    const size_t dim_beta;  // the dimension of the beta CI space

    const size_t K;  // number of spatial orbitals

    const size_t N_alpha;  // number of alpha electrons
    const size_t N_beta;  // number of beta electrons

    const bmqc::AddressingScheme addressing_scheme_alpha;
    const bmqc::AddressingScheme addressing_scheme_beta;

    Eigen::MatrixXd oei;
    Eigen::Tensor<double,4> tei;
    Eigen::MatrixXd oi;

    Eigen::MatrixXd hamiltonian_alpha;
    Eigen::MatrixXd hamiltonian_beta;
    Eigen::MatrixXd hamiltonian;

    Eigen::MatrixXd overlap_alpha;
    Eigen::MatrixXd overlap_beta;
    Eigen::MatrixXd overlap;

    struct OneElectronCoupling {
        int sign;
        size_t p;
        size_t q;
        double overlap;
        size_t address_target;
    };

    std::vector<size_t> as;
    std::vector<size_t> bs;
    std::vector<std::vector<OneElectronCoupling>> aoec;
    std::vector<std::vector<OneElectronCoupling>> boec;



    void calculate_separated_elements(size_t one_string, size_t two_string,size_t index_one, size_t index_two, Eigen::MatrixXd& overlap, Eigen::MatrixXd& hamiltonian, std::vector<std::vector<OneElectronCoupling>>& coupling);
    void calculate_mixed_elements(size_t aone_string, size_t atwo_string, size_t bone_string, size_t btwo_string,
                                  size_t index_one, size_t index_two);
    void mixer();
    double calculate_overlap(size_t aone_string, size_t atwo_string);

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

    /*
     *  Recombines alpha and beta separated calculations
     */
    void finish_off();

    /*
     * Solves
     */
    double solve();
};

}  // namespace vb
#endif //VB_COMPLETEVB_HPP
