#ifndef VB_COMPLETEVB_HPP
#define VB_COMPLETEVB_HPP

#include "libwint.hpp"
#include "common.hpp"
#include <boost/numeric/conversion/converter.hpp>
#include <boost/math/special_functions.hpp>


namespace vb {




class CompleteVB {
private:
    libwint::AOBasis ao_basis;  // The atomic basis-set
    const size_t dim;
    const size_t dim_alpha;  // the dimension of the alpha CI space
    const size_t dim_beta;  // the dimension of the beta CI space

    const size_t K;  // number of spatial orbitals

    const size_t N_alpha;  // number of alpha electrons
    const size_t N_beta;  // number of beta electrons

    Eigen::MatrixXd oei;  // one electeron integrals
    Eigen::Tensor<double,4> tei;  // two electron integrals
    Eigen::MatrixXd oi;  // overlap integrals

    Eigen::MatrixXd hamiltonian_alpha; // Alpha exclusive values (4-alpha two electron integrals, 2-alpha one electron integrals)
    Eigen::MatrixXd hamiltonian_beta;  // Beta exclusive values
    Eigen::MatrixXd hamiltonian;  // Combinations of all alpha and beta exclusive values + 2 * (alpha beta beta alpha two electron integrals)

    Eigen::MatrixXd overlap_alpha;
    Eigen::MatrixXd overlap_beta;
    Eigen::MatrixXd overlap;

    std::vector<size_t> as;  // Alpha bitstrings
    std::vector<size_t> bs;  // Beta bitstrings

    std::vector<std::vector<OneElectronCoupling>> aoec;  // overlap after a one electron evaluation
    std::vector<std::vector<OneElectronCoupling>> boec;

    Eigen::MatrixXd eigenvectors;
    Eigen::VectorXd eigenvalues;


    /**
     *  All calculations from the perspective of one of the spin functions (no mixing)
     */
    void calculate_separated_elements(size_t string_state_one, size_t string_state_two,size_t address_state_one, size_t address_state_two, Eigen::MatrixXd& overlap, Eigen::MatrixXd& hamiltonian, std::vector<std::vector<OneElectronCoupling>>& coupling);

    /**
     *  All calculations exclusive to both the perspectives of the spin functions (mixing, only in two-electron terms)
     */
    void calculate_two_electron_mixed_elements();

    /**
     *  Calculates the overlap for two given states (possibly after operator evaluation)
     */
    double calculate_overlap(size_t string_state_one, size_t string_state_two);

    /**
     *  Calculates the overlap Matrix (overlap) and the (non-orthogonalized) Hamiltonian (hamiltonian)
     */
    void calculate_matrices();


public:
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
     * Solves the complete VB (non-orthognal fci) problem.
     */
    double solve();

    // GETTERS

    const Eigen::MatrixXd &get_eigenvectors() const {
        return this->eigenvectors;
    }

    const Eigen::VectorXd &get_eigenvalues() const {
        return this->eigenvalues;
    }

    void set_stuff(Eigen::MatrixXd oei, Eigen::Tensor<double,4> tei, Eigen::MatrixXd oi, std::vector<size_t> alpha,std::vector<size_t> beta){
        this->oei = oei;
        this->tei = tei;
        this->oi = oi;
        this->as = alpha;
        this->bs = beta;

    }
};

}  // namespace vb



#endif //VB_COMPLETEVB_HPP
