#ifndef VB_SELECTIVEVB_HPP
#define VB_SELECTIVEVB_HPP

#include "libwint.hpp"
#include "common.hpp"

namespace vb {

class SelectiveVB {
protected:
    size_t dim;
    size_t dim_alpha;  // the dimension of the alpha CI space
    size_t dim_beta;  // the dimension of the beta CI space

    size_t K;  // number of spatial orbitals

    size_t N_alpha;  // number of alpha electrons
    size_t N_beta;  // number of beta electrons

    Eigen::MatrixXd oei;
    Eigen::Tensor<double,4> tei;
    Eigen::MatrixXd oi;

    std::vector<State> states;

    Eigen::MatrixXd hamiltonian;
    Eigen::MatrixXd overlap;

    Eigen::MatrixXd eigenvectors;
    Eigen::VectorXd eigenvalues;








    /**
     *  All calculations from the perspective of one of the spin functions (no mixing)
     */
    double calculate_separated_elements(size_t string_state_one, size_t string_state_two);

    /**
     *  All calculations exclusive to both the perspectives of the spin functions (mixing, only in two-electron terms)
     */
    double calculate_mixed_elements(size_t aone_string, size_t atwo_string, size_t bone_string, size_t btwo_string);

    /**
     *  Calculates the overlap for two given states (possibly after operator evaluation)
     */
    double calculate_overlap(size_t string_state_one, size_t string_state_two);

     /**
     *  Calculates the overlap Matrix (overlap) and the (non-orthogonalized) Hamiltonian (hamiltonian)
     */
    virtual void calculate_matrices();


public:
    // CONSTRUCTORS

    /**
     *  Constructor based on a given @param ao_basis, a number of alpha electrons @param N_alpha and a number of beta electrons
     *  @param N_beta.
     */
    SelectiveVB(Eigen::MatrixXd oei, Eigen::Tensor<double,4> tei, Eigen::MatrixXd oi, std::vector<State> states);
    SelectiveVB()=default;

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

};

}
#endif //VB_SELECTIVEVB_HPP
