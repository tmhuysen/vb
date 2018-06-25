#ifndef VB_SelectiveVBweight_HPP
#define VB_SelectiveVBweight_HPP

#include "SelectiveVB.hpp"
#include "libwint.hpp"
#include "common.hpp"

namespace vb {
/*
 * SelectiveVB derivation where you can add fixed custom weights.
 */

class SelectiveVBweight : public SelectiveVB {
private:

    std::vector<StateW> states;

    /**
     *  Calculates the overlap Matrix (overlap) and the (non-orthogonalized) Hamiltonian (hamiltonian)
     */
    void calculate_matrices() override;


public:
    // CONSTRUCTORS

    /**
     *  Constructor based on a given @param ao_basis, a number of alpha electrons @param N_alpha and a number of beta electrons
     *  @param N_beta.
     */
    SelectiveVBweight(Eigen::MatrixXd oei, Eigen::Tensor<double,4> tei, Eigen::MatrixXd oi, std::vector<StateW> states);


};

}
#endif //VB_SelectiveVBweight_HPP
