#include "SelectiveVBweight.hpp"


namespace vb{




void SelectiveVBweight::calculate_matrices() {
    this->hamiltonian = Eigen::MatrixXd::Zero(this->dim,this->dim);
    this->overlap = Eigen::MatrixXd::Zero(this->dim,this->dim);
    for(int index_row = 0; index_row<dim;index_row++){
        for(int index_col = index_row; index_col<dim;index_col++){
            for(int i = 0; i<states[index_row].alpha.size();i++){
                for(int j = 0; j<states[index_col].alpha.size();j++){
                    double weights = states[index_row].weights[i]*states[index_col].weights[j];
                    double overlap_alpha = calculate_overlap(states[index_row].alpha[i],states[index_col].alpha[j]);
                    double overlap_beta = calculate_overlap(states[index_row].beta[i],states[index_col].beta[j]);

                    this->overlap(index_row,index_col) += overlap_alpha*overlap_beta*weights;
                    double hamvalue = 0;
                    hamvalue += calculate_separated_elements(states[index_row].alpha[i],states[index_col].alpha[j])*overlap_beta;
                    hamvalue += calculate_separated_elements(states[index_row].beta[i],states[index_col].beta[j])*overlap_alpha;
                    hamvalue += calculate_mixed_elements(states[index_row].alpha[i],states[index_col].alpha[j],states[index_row].beta[i],states[index_col].beta[j]);

                    this->hamiltonian(index_row,index_col) += hamvalue*weights;
                    this->overlap(index_col,index_row)=this->overlap(index_row,index_col);
                    this->hamiltonian(index_col,index_row)=this->hamiltonian(index_row,index_col);

                }
            }
        }
    }

}

SelectiveVBweight::SelectiveVBweight(Eigen::MatrixXd oei, Eigen::Tensor<double, 4> tei, Eigen::MatrixXd oi,
                                     std::vector<StateW> states) : SelectiveVB(oei,tei,oi,std::vector<State>{}),
                                     states(states){
    this->dim = this->states.size();

}

}  // namespace vb