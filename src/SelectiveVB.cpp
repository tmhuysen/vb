#include "SelectiveVB.hpp"


namespace vb{

SelectiveVB::SelectiveVB(Eigen::MatrixXd oei, Eigen::Tensor<double,4> tei, Eigen::MatrixXd oi, std::vector<State> states):
        oei(std::move(oei)),tei(std::move(tei)),oi(std::move(oi)),states(std::move(states)) {
    this->dim = this->states.size();
}


void SelectiveVB::calculate_matrices() {
    this->hamiltonian = Eigen::MatrixXd::Zero(this->dim,this->dim);
    this->overlap = Eigen::MatrixXd::Zero(this->dim,this->dim);
    for(int index_row = 0; index_row<dim;index_row++){
        for(int index_col = index_row; index_col<dim;index_col++){
            for(int i = 0; i<states[index_row].alpha.size();i++){
                for(int j = 0; j<states[index_col].alpha.size();j++){
                    double overlap_alpha = calculate_overlap(states[index_row].alpha[i],states[index_col].alpha[j]);
                    double overlap_beta = calculate_overlap(states[index_row].beta[i],states[index_col].beta[j]);

                    this->overlap(index_row,index_col) += overlap_alpha*overlap_beta;
                    this->hamiltonian(index_row,index_col) += calculate_separated_elements(states[index_row].alpha[i],states[index_col].alpha[j])*overlap_beta;
                    this->hamiltonian(index_row,index_col) += calculate_separated_elements(states[index_row].beta[i],states[index_col].beta[j])*overlap_alpha;
                    this->hamiltonian(index_row,index_col) += calculate_mixed_elements(states[index_row].alpha[i],states[index_col].alpha[j],states[index_row].beta[i],states[index_col].beta[j]);
                    this->overlap(index_col,index_row)=this->overlap(index_row,index_col);
                    this->hamiltonian(index_col,index_row)=this->hamiltonian(index_row,index_col);

                }
            }
        }
    }

}

double SelectiveVB::calculate_overlap(size_t string_state_one, size_t string_state_two) {
    double overlaps = 1;
    size_t N = __builtin_popcount(string_state_one);
    Eigen::MatrixXd determinant;
    if(N != 0){
        determinant = Eigen::MatrixXd::Zero(N,N);
        size_t copy_one = string_state_one;
        size_t outerdex = 0;
        while (copy_one != 0) {
            size_t innerdex = 0;
            size_t p = __builtin_ctzl(copy_one);
            size_t copy_two = string_state_two;
            while (copy_two != 0){
                size_t q = __builtin_ctzl(copy_two);
                determinant(innerdex,outerdex) = this->oi(p,q);
                innerdex++;
                copy_two ^= (copy_two & -copy_two);
            }
            outerdex++;
            copy_one ^= (copy_one & -copy_one);
        }
        overlaps = determinant.determinant();
    }
    return overlaps;
}

double SelectiveVB::calculate_separated_elements(size_t string_state_one, size_t string_state_two) {

    double results = 0;
    // Calculate the amount of electrons for the string
    size_t N = __builtin_popcount(string_state_one);

    // Initialize the one-electron integral determinants (one determinant for each electron)
    std::vector<Eigen::MatrixXd> determinants(N);
    for(int i = 0;i<N;i++){
        determinants[i] = Eigen::MatrixXd::Zero(N,N);
    }

    size_t row_index = 0;  // Index related the determinant matrices
    size_t copy_one = string_state_one;  // copy bitstring
    while (copy_one != 0) {
        size_t column_index = 0;  // Index related the determinant matrices
        size_t p = __builtin_ctzl(copy_one);  // Index of the overlap or electron operator
        size_t copy_two = string_state_two;
        while (copy_two != 0){
            size_t q = __builtin_ctzl(copy_two);  // Index of the overlap or electron operator
            for(int i = 0;i<N;i++){
                if(row_index == i){
                    determinants[i](row_index,column_index) += this->oei(p,q);  // Fill one row with one electron-operator evaluations.
                }else{
                    determinants[i](row_index,column_index) += this->oi(p,q);  // Fill the rest of the rows with overlap terms
                }
            }

            // Resolve loop
            column_index++;
            copy_two ^= (copy_two & -copy_two);  // least significant bit removal.
        }
        row_index++;
        copy_one ^= (copy_one & -copy_one);  // least significant bit removal.
    }
    // Calculate the determinants and add them to the corresponding indexes
    for(const Eigen::MatrixXd &x:determinants){
        results += x.determinant();
    }
    // Two-electron part perform two annihilations on one_string and two_string, we can do this because all strings will couple (non-orthogonality).
    size_t copy1 = string_state_one;
    while (copy1 != 0) {
        int sign2 = 1; // first sign is always positive because sign with first annihilation is always the same for the annihilation directly afterwards (-1*-1 or 1*1 = 1)
        size_t p = __builtin_ctzl(copy1);  // Get operator index
        copy1 ^= (copy1 & -copy1);
        size_t copy2 = copy1;  // Annihilate only indexes bigger than previous annihilation
        while (copy2 != 0){
            int sign3 = sign2;  // Take over current sign (again no sign change for the first anni on the two_string)
            size_t q = __builtin_ctzl(copy2);  // Get operator index
            size_t copy3 = string_state_two;
            while(copy3 !=0){
                int sign4 = sign3;
                size_t r = __builtin_ctzl(copy3);  // Get operator index
                copy3 ^= (copy3 & -copy3);
                size_t copy4 = copy3;
                while(copy4 != 0){
                    size_t s = __builtin_ctzl(copy4);  // Get operator index

                    // Remove annihilated indexes for the remaining overlap strings
                    size_t aone = (string_state_one - (1<<p)) - (1<<q);
                    size_t atwo = (string_state_two - (1<<s)) - (1<<r);
                    double overlaps = calculate_overlap(aone, atwo);

                    // We only do one set of secondary annihilations, the opposite order results in a sign change (3rd and 4th term include one switch)
                    // double opposite switch for both strings results in a positive again (2nd term).
                    double term = tei(p,r,q,s) + tei(q,s,p,r) - tei(q,r,p,s) - tei(p,s,q,r);
                    results += sign4*overlaps*term/2;

                    copy4 ^= (copy4 & -copy4);
                    sign4 *=-1;  // For subsequent secondary annihilations change the sign.
                }
            }
            sign2 *= -1; // For subsequent secondary annihilations change the sign.
            copy2 ^= (copy2 & -copy2);
        }
    }
    return results;
}

double SelectiveVB::calculate_mixed_elements(size_t aone_string, size_t atwo_string, size_t bone_string, size_t btwo_string) {
    double results = 0;

    size_t copy1 = aone_string;
    int sign1 = 1;
    while (copy1 != 0) {
        int sign2 = sign1;
        size_t p = __builtin_ctzl(copy1);
        size_t copy2 = bone_string;
        while (copy2 != 0){
            int sign3 = sign2;
            size_t q = __builtin_ctzl(copy2);
            size_t copy3 = btwo_string;
            while(copy3 !=0){
                int sign4 = sign3;
                size_t r = __builtin_ctzl(copy3);
                size_t copy4 = atwo_string;
                while(copy4 != 0){
                    size_t s = __builtin_ctzl(copy4);



                    size_t aone = aone_string - (1<<p);
                    size_t bone = bone_string - (1<<q);
                    size_t btwo = btwo_string - (1<<r);
                    size_t atwo = atwo_string - (1<<s);

                    double overlaps = calculate_overlap(aone, atwo)* calculate_overlap(bone, btwo);
                    results += sign4*overlaps*tei(p,s,q,r);




                    copy4 ^= (copy4 & -copy4);
                    sign4 *=-1;

                }

                sign3 *=-1;
                copy3 ^= (copy3 & -copy3);
            }

            sign2 *= -1;
            copy2 ^= (copy2 & -copy2);

        }
        sign1 *= -1;
        copy1 ^= (copy1 & -copy1);
    }
    return results;
}

double SelectiveVB::solve() {
    this->calculate_matrices();
    if(dim==1){
        return hamiltonian(0)/overlap(0);
    }

    // We reduced the problem to an eigenvalue problem, here we solve it.
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> eigen_solver (this->hamiltonian, this->overlap);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver2(this->overlap);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver3(solver2.operatorInverseSqrt()*this->hamiltonian*solver2.operatorInverseSqrt());

    this->eigenvectors = solver3.eigenvectors();
    this->eigenvalues = solver3.eigenvalues();
    return eigen_solver.eigenvalues()(0);
}


}  // namespace vb