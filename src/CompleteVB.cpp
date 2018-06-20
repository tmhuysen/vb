#include "CompleteVB.hpp"
#include "bits.hpp"
#include <chrono>

namespace vb{

CompleteVB::CompleteVB(libwint::AOBasis& ao_basis, size_t N_A, size_t N_B):
        ao_basis(ao_basis), N_alpha(N_A), N_beta(N_B),K(ao_basis.calculateNumberOfBasisFunctions()),
        dim_alpha (CompleteVB::calculateDimension(ao_basis.calculateNumberOfBasisFunctions(), N_A, 0)),
        dim_beta (CompleteVB::calculateDimension(ao_basis.calculateNumberOfBasisFunctions(), 0, N_B)),
        dim(CompleteVB::calculateDimension(ao_basis.calculateNumberOfBasisFunctions(), N_A, N_B)) {

    this->oei = ao_basis.get_T()+ao_basis.get_V();
    this->oi = ao_basis.get_S();
    this->tei = ao_basis.get_g();

    this->hamiltonian_beta = Eigen::MatrixXd::Zero(dim_beta,dim_beta);
    this->hamiltonian_alpha = Eigen::MatrixXd::Zero(dim_alpha,dim_alpha);
    this->overlap_alpha = Eigen::MatrixXd::Zero(dim_alpha,dim_alpha);
    this->overlap_beta = Eigen::MatrixXd::Zero(dim_beta,dim_beta);

    this->as = std::vector<size_t>(dim_alpha);
    size_t ssa_a = vb::smallest_bitset(this->N_alpha);
    for (size_t i = 0; i < this->dim_alpha; i++) {
        this->as[i] = ssa_a;
        ssa_a = next_bitset_permutation(ssa_a);
    }

    this->bs = std::vector<size_t>(dim_beta);
    size_t ssa_b = vb::smallest_bitset(this->N_beta);
    for (size_t i = 0; i < this->dim_beta; i++) {
        bs[i] = ssa_b;
        ssa_b = next_bitset_permutation(ssa_b);
    }
    this->aoec = std::vector<std::vector<OneElectronCoupling>>(dim_alpha);
    this->boec = std::vector<std::vector<OneElectronCoupling>>(dim_beta);

}

void CompleteVB::calculate_matrices() {
    overlap = Eigen::MatrixXd::Zero(dim,dim);
    hamiltonian = Eigen::MatrixXd::Zero(dim,dim);

    for (size_t I_alpha = 0; I_alpha < this->dim_alpha; I_alpha++) {
        for (size_t J_alpha = 0; J_alpha < this->dim_alpha; J_alpha++){
            calculate_separated_elements(as[I_alpha], as[J_alpha], I_alpha, J_alpha,this->overlap_alpha,this->hamiltonian_alpha,aoec);
        }
    }

    for (size_t I_beta = 0; I_beta < this->dim_beta; I_beta++) {
        for (size_t J_beta = 0; J_beta < this->dim_beta; J_beta++){
            calculate_separated_elements(bs[I_beta], bs[J_beta], I_beta, J_beta,this->overlap_beta,this->hamiltonian_beta,boec);
        }
    }

    for(int i = 0;i<this->dim_alpha;i++){
        for(int j = 0;j<this->dim_alpha;j++) {
            overlap.block(i*dim_beta, j*dim_beta,dim_beta,dim_beta)    += overlap_alpha(i,j)*overlap_beta;
            hamiltonian.block(i*dim_beta,j*dim_beta,dim_beta,dim_beta) += overlap_alpha(i,j)*hamiltonian_beta;
            hamiltonian.block(i*dim_beta,j*dim_beta,dim_beta,dim_beta) += hamiltonian_alpha(i,j)*overlap_beta;
        }
    }
    this->calculate_two_electron_mixed_elements();

}

size_t CompleteVB::calculateDimension(size_t K, size_t N_alpha, size_t N_beta) {
    // K, N_alpha, N_beta are expected to be small, so static-casting them to unsigned (what boost needs) is permitted.
    auto dim_double_alpha = boost::math::binomial_coefficient<double>(static_cast<unsigned>(K), static_cast<unsigned>(N_alpha));
    auto dim_double_beta = boost::math::binomial_coefficient<double>(static_cast<unsigned>(K), static_cast<unsigned>(N_beta));
    auto dim_double = dim_double_alpha * dim_double_beta;

    // Check if the resulting dimension is appropriate to be stored in size_t
    return boost::numeric::converter<double, size_t>::convert(dim_double);
}

double CompleteVB::solve() {
    this->calculate_matrices();

    // We reduced the problem to an eigenvalue problem, here we solve it.
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> eigen_solver (this->hamiltonian, this->overlap);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver2(this->overlap);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver3(solver2.operatorInverseSqrt()*this->hamiltonian*solver2.operatorInverseSqrt());

    return eigen_solver.eigenvalues()(0);
}

void CompleteVB::calculate_two_electron_mixed_elements(){
    for(int a = 0;a<dim_alpha;a++){
        for(OneElectronCoupling x: aoec[a]){
            for(int b = 0;b<dim_beta;b++){
                for(OneElectronCoupling y: boec[b]){
                    int sign = x.sign*y.sign;
                    double overlaps = x.overlap*y.overlap;
                    hamiltonian(a*dim_beta+b,x.address_target*dim_beta+y.address_target) += sign*overlaps*tei(x.p,x.q,y.p,y.q);
                }
            }
        }
    }
}

double CompleteVB::calculate_overlap(size_t string_state_one, size_t string_state_two) {
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

void CompleteVB::calculate_separated_elements(size_t string_state_one, size_t string_state_two, size_t address_state_one, size_t address_state_two,
                                              Eigen::MatrixXd& overlap, Eigen::MatrixXd& hamiltonian, std::vector<std::vector<OneElectronCoupling>>& coupling_vector) {
    
    // Calculate the amount of electrons for the string
    size_t N = __builtin_popcount(string_state_one);
    
    // Initialize the overlap determinant
    Eigen::MatrixXd determinant_overlap = Eigen::MatrixXd::Zero(N,N);
    
    // Initialize the one-electron integral determinants (one determinant for each electron)
    std::vector<Eigen::MatrixXd> determinants(N);
    for(int i = 0;i<N;i++){
        determinants[i] = Eigen::MatrixXd::Zero(N,N);
    }

    size_t row_index = 0;  // Index related the determinant matrices
    size_t copy_one = string_state_one;  // copy bitstring
    int signin =1;
    while (copy_one != 0) {
        int signout = signin;
        size_t column_index = 0;  // Index related the determinant matrices
        size_t p = __builtin_ctzl(copy_one);  // Index of the overlap or electron operator
        size_t copy_two = string_state_two;
        while (copy_two != 0){
            size_t q = __builtin_ctzl(copy_two);  // Index of the overlap or electron operator
            determinant_overlap(row_index,column_index) += this->oi(p,q);  // The determinant for the overlap matrix
            for(int i = 0;i<this->N_alpha;i++){
                if(row_index == i){
                    determinants[i](row_index,column_index) += this->oei(p,q);  // Fill one row with one electron-operator evaluations.
                }else{
                    determinants[i](row_index,column_index) += this->oi(p,q);  // Fill the rest of the rows with overlap terms
                }
            }
            // Important step for two-electron mixed calculations (stores the overlap for a given operator combination)
            double overlap_coupling = calculate_overlap(string_state_one-(1<<p),string_state_two-(1<<q));
            coupling_vector[address_state_one].emplace_back(OneElectronCoupling{signout,p,q,overlap_coupling,address_state_two});

            // Resolve loop
            column_index++;
            copy_two ^= (copy_two & -copy_two);  // least significant bit removal.
            signout *=-1;
        }
        row_index++;
        copy_one ^= (copy_one & -copy_one);  // least significant bit removal.
        signin *= -1;
    }
    // Calculate the determinants and add them to the corresponding indexes
    for(const Eigen::MatrixXd &x:determinants){
        hamiltonian(address_state_one,address_state_two) += x.determinant();
    }
    // Add the overlap
    overlap(address_state_one,address_state_two) = determinant_overlap.determinant();
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
                    hamiltonian(address_state_one,address_state_two) += sign4*overlaps*term/2;

                    copy4 ^= (copy4 & -copy4);
                    sign4 *=-1;  // For subsequent secondary annihilations change the sign.
                }
            }
            sign2 *= -1; // For subsequent secondary annihilations change the sign.
            copy2 ^= (copy2 & -copy2);
        }
    }
}
}  // namespace vb