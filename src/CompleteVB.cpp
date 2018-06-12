#include "CompleteVB.hpp"

namespace vb{

CompleteVB::CompleteVB(libwint::AOBasis& ao_basis, size_t N_A, size_t N_B):
        ao_basis(ao_basis), N_alpha(N_A), N_beta(N_B),K(ao_basis.calculateNumberOfBasisFunctions()),
        addressing_scheme_alpha (bmqc::AddressingScheme(ao_basis.calculateNumberOfBasisFunctions(), N_A)),
        addressing_scheme_beta (bmqc::AddressingScheme(ao_basis.calculateNumberOfBasisFunctions(), N_A)),
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
}

void CompleteVB::calculate_matrices() {
    overlap = Eigen::MatrixXd::Zero(dim,dim);
    hamiltonian = Eigen::MatrixXd::Zero(dim,dim);


    bmqc::SpinString<unsigned long> ssa_one (0, this->addressing_scheme_alpha);  // alpha spin string with address 0
    for (size_t I_alpha = 0; I_alpha < this->dim_alpha; I_alpha++) {
        bmqc::SpinString<unsigned long> ssa_two (0, this->addressing_scheme_alpha);  // alpha spin string with address 0
        for (size_t J_alpha = 0; J_alpha < this->dim_alpha; J_alpha++){
            calculate_element_alpha(ssa_one.get_representation(), ssa_two.get_representation(), I_alpha, J_alpha);


            // BETA-ALPHA
            bmqc::SpinString<unsigned long> ssa_2one (0, this->addressing_scheme_beta);  // beta spin string with address 0
            for (size_t I_beta = 0; I_beta < this->dim_beta; I_beta++) {
                bmqc::SpinString<unsigned long> ssa_2two (0, this->addressing_scheme_beta);  // beta spin string with address 0
                for (size_t J_beta = 0; J_beta < this->dim_beta; J_beta++) {
                    calculate_element(ssa_one.get_representation(), ssa_two.get_representation(),ssa_2one.get_representation(), ssa_2two.get_representation(), this->dim_beta*I_alpha + I_beta, this->dim_beta*J_alpha+J_beta);
                    ssa_2two.nextPermutation();
                }
                ssa_2one.nextPermutation();
            }





            ssa_two.nextPermutation();



        }
        ssa_one.nextPermutation();
    }


    bmqc::SpinString<unsigned long> ssa_2one (0, this->addressing_scheme_beta);  // beta spin string with address 0
    for (size_t I_beta = 0; I_beta < this->dim_beta; I_beta++) {
        bmqc::SpinString<unsigned long> ssa_2two (0, this->addressing_scheme_beta);  // beta spin string with address 0
        for (size_t J_beta = 0; J_beta < this->dim_beta; J_beta++){
            calculate_element_beta(ssa_2one.get_representation(), ssa_2two.get_representation(), I_beta, J_beta);
            ssa_2two.nextPermutation();



        }
        ssa_2one.nextPermutation();
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

void CompleteVB::calculate_element_alpha(size_t one_string, size_t two_string, size_t index_one, size_t index_two) {
    // Initialize a matrix t
    Eigen::MatrixXd determinant_overlap = Eigen::MatrixXd::Zero(this->N_alpha,this->N_alpha);
    Eigen::MatrixXd determinant_one = Eigen::MatrixXd::Zero(this->N_alpha,this->N_alpha);
    Eigen::MatrixXd determinant_two = Eigen::MatrixXd::Zero(this->N_alpha,this->N_alpha);
    size_t outerdex = 0;
    size_t copy_one = one_string;
    while (copy_one != 0) {
        size_t innerdex = 0;
        size_t p = __builtin_ctzl(copy_one);
        size_t copy_two = two_string;
        while (copy_two != 0){
            size_t q = __builtin_ctzl(copy_two);
            determinant_overlap(innerdex,outerdex) += this->oi(p,q);
            if(outerdex==0){
                determinant_one(innerdex,outerdex) += this->oei(p,q);

            }else{
                determinant_one(innerdex,outerdex) += this->oi(p,q);
            }

            innerdex++;
            copy_two ^= (copy_two & -copy_two);

        }
        outerdex++;
        copy_one ^= (copy_one & -copy_one);
    }

    this->hamiltonian_alpha(index_one,index_two) = determinant_one.determinant();
    this->overlap_alpha(index_one,index_two) = determinant_overlap.determinant();


    size_t copy1 = one_string;
    int sign1 = 1;
    while (copy1 != 0) {
        int sign2 = sign1;
        size_t p = __builtin_ctzl(copy1);
        size_t copy2 = one_string;
        while (copy2 != 0){
            int sign3 = sign2;
            size_t q = __builtin_ctzl(copy2);
            size_t copy3 = two_string;
            while(copy3 !=0){
                int sign4 = sign3;
                size_t r = __builtin_ctzl(copy3);
                size_t copy4 = two_string;
                while(copy4 != 0){
                    size_t s = __builtin_ctzl(copy4);


                    if(p == q || r == s){

                    }else{
                        size_t aone = one_string - (1<<p) - (1<<q);
                        size_t atwo = two_string - (1<<s) - (1<<r);
                        double overlaps = calculate_overlap(aone,atwo);
                        hamiltonian_alpha(index_one,index_two) += sign4*overlaps*tei(p,s,q,r)/2*0;

                    }
                    copy4 ^= (copy4 & -copy4);
                    if(s<r){
                        sign4 *=-1;
                    }

                }

                sign3 *=-1;
                copy3 ^= (copy3 & -copy3);
            }
            if(p<q){
                sign2 *= -1;
            }
            copy2 ^= (copy2 & -copy2);

        }
        sign1 *= -1;
        copy1 ^= (copy1 & -copy1);
    }



}


void CompleteVB::calculate_element_beta(size_t one_string, size_t two_string, size_t index_one, size_t index_two) {
    // Initialize a matrix t
    Eigen::MatrixXd determinant_overlap = Eigen::MatrixXd::Zero(this->N_beta,this->N_beta);
    Eigen::MatrixXd determinant_one = Eigen::MatrixXd::Zero(this->N_beta,this->N_beta);
    Eigen::MatrixXd determinant_two = Eigen::MatrixXd::Zero(this->N_beta,this->N_beta);
    size_t outerdex = 0;
    size_t copy_one = one_string;
    while (copy_one != 0) {
        size_t innerdex = 0;
        size_t p = __builtin_ctzl(copy_one);
        size_t copy_two = two_string;
        while (copy_two != 0){
            size_t q = __builtin_ctzl(copy_two);
            determinant_overlap(innerdex,outerdex) = this->oi(p,q);
            if(outerdex==0){
                determinant_one(innerdex,outerdex) = this->oei(p,q);

            }else{
                determinant_one(innerdex,outerdex) = this->oi(p,q);
            }

            innerdex++;
            copy_two ^= (copy_two & -copy_two);

        }
        outerdex++;
        copy_one ^= (copy_one & -copy_one);
    }

    this->hamiltonian_beta(index_one,index_two) = determinant_one.determinant();
    this->overlap_beta(index_one,index_two) = determinant_overlap.determinant();


    size_t copy1 = one_string;
    int sign1 = 1;
    while (copy1 != 0) {
        int sign2 = sign1;
        size_t p = __builtin_ctzl(copy1);
        size_t copy2 = one_string;
        while (copy2 != 0){
            int sign3 = sign2;
            size_t q = __builtin_ctzl(copy2);
            size_t copy3 = two_string;
            while(copy3 !=0){
                int sign4 = sign3;
                size_t r = __builtin_ctzl(copy3);
                size_t copy4 = two_string;
                while(copy4 != 0){
                    size_t s = __builtin_ctzl(copy4);


                    if(p == q || r == s){

                    }else{
                        size_t aone = (one_string - (1<<p)) - (1<<q);
                        size_t atwo = (two_string - (1<<s)) - (1<<r);
                        double overlaps = calculate_overlap(aone,atwo);
                        std::cout<<-sign4<<" "<<p<<" "<<q<<" "<<r<<" "<<s<<std::endl;
                        hamiltonian_beta(index_one,index_two) += -sign4*overlaps*tei(s,p,r,q)/2;
                        hamiltonian_alpha(index_one,index_two) += -sign4*overlaps*tei(p,s,q,r)/2;

                    }
                    copy4 ^= (copy4 & -copy4);
                    if(s>r){
                        sign4 *=-1;

                    }else{
                        //sign4 *= sign3;
                    }

                }

                sign3 *=-1;
                copy3 ^= (copy3 & -copy3);
            }
            if(p<q){
                sign2 *= -1;
            }else{
                //sign2 = sign2*sign1;
            }
            copy2 ^= (copy2 & -copy2);

        }
        sign1 *= -1;
        copy1 ^= (copy1 & -copy1);
    }





}

void CompleteVB::finish_off() {
    for(int i = 0;i<this->dim_alpha;i++){
        for(int j = 0;j<this->dim_alpha;j++) {
            overlap.block(i*dim_beta, j*dim_beta,dim_beta,dim_beta)    += overlap_alpha(i,j)*overlap_beta;
            hamiltonian.block(i*dim_beta,j*dim_beta,dim_beta,dim_beta) += overlap_alpha(i,j)*hamiltonian_beta;
            hamiltonian.block(i*dim_beta,j*dim_beta,dim_beta,dim_beta) += hamiltonian_alpha(i,j)*overlap_beta;
        }
    }
}

double CompleteVB::solve() {
    calculate_matrices();
    finish_off();
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> eigen_solver (this->hamiltonian, this->overlap);

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver2(this->overlap);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver3(solver2.operatorInverseSqrt()*this->hamiltonian*solver2.operatorInverseSqrt());

    std::cout<<eigen_solver.eigenvalues().rows();
    return eigen_solver.eigenvalues()(0);
}

void CompleteVB::calculate_element(size_t aone_string, size_t atwo_string, size_t bone_string, size_t btwo_string,
                                   size_t index_one, size_t index_two) {

    // Initialize a matrix t
    // Eigen::MatrixXd determinant_two = Eigen::MatrixXd::Zero(this->N_beta,this->N_beta);
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

                    double overlaps = calculate_overlap(aone,atwo)*calculate_overlap(bone,btwo);
                    hamiltonian(index_one,index_two) += sign4*overlaps*tei(p,s,q,r);




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



}

double CompleteVB::calculate_overlap(size_t one_string, size_t two_string) {
    double overlaps = 1;
    size_t N = __builtin_popcount(one_string);
    Eigen::MatrixXd determinant;
    if(N != 0){
        determinant = Eigen::MatrixXd::Zero(N,N);
        size_t copy_one = one_string;
        size_t outerdex = 0;
        while (copy_one != 0) {
            size_t innerdex = 0;
            size_t p = __builtin_ctzl(copy_one);
            size_t copy_two = two_string;
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
}  // namespace vb