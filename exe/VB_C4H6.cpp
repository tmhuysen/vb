/*
 * Hard-coded driver for C4H6 VB calculations
 */

#include <libwint.hpp>
#include "common.hpp"
#include <hf.hpp>
#include <iomanip>
#include "SelectiveVB.hpp"

int main() {

    std::cout<<std::setprecision(12);  // Set print precision to 12

    //  Retrieve geomrty, caclulate integrals etc.
    std::string filename = "../../exe/C4H6.xyz";
    libwint::Molecule C4H6 (filename);
    libwint::AOBasis ao_basis (C4H6, "STO-3G");  // init Atomic Orbital basis
    double repulsion = C4H6.calculateInternuclearRepulsionEnergy();  // repulsion
    ao_basis.calculateIntegrals();  // Calculate integrals
    hf::rhf::RHF rhf (C4H6, ao_basis, 1.0e-013);  // init RHF
    rhf.solve( hf::rhf::solver::SCFSolverType::DIIS);  // Solve RHF
    std::cout<<std::endl<<"RHF (30 electrons):"<<rhf.get_electronic_energy()+repulsion<<std::endl;  // print RHF energy
    libwint::SOBasis so_basis(ao_basis,rhf.get_C_canonical());  // Generate orthogal basis using the canonical matrix

    Eigen::MatrixXd one_int_so = so_basis.get_h_SO();  // Retrieve one electron MO integrals
    Eigen::Tensor<double,4> two_int_so = so_basis.get_g_SO();  // Retrieve two electron MO integrals

    // Calculate 26 doubly occupied single slater (test) with the 30 electron RHF basis.
    int N = 26;
    double energy = 0;
    for(size_t i = 0;i<N/2;i++){
        energy += 2*one_int_so(i,i);
        energy += two_int_so(i,i,i,i);
        for(size_t j = 0;j<i;j++){
            energy += 4*two_int_so(i,i,j,j);
            energy -= 2*two_int_so(i,j,j,i);
        }
    }
    energy += repulsion;
    std::cout<<"RHF (26 electrons):"<<energy<<std::endl;  // print RHF energy with 26 electrons
    int num_vb_bf = 17;  // number of basis functions used in the VB

    // Prepare the custom_basis for our calculation (rhf+ao mix)
    Eigen::MatrixXd oei = Eigen::MatrixXd::Zero(num_vb_bf,num_vb_bf);  // oei = one electron integrals
    Eigen::MatrixXd oi = Eigen::MatrixXd::Identity(num_vb_bf,num_vb_bf);  // overlap integrals
    Eigen::MatrixXd one_inter = ao_basis.get_V()+ao_basis.get_T();  // intermediate save of atomic one electron orbitals (kinetic + potential)
    Eigen::Tensor<double,4> two_ei(num_vb_bf,num_vb_bf,num_vb_bf,num_vb_bf);  // two electron integrals for the custom.
    two_ei.setZero();
    int core_number = 13;  // RHF basis functions used
    int pz_number = 4;  // Pz basis functions used

    // Copying MO basis in to custom basis
    oei.block(0,0,core_number,core_number) = so_basis.get_h_SO().block(0,0,core_number,core_number);
    Eigen::Tensor<double,4> two_electron_core = so_basis.get_g_SO();
    for(int i = 0;i<core_number;i++){
        for (int j = 0;j<core_number;j++){
            for (int k =0;k<core_number;k++){
                for(int l =0;l<core_number;l++){
                    two_ei(i,j,k,l) = two_electron_core(i,j,k,l);
                }
            }
        }
    }

    // Copying AO basis into custom
    std::vector<int> pz_indexes = {10,15,20,25};
    Eigen::Tensor<double,4> ao_two_electron = ao_basis.get_g();
    for(int i = 0;i<pz_number;i++){
        for(int j = 0;j<pz_number;j++){
            oei(i+core_number,j+core_number) = one_inter(pz_indexes[i],pz_indexes[j]);
            oi(i+core_number,j+core_number) = ao_basis.get_S()(pz_indexes[i],pz_indexes[j]);
            for (int k =0;k<pz_number;k++){
                for(int l =0;l<pz_number;l++){
                    two_ei(i+core_number,j+core_number,k+core_number,l+core_number) = ao_two_electron(pz_indexes[i],pz_indexes[j],pz_indexes[k],pz_indexes[l]);

                }
            }
        }
    }

    // Calculating two-electron integrals between RHF and AO basis for custom basis

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solvy(ao_basis.get_S());
    //Eigen::MatrixXd C = solvy.operatorSqrt()*rhf.get_C_canonical();
    Eigen::MatrixXd C = rhf.get_C_canonical();
    for(int i = 0;i<pz_number;i++){
        for(int j = 0; j<pz_number;j++){
            for(int k = 0;k<core_number;k++){
                for (int l = 0;l<core_number;l++){
                    for(int ao_index1 = 0;ao_index1<26;ao_index1++){
                        for(int ao_index2 = 0; ao_index2<26;ao_index2++){
                            // all (Pz,Pz,Pz,Rhf) or (Pz,Rhf,Rhf,Rhf) are 0.
                            two_ei(i+core_number,j+core_number,k,l) += C(ao_index1,k)*C(ao_index2,l)*ao_two_electron(pz_indexes[i],pz_indexes[j],ao_index1,ao_index2);
                            two_ei(k,l,i+core_number,j+core_number) =  two_ei(i+core_number,j+core_number,k,l);

                            two_ei(i+core_number,k,j+core_number,l) += C(ao_index1,k)*C(ao_index2,l)*ao_two_electron(pz_indexes[i],ao_index1,pz_indexes[j],ao_index2);
                            two_ei(i+core_number,k,l,j+core_number) = two_ei(i+core_number,k,j+core_number,l);
                            two_ei(k,i+core_number,j+core_number,l) = two_ei(i+core_number,k,j+core_number,l);
                            two_ei(k,i+core_number,l,j+core_number) = two_ei(i+core_number,k,j+core_number,l);


                        }
                    }
                }
            }
        }
    }

    // Core electrons indexes
    std::vector<int> baseline = {0,1,2,3,4,5,6,7,8,9,10,11,12};
    size_t base_det = vb::bitstring(baseline);

    // Butadiene Pz orbitals are indexed as follows:
    //
    //
    //                    C4
    //       C1--_____   //
    //     //         --C2
    //    C3
    //
    // C1 = 13
    // C2 = 14
    // C3 = 15
    // C4 = 16

    // Possible determinants:

    size_t det1 = base_det + vb::bitstring({13,14});
    size_t det2 = base_det + vb::bitstring({13,16});
    size_t det3 = base_det + vb::bitstring({14,15});
    size_t det4 = base_det + vb::bitstring({15,16});
    size_t det5 = base_det + vb::bitstring({13,15});
    size_t det6 = base_det + vb::bitstring({14,16});
    std::vector<size_t> slaters = {det1,det2,det3,det4,det5,det6};
    vb::State cov {{det1,det2,det3,det4},{det4,det3,det2,det1},{-1,1,1,-1}};

    vb::SelectiveVB selectiveVB(oei,two_ei,oi,{cov});
    selectiveVB.orthogonality_set = 13;
    std::cout<<std::endl<<" complete pz "<<" : "<<selectiveVB.solve()+repulsion<<std::endl;
    //igen::VectorXd veck = selectiveVB.get_eigenvectors().col(0);
}
