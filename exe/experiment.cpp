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
    int num_vb_bf = 5;  // number of basis functions used in the VB

    // Prepare the custom_basis for our calculation (rhf+ao mix)
    Eigen::MatrixXd oei = Eigen::MatrixXd::Zero(num_vb_bf,num_vb_bf);  // oei = one electron integrals
    Eigen::MatrixXd oi = Eigen::MatrixXd::Identity(num_vb_bf,num_vb_bf);  // overlap integrals
    Eigen::MatrixXd one_inter = ao_basis.get_V()+ao_basis.get_T();  // intermediate save of atomic one electron orbitals (kinetic + potential)
    Eigen::Tensor<double,4> two_ei(num_vb_bf,num_vb_bf,num_vb_bf,num_vb_bf);  // two electron integrals for the custom.
    two_ei.setZero();
    int core_number = 13;  // RHF basis functions used
    int pz_number = 4;  // Pz basis functions used

    // Copying MO basis in to custom basis
    for(int i =0; i<core_number; i++){
        oei(0,0) += so_basis.get_h_SO()(i,i);
    }
    oi(0,0) = 13;
    Eigen::Tensor<double,4> two_electron_core = so_basis.get_g_SO();
    for(int i = 0;i<core_number;i++){
        for (int j = 0;j<core_number;j++){
            for (int k =0;k<core_number;k++){
                for(int l =0;l<core_number;l++){
                    two_ei(0,0,0,0) += two_electron_core(i,j,k,l);
                }
            }
        }
    }
    int core_total=1;
    // Copying AO basis into custom
    std::vector<int> pz_indexes = {10,15,20,25};
    Eigen::Tensor<double,4> ao_two_electron = ao_basis.get_g();
    for(int i = 0;i<pz_number;i++){
        for(int j = 0;j<pz_number;j++){
            oei(i+core_total,j+core_total) = one_inter(pz_indexes[i],pz_indexes[j]);
            oi(i+core_total,j+core_total) = ao_basis.get_S()(pz_indexes[i],pz_indexes[j]);
            for (int k =0;k<pz_number;k++){
                for(int l =0;l<pz_number;l++){
                    two_ei(i+core_total,j+core_total,k+core_total,l+core_total) = ao_two_electron(pz_indexes[i],pz_indexes[j],pz_indexes[k],pz_indexes[l]);
                }
            }
        }
    }

    // Calculating two-electron integrals between RHF and AO basis for custom basis
    Eigen::MatrixXd C = rhf.get_C_canonical();

    for(int i = 0;i<pz_number;i++){
        for(int j = 0; j<pz_number;j++){
            for(int k = 0;k<core_number;k++){
                for (int l = 0;l<core_number;l++){
                    for(int ao_index1 = 0;ao_index1<26;ao_index1++){
                        for(int ao_index2 = 0; ao_index2<26;ao_index2++){
                            // Only have to calculate (Pz,Pz,Rhf,Rhf) because all (Pz,Rhf,Pz,Rhf) are 0.
                            two_ei(i+core_total,j+core_total,0,0) += C(ao_index1,k)*C(ao_index2,l)*ao_two_electron(pz_indexes[i],pz_indexes[j],ao_index1,ao_index2);
                            two_ei(0,0,j+core_total,i+core_total) =  two_ei(i+core_number,j+core_number,k,l);
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

    size_t det1 = vb::bitstring({0,3,4});
    size_t det2 = vb::bitstring({0,1,2});
    size_t det3 = vb::bitstring({0,2,3});
    size_t det4 = vb::bitstring({0,1,4});

    vb::State onevibe{{det1,det2,det3,det4},{det2,det1,det4,det3}};
    try{
        vb::SelectiveVB selectiveVB(oei,two_ei,oi,{onevibe});
        std::cout<<std::endl<<"vb:"<<selectiveVB.solve()+repulsion<<std::endl;
    }catch (const std::overflow_error& e){

    }

    /*
    try{
        vb::State onevibe{{det1,det2,det3,det4},{det2,det1,det4,det3}};
        vb::SelectiveVB selectiveVB(oei,g,oi,{onevibe});
        std::cout<<std::endl<<"vb:"<<selectiveVB.solve()+repulsion<<std::endl;
    }catch (const std::overflow_error& e){

    }




    std::vector<size_t> becky = {det1,det2,det3,det4};
    int myints[] = {0,1,2,3};
    std::sort (myints,myints+4);
    std::cout<<std::setprecision(12);
    for(int i = 0;i<24;i++){
        vb::State onevibe{{det1,det2,det3,det4},{becky[myints[0]],becky[myints[1]],becky[myints[2]],becky[myints[3]]}};
        vb::SelectiveVB selectiveVB(oei,g,oi,{onevibe});
        std::cout<<std::endl<<"vb:"<<selectiveVB.solve()+repulsion<<std::endl;
        std::cout<<myints;
        std::next_permutation(myints,myints+4);

    }

    */
}
