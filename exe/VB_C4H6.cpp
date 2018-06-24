/*
 * Hard-code exe cause i'm lazy
 */

#include <libwint.hpp>
#include "common.hpp"
#include <hf.hpp>
#include <iomanip>
#include "SelectiveVB.hpp"

int main() {
    std::string filename = "../../exe/C4H6.xyz";
    libwint::Molecule C4H6 (filename);
    libwint::AOBasis ao_basis (C4H6, "STO-3G");
    double repulsion = C4H6.calculateInternuclearRepulsionEnergy();
    ao_basis.calculateIntegrals();
    hf::rhf::RHF rhf (C4H6, ao_basis, 1.0e-013);
    rhf.solve( hf::rhf::solver::SCFSolverType::DIIS);
    std::cout<<std::endl<<"RHF :"<<rhf.get_electronic_energy()+repulsion<<std::endl;
    libwint::SOBasis so_basis(ao_basis,rhf.get_C_canonical());

    Eigen::MatrixXd one_int_so = so_basis.get_h_SO();
    Eigen::Tensor<double,4> two_int_so = so_basis.get_g_SO();
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
    std::cout<<"good"<<energy;






    std::cout<<ao_basis.calculateNumberOfBasisFunctions()<<std::endl;
    int number = 17;
    Eigen::MatrixXd oei = Eigen::MatrixXd::Zero(number,number);
    Eigen::MatrixXd oi = Eigen::MatrixXd::Identity(number,number);
    Eigen::MatrixXd one = ao_basis.get_V()+ao_basis.get_T();
    Eigen::Tensor<double,4> g(number,number,number,number);
    g.setZero();
    int n2 = 13;
    int nr = 4;
    oei.block(0,0,n2,n2) = so_basis.get_h_SO().block(0,0,n2,n2);
    Eigen::Tensor<double,4> yyyy = so_basis.get_g_SO();
    for(int i = 0;i<n2;i++){
        for (int j = 0;j<n2;j++){
            for (int k =0;k<n2;k++){
                for(int l =0;l<n2;l++){
                    g(i,j,k,l) = yyyy(i,j,k,l);
                }
            }
        }
    }
    std::vector<int> pz = {10,15,20,25};
    Eigen::Tensor<double,4> xxx = ao_basis.get_g();
    for(int i = 0;i<nr;i++){
        for(int j = 0;j<nr;j++){
            oei(i+n2,j+n2) = one(pz[i],pz[j]);
            oi(i+n2,j+n2) = ao_basis.get_S()(pz[i],pz[j]);
            for (int k =0;k<nr;k++){
                for(int l =0;l<nr;l++){
                    g(i+n2,j+n2,k+n2,l+n2) = xxx(pz[i],pz[j],pz[k],pz[l]);


                }

            }

        }
    }
    Eigen::MatrixXd C = rhf.get_C_canonical();
    for(int i = 0;i<nr;i++){
        for(int j = 0; j<nr;j++){
            for(int k = 0;k<n2;k++){
                for (int l = 0;l<n2;l++){
                    for(int lo = 0;lo<26;lo++){
                        for(int zod = 0; zod<26;zod++){
                            g(i+n2,j+n2,k,l) += C(lo,k)*C(zod,l)*xxx(pz[i],pz[j],lo,zod);
                            g(l,k,j+n2,i+n2) =  g(i+n2,j+n2,k,l);

                        }
                    }
                }
            }
        }

    }



    //std::cout<<rhf.get_C_canonical();
    //std::cout<<oei;
    //std::cout<<std::endl;
    //std::cout<<oi;
    std::vector<int> baseline = {0,1,2,3,4,5,6,7,8,9,10,11,12};



    size_t det1 = vb::bitstring({0,1,2,3,4,5,6,7,8,9,10,11,12,15,16});
    size_t det2 = vb::bitstring({0,1,2,3,4,5,6,7,8,9,10,11,12,13,14});
    size_t det3 = vb::bitstring({0,1,2,3,4,5,6,7,8,9,10,11,12,14,15});
    size_t det4 = vb::bitstring({0,1,2,3,4,5,6,7,8,9,10,11,12,13,16});
//    size_t det1= vb::bitstring({0,1,2,3,4,5,6,7,8,9,10,11,12});
//    size_t det2= vb::bitstring({0,1,2,3,4,5,6,7,8,9,10,11,12});
//    size_t det3= vb::bitstring({0,1,2,3,4,5,6,7,8,9,10,11,12});
//    size_t det4= vb::bitstring({0,1,2,3,4,5,6,7,8,9,10,11,12});
    vb::State onevibe{{det1,det2,det3,det4},{det2,det1,det4,det3}};
    det1 = vb::bitstring({0,1,2,3,4,5,6,7,8,9,10,11,12,13,15});
    det2 = vb::bitstring({0,1,2,3,4,5,6,7,8,9,10,11,12,14,16});
    det3 = vb::bitstring({0,1,2,3,4,5,6,7,8,9,10,11,12,13,16});
    det4 = vb::bitstring({0,1,2,3,4,5,6,7,8,9,10,11,12,14,15});
    vb::State twovibe{{det1,det2,det3,det4},{det2,det1,det4,det3}};
    size_t det11 = vb::bitstring({15,16});
    size_t det22 = vb::bitstring({13,14});
    size_t det33 = vb::bitstring({14,15});
    size_t det44 = vb::bitstring({13,16});

    try{
        vb::SelectiveVB selectiveVB(oei,g,oi,{onevibe,twovibe});
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
