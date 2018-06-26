/*
 * Hard-coded driver for C4H6 VB calculations
 */

#include <libwint.hpp>
#include <iomanip>

int main() {

    std::cout << std::setprecision(12);  // Set print precision to 12

    //  Retrieve geomrty, caclulate integrals etc.
    std::string filename = "../../exe/C4H6.xyz";
    libwint::Molecule C4H6(filename);
    libwint::AOBasis ao_basis(C4H6, "STO-3G");  // init Atomic Orbital basis
    double repulsion = C4H6.calculateInternuclearRepulsionEnergy();  // repulsion
    ao_basis.calculateIntegrals();  // Calculate integrals
    int zen = 26;
    Eigen::Tensor<double,4> two_ei(zen,zen,zen,zen);  // two electron integrals for the custom.

    std::ofstream outfile("ggg");
    for(int i = 0;i<zen;i++){
        for(int j =0;j<zen;j++){
            for(int q =0;q<zen;q++){
                for(int k  = 0;k<zen;k++){
                    outfile<<two_ei(i,j,q,k);
                    outfile<<"\t";

                }
                outfile<<"\n";
            }
            outfile<<std::endl;
        }
        outfile<<std::endl;
    }

    std::cout<<two_ei(0,10,0,10)<<"wtf"<<std::endl;
    std::cout<<two_ei(0,10,10,10)<<"wtf"<<std::endl;
    std::cout<<two_ei(0,0,0,10)<<"wtf"<<std::endl;
    std::cout<<two_ei(10,10,0,0)<<"wtf"<<std::endl;
    std::cout<<two_ei(10,0,0,10)<<"wtf"<<std::endl;
    std::cout<<two_ei(10,1,0,10)<<"wtf"<<std::endl;

}


