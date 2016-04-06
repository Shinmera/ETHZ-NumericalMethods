#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>

double pi_approx(int k) {
    assert(k >= 1 && k < 8 && "k must be in {0,...,7}");
    
    // Store initial values of t = 1/n (nodes) and U (values at t)
    std::vector< double > t = {1/2., 1/3., 1/4., 1/5., 1/6., 1/8., 1/10.};
    std::vector< double > U = {
                                2.,
                                3/2.*std::sqrt(3),
                                2.*std::sqrt(2),
                                5/4.*std::sqrt(10 - 2.*std::sqrt(5)),
                                3.,
                                4.*std::sqrt(2 - sqrt(2)),
                                5/2.*(std::sqrt(5) - 1)
                            };
                            
    // Loop over each degree (tirangle scheme: left to right)
    for(int l = 0; l < k-1; ++l) {
        // Loop over each node (triangle scheme: top to bottom)
        for(int i = 0; i < k-l-1; ++i) {
            // Store the values at the beginning of U, forget about previous
            // value of U
            // t has to be offset by l
            U[i] = (-t[i]*U[i+1] + t[l+i+1]*U[i])/(t[l+i+1] - t[i]);
        }
    }
    // At the end U[0] contains the value of p at 0, i.e. almost \pi
    return U[0];
}

int main() {
    // Test routine
    for(int k = 1; k < 8; ++k) {
        std::cout << "k = " << k << ", pi_k = " << pi_approx(k) << std::endl;
    }
    
}