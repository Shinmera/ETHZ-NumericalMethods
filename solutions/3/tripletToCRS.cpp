#include <iostream>
#include <vector>
#include <algorithm>

#include <Eigen/Dense>

// Uncomment if you want to use this structure for TripletMatrix
template <class scalar>
struct Triplet {
    
};

template <class scalar>
struct TripletMatrix {
    // TODO: insert here members and methods to TripletMatrix
    
    // densify() prototype
    Eigen::Matrix<scalar, -1, -1> densify() const;
};

// Uncomment if you want to use this structure for CRSMatrix
// template <class scalar>
// struct ColValPair {
//     
// };

template <class scalar>
struct CRSMatrix {
    // TODO: insert here members and methods to CRSMatrix
    
    // densify() prototype
    Eigen::Matrix<scalar, -1, -1> densify() const;
};

template <class scalar>
void tripletToCRS(const TripletMatrix<scalar> & T, CRSMatrix<scalar> & C) {
    // TODO: conversion function
}

//! \brief overload of operator << for output of Triplet Matrix (debug).
//! WARNING: uses densify() so there may be a lot of fill-in
//! this allows something like std::cout << S
//! \param o standard output stream
//! \param S matrix in Triplet matrix format
//! \return a ostream o, s.t. you can write o << A << B;
std::ostream & operator<<(std::ostream & o, const TripletMatrix<double> & S) {
    return o << S.densify();
}

//! \brief overload of operator << for output of CRS Matrix (debug).
//! WARNING: uses densify() so there may be a lot of fill-in
//! this allows something like std::cout << S
//! \param o standard output stream
//! \param S matrix in CRS matrix format
//! \return a ostream o, s.t. you can write o << A << B;
std::ostream & operator<<(std::ostream & o, const CRSMatrix<double> & S) {
    return o << S.densify();
}

int main() {
    //// Correctness test
    std::size_t nrows = 7, ncols = 5, ntriplets = 9;
    
    TripletMatrix<double> T;
    CRSMatrix<double> C;
    
    // TODO: contrtuct T here
    // TODO: Use this loop to push back triplets in your matrix
    for(auto i = 0u; i < ntriplets; ++i) {
        // TODO: Insert triplet (rand() % nrows, rand() % ncols, rand() % 1000))
    }
    
    std::cout << "***Test conversion with random matrices***" << std::endl;
    tripletToCRS(T, C);
    // TODO: Uncomment if you implemented densify()
//     std::cout << "--> Frobenius norm of T - C: " << (T.densify()-C.densify()).norm() << std::endl;
    std::cout << "T = " << std::endl << T << std::endl;
    std::cout << "C = " << std::endl << C << std::endl;
}
