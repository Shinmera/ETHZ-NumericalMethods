#pragma once // header guard

#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>

#include <cmath>

//! Shortcut for container vector of nodes/weights
using vector = Eigen::VectorXd;

//! Structure containing a Quadrature rule on [-1,1], comprised of weights and nodes
struct QuadRule {
    vector weights;
    vector nodes;
};

//! \brief Golub-Welsh implementation 5.3.35
//! \param[in] n number of Gauss nodes
//! \return structure QuadRule containing nodes and wights of gauss Quadrature
QuadRule gauleg(unsigned int n);
