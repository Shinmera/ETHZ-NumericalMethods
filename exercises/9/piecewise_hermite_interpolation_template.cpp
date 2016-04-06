#include <vector>
#include <cassert>

#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Sparse>

// Flag for slope reconstruction type
enum class Slope { Zero, Reconstructed };

//! \breif Implements a piecewise cubic Herite interpolation on equidistant meshes.
class PCHI {
public:
    //! \brief Construct the slopes frome the data, either usinf dinite-differences or assuming f'(x_j) = 0
    //! \param[in] t vector of nodes (assumed equidistant and sorted)
    //! \param[in] y vector of values at nodes t
    //! \param[in] s Flag to set if you want to reconstruct or set slopes to zero
    PCHI(const Eigen::VectorXd & t, const Eigen::VectorXd & y, Slope s = Slope::Reconstructed)
        : t(t), y(y), c(t.size()) {
        // Sanity check
        n = t.size();
        assert( n == y.size() && "t and y must have same dimension." );
        assert( n >= 3 && "need at least two nodes." );
        h = t(1) - t(0);
        
        switch(s) {
            //// CASE: reconstruction of the slope, assuming f'(x_j) = 0 (O(1))
            case Slope::Zero:
                // TODO: reconstruct zero slope
                break;
            //// CASE: reconstruction of the slope using a second order finite difference (O(h^2))
            case Slope::Reconstructed:
            default:
                // TODO: reconstruct finite-difference slope
        }
    }
    
    //! \brief Evaluate the intepolant at the nodes x
    //! Input assumed sorted, unique and inside the interval
    //! \param[in] x vector of points t where to compute s(t)
    //! \return values of interpolant at x (vector)
    Eigen::VectorXd operator() (Eigen::VectorXd x) const {
        
        Eigen::VectorXd ret(x.size());
        // Stores the current interval index and some temporary variable
        
        // TODO: evaluate interpolant at x
        
        return ret;
    }
    
private:
    // Provided nodes and values (t,y) to compute spline, same size Eigen vectors, c contains slopes
    Eigen::VectorXd t, y, c;
    // Difference t(i)-t(i-1) and coefficients of spline (s'(t_j))
    double h;
    // Size of t, y and c.
    int n;
};

// Interpoland
auto f = [] (double x) { return 1. / (1. + x*x); };

int main() {
    
    double a = 5; // Interval (-a,a) bounds
    int M = 1000; // Number of  points in which to evaluate
    
    // Number of subintervals for each test
    std::vector<int> N = {4,8,16,32,64,128,256,512};
    
    // Precompute values at which evaluate f
    Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(M, -a, a);
    Eigen::VectorXd fx(x.size());
    for(int i = 0; i < x.size(); ++i) {
        fx(i) = f(x(i));
    }
    
    // TODO: error and rates and print
}
