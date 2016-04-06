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
                c = Eigen::VectorXd::Zero(n);
                break;
            //// CASE: reconstruction of the slope using a second order finite difference (O(h^2))
            case Slope::Reconstructed:
            default:
//                 c(0) = ( y(1) - y(0) ) / h; // First order
                c(0) = ( -1*y(2) + 4*y(1) - 3*y(0) ) / 2 / h;
                for(int i = 1; i < n-1; ++i) {
                    c(i) = ( y(i+1) - y(i-1) ) / 2 / h;
                }
//                 c(n-1) = ( y(n-1) - y(n-2) ) / h; // First order
                c(n-1) = ( 3*y(n-1) - 4*y(n-2) + 1*y(n-3) ) / 2 / h;
                break;
        }
    }
    
    //! \brief Evaluate the intepolant at the nodes x
    //! Input assumed sorted, unique and inside the interval
    //! \param[in] x vector of points t where to compute s(t)
    //! \return values of interpolant at x (vector)
    Eigen::VectorXd operator() (Eigen::VectorXd x) const {
        
        Eigen::VectorXd ret(x.size());
        // Stores the current interval index and some temporary variable
        int i_star = 0;
        double tmp,t1,y1,y2,c1,c2,a1,a2,a3;
        for(int j = 0; j < x.size(); ++j) {
            // Find interval and porting of hermloceval Matlab code 3.4.6
            if( t(i_star) < x(j) || i_star == 0) { 
                ++i_star;
                t1 = t(i_star-1);
                y1 = y(i_star-1);
                y2 = y(i_star);
                c1 = c(i_star-1);
                c2 = c(i_star);
                a1 = y2 - y1;
                a2 = a1 - h*c1;
                a3 = h*c2 - a1 - a2;
            }
            // Compute s(x(j))
            tmp = ( x(j) - t1 ) / h;
            ret(j) = y1 + ( a1 + ( a2+a3*tmp ) * ( tmp - 1. ) ) * tmp;
        }
        
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
// auto f = [] (double x) {return cos(x); };

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
    
    // Store error and rates
    std::vector<double> err, err_zero, rate, rate_zero;
    
    std::cout << "L^infty-error [reconstruction, zero] (rate / rate)" << std::endl;
    for(int n : N) {
        // Define subintervals and evaluate f there (find pairs (t,y))
        Eigen::VectorXd t = Eigen::VectorXd::LinSpaced(n, -a, a);
        Eigen::VectorXd y(t.size());
        for(int i = 0; i < t.size(); ++i) {
            y(i) = f(t(i));
        }
        
        // Construct PCHI with zero and reconstructed slopes
        PCHI P(t,y), Pz(t,y,Slope::Zero);
        
        // Compute infinity norm of error
        err.push_back((P(x) - fx).lpNorm<Eigen::Infinity>());
        err_zero.push_back((Pz(x) - fx).lpNorm<Eigen::Infinity>());
        
        // Store errors and rates
        std::cout << err.back()  << "    " << err_zero.back() ;
        if( err.size() > 1 ) {
            rate.push_back( log( *(err.end() - 2) / err.back()  ) / log(2) );
            rate_zero.push_back( log( *(err_zero.end() - 2) / err_zero.back() ) / log(2) );
            std::cout << " (" << rate.back() << " / " << rate_zero.back() << ")";
        }
        std::cout << std::endl;
        
    }
}
