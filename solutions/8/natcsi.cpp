#include <vector>
#include <cassert>

#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Sparse>

class NatCSI {
public:
    NatCSI(const std::vector<double> & t, const std::vector<double> & y)
        : t(t), y(y), h(t.size()-1), c(t.size()) {
        // Size check
        assert( ( t.size() == y.size() ) && "Error: mismatched size of t and y!");
        
        // m is the number of conditions (t goes from t_0 to t_n)
        // WARNING: m = n+1
        m = t.size();
        
        // Vector containing increments (from the right)
        for(int i = 0; i < (m - 1); ++i) {
            h(i) = t[i+1] - t[i];
            // Check that t is sorted
            assert( ( h(i) > 0 ) && "Error: array t must be sorted!");
        }
        
        // System matrix and rhs as in 3.5.9, we remove first and last row (known by natural contition)
        Eigen::SparseMatrix<double> A(m,m);
        Eigen::VectorXd b(m);
        
        // WARNING: sparse reserve space
        A.reserve(3);
        
        // Fill in natural conditions 3.5.10 for matrix
        A.coeffRef(0,0) = 2 / h(0);
        A.coeffRef(0,1) = 1 / h(0);
        A.coeffRef(m-1,m-2) = 1 / h(m-2);
        A.coeffRef(m-1,m-1) = 2 / h(m-2);
        
        // Reuse computation for rhs
        double bold = (y[1] - y[0]) / (h(0)*h(0));
        b(0) = 3*bold; // Fill in natural conditions 3.5.10
        // Fill matrix A and rhs b
        for(int i = 1; i < m-1; ++i) {
            // PRecompute b_i
            double hinv = 1./h(i);
            
            // Fill in a
            A.coeffRef(i,i-1) = hinv;
            A.coeffRef(i,i) = 2./h(i) + 2./h(i-1);
            A.coeffRef(i,i+1) = hinv;
            
            // Reuse computation for rhs b
            double bnew = (y[i+1] - y[i]) / (h(i)*h(i));
            b(i) = 3. * (bnew + bold);
            bold = bnew;
        }
        b(m-1) = 3*bold; // Fill in natural conditions 3.5.10
        // Compress the matrix
        A.makeCompressed();
        std::cout << A;
        
        // Factorize A and solve system A*c(1:end) = b
        Eigen::SparseLU<Eigen::SparseMatrix<double>> lu;
        lu.compute(A);
        c = lu.solve(b);
        
    }
    
    double operator() (double x) const {
        // Assert that x \in (t_0, t_n)
        assert( ( t[0] <= x && x <= t[m-1] ) && "Error: x must be in (t_0, t_n)");
        
        // Find j s.t. x \in (t_j,t_j+1)
        // as with the previous sheet (PS 6), we perform an efficient binary search on the data (log. complexity)
        auto j = (std::lower_bound(t.begin(), t.end(), x) - t.begin()) - 1;
        if( j == -1 ) j++;
        
        
        // Precompute tau and evaluate spline (3.5.5)
        double tau = (x - t[j]) / h(j);
        double tau2 = tau*tau;
        double tau3 = tau2*tau;
        return y[j] *         ( 1. - 3. * tau2 + 2 * tau3 ) +
               y[j+1] *           ( 3 * tau2 - 2 * tau3 ) +
               h(j) * c(j) *  ( tau - 2 * tau2 + tau3 ) +
               h(j) * c(j+1) *    ( - tau2 + tau3 );
               
    }
    
private:
    // Provided nodes and values to compute spline, same size Eigen vectors
    std::vector<double> t, y;
    // Difference t(i)-t(i-1) and coefficients of spline (s'(t_j))
    Eigen::VectorXd h, c;
    // Size of t, y and c  = m. h has size m-1
    int m;
};

int main() {
    
    int n = 8, m = 100;
    std::vector<double> t;
    std::vector<double> y;
    t.resize(n);
    y.resize(n);
    for(int i = 0; i < n; ++i) {
        t[i] = (double) i / (n-1);
        y[i] = cos(t[i]);
    }
    
    NatCSI N(t,y);
    
    Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(m, 0, 1);
    for(int i = 0; i < x.size(); ++i) {
        x(i) = N(x(i));
    }
    std::cout << x << std::endl;
}
