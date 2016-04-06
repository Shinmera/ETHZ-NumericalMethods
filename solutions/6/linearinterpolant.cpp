#include <vector>
#include <algorithm>
#include <iostream>
#include <cassert>

class LinearInterpolant {
    public:
        //! pair holds the pair (t_i, y_i) with the value of the interpolant y_i at the node t_i
        using pair = std::pair<double, double>;
        //! data holds the list of pairs, may be ordered or not, but *cannot* contain duplicate t_i's
        using data = std::vector< pair >;

        //! Constructor, builds the (t_i,y_i) Pairs from the data provided, i.e. defines I = \sum c_i b_i
        //! Sort the array for the first time, the data is not assumed to be sorted, sorting is necessary for binary search
        LinearInterpolant(data && Pairs_) {
            // WARNING: assignment of std::vectors in C++ performs deep copy
            Pairs = Pairs_;
            assert( Pairs.size() > 1 && "Must specify at least two nodes");
            
            // This is needed in std::sort
            auto Order = [] (const pair & P, const pair & Q) -> bool {
                return P.first < Q.first;
            };
            // Sort the array
            std::sort(Pairs.begin(), Pairs.end(), Order);
        }
        
        //! Evaluation operator, value of I at x, i.e. I(x).
        //! Performs bound checks (i.e. if x < t_0 or x >= t_n )
        double operator() (double x) {
            
            // Lambdas for compasiron of pair type, we want to compare only the 
            // first element of the pair, i.e. the t_i's
            auto Compare = [] (const pair & P, double V) -> bool { return P.first < V; };
            
            // Find the place i (as iterator), wehre t_i >= x, notice that Pairs 
            // *must* be sorted (This is needed in std::lower_bound)
            // IMPORTANT: lower_bound performs a binary_search on the data, provided the data has a random access iterator 
            // http://www.cplusplus.com/reference/iterator/RandomAccessIterator/$
            // this is crucial here, since evaluation operators must be fast
            // this allows to reduce the complexity, from a linear complexity to a logarithmic one
            auto it = std::lower_bound( Pairs.begin(), Pairs.end(), x, Compare );
            // Bound checks, if i = 0 and t_0 != x (we are before the first node) or if x > t_n (we are after the last node)
            // In such cases return 0 (*it and *(it-1) would be undefined in such cases)
            if( ( it == Pairs.begin() && it->first != x ) || it == Pairs.end() ) return 0;
            
            // Actually compute the interpolated value, dist_rato contains the distance from t_{i-1} to x (as a ratio)
            double dist_ratio = (x - (it-1)->first) / (it->first - (it-1)->first);
            return (it-1)->second * (1 - dist_ratio) + it->second * dist_ratio; 
        }
    private:
        data    Pairs; //! vector of pairs (t_i, y_i), nodes and value at nodes
};

int main(void) {
    // Test the class with the basis with nodes (-1,1,2,4) and interpolant with values (-1,2,3,4) at said nodes.
    LinearInterpolant I = LinearInterpolant({{1,2},{2,3},{4,4},{-1,-1}});
    // Output value at the specified point
    std::cout << I(-2) << " " << I(-1) << " " << I(1) << " " << I(2) << " " << I(4) << " " << I(5) << std::endl;
    std::cout << I(0.5) << " " << I(1) << " " << I(1.5) << " " << I(2.1) << " " << I(3) << " " << I(3.1) << " " << I(4) << std::endl;
}
