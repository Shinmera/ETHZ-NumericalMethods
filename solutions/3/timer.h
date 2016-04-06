#include <chrono>

typedef std::chrono::high_resolution_clock::time_point time_t_;

//! Class implementing a simple chronometer, similar to tic-toc of Matlab
//! Just hit start() to start the timer, stop() to increment the timer (with time since last start/stop or reset)
//! Hit reset() to reset to 0 elapsed time
//! Get elapsed time with elapsed() or averaged time with avg()
template<typename duration_t_ = std::chrono::nanoseconds>
class timer {
public:
//     timer() {
//         elapsed_ = duration_t_::zero();
//         leaps_ = 0;
//     }

    //! Just save current time
    void start() {
        start_ = std::chrono::high_resolution_clock::now();
    }
    
    //! Increment timer by tme since last start, calls start()
    void stop() {
        auto tmp = std::chrono::duration_cast<duration_t_>( std::chrono::high_resolution_clock::now() - start_ );
        if(tmp < min_ or leaps_ == 0) min_ = tmp;
        elapsed_ += tmp;
        ++leaps_;
        start();
    }
    
    //! Reset timer to 0, calls start()
    void reset() {
        leaps_ = 0;
        elapsed_ = duration_t_::zero();
        start();
    }
    
    //! Get sum of all start()/stop() cycles
    duration_t_ elapsed() const {
        return elapsed_;
    }
    
    //! Get min of all start()/stop() cycles
    duration_t_ min() const {
        return min_;
    }
    
    //! Get average of all start()/stop() cycles
    duration_t_ avg() const {
        return elapsed_ / leaps_;
    }
    
    //! Print elapsed time
    void print() const {
        std::cout << elapsed().count() << std::endl;
    }
    
    //! Print averaged elapsed time
    void print_avg() const {
        std::cout << avg().count() << std::endl;
    }
    
    //! Print min elapsed time
    void print_min() const {
        std::cout << min().count() << std::endl;
    }
    
private:
    //! Number of cycles
    unsigned int    leaps_;
    //! Time of start()
    time_t_         start_;
    //! Elapsed until last stop() and minimal duration
    duration_t_     elapsed_, min_;

};
