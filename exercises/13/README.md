# 1. Usage #

ode45 is header-only, just include and enjoy. defines class ode45 implementing emulation of Matlab's RK45 ode45.

Instructions:
1. construct ode45 class instance
    Create an instance of the class, pass r.h.s. function f to constructor:
        ode45<StateType> = O(f);
    
    Template parameters:
        StateType:  type of initial data and solution (state space)
        RhsType:    type of rhs function (automatically deduced)
        
    Arguments:
        f:          rhs function handle with operator()(const StateType & vec) -> StateType
        
2. (optional) Set options
    Set members of struct ode45.options to configure the solver:
        O.options.<option_you_want_to_set> = <value>
    examples:
        rtol:       relative tolerance for error control (default = 10e-6)
        atol:       absolute tolerance for error control (default = 10e-8)
    e.g.:
        O.options.rtol = 10e-5;
        
    ode45::Options struct public members (set them to wanted value):
        bool         save_init         = true;  //!< Set true if you want to save the initial data
        bool         fixed_stepsize    = false; //!< Set true if you want a fixed step size
        unsigned int max_iterations    = 5000;  //!< Set the maximum number of rejected iterations
        double       min_dt            = -1.;   //!< Set the minimum step size
        double       max_dt            = -1.;   //!< Set the maximum step size
        double       initial_dt        = -1.;   //!< Set an initial step size
        double       start_time        = 0;     //!< Set a starting time
        double       atol              = 10e-8;
        double       rtol              = 10e-6;
        
3. Solve
    Call the solver:
        std::vectro<std::pair<StateType> sol = O.solve(y0, T, norm)
        
    Template parameters:
        NormType:   type of norm function, automatically deduced
        
    Arguments:
        y0:         initial value in StateType (y(0) = y0)
        T:          final time of integration
        norm:       (optionaL) norm function to call on member of StateType, for computation of error
    
    Return:
        ODE solution std::vector of std::pair (y(t), t) for every snapshot
        
4. Get statistics
    Member statistics of ode45 of type Statistics contains some useful statistic for the solver
    
5. Print
    Print additional info ad the end, use O.print()

# 2. Sources/References: #
 - Emulates http://ch.mathworks.com/help/matlab/ref/ode45.html
 - Ported from ode45.py: https://github.com/rngantner/
 - See also: http://sourceforge.net/p/octave/odepkg/

# 3. Changelog: #
 - initial release
 
# 4. TODOs: #
 - non-autonomous ODEs
 - backward time integration
 - fixed time step
 - selective snapshot saving
 - componentwise norm control
 
# 5. License: #
  GNU GPL v.3
