/*==========================================================
 * tCGGPsumrnd.cpp
 *
 * The calling syntax is:
 *
 *		S_w = tCGGPsumrnd(p, logalpha, sigma, tau, a, b, gamma, T)
 *		S_w = tCGGPsumrnd(p, logalpha, sigma, tau, a, b, gamma, T, seed)
 *
 * Compilation command (from MATLAB console)
 * -----------------------------------------
 * Using gcc (Linux, Mac):
 * - using a C++11 compiler (gcc>=4.3):
 *
 *      mex tCGGPsumrnd.cpp CXXFLAGS='$CXXFLAGS -std=c++0x -fopenmp' LDFLAGS='$LDFLAGS -fopenmp'
 *
 * - using Boost installed in a standard path:
 *
 *      mex tCGGPsumrnd.cpp CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp'
 *
 * - using custom Boost path:
 *
 *      mex tCGGPsumrnd.cpp -Ipath/to/boost CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp'
 *
 * Using Microsoft Visual C++ (Windows) and Boost:
 * 
 *      mex tCGGPsumrnd.cpp -Ipath/to/boost COMPFLAGS="$COMPFLAGS /openmp" LINKFLAGS="$LINKFLAGS /openmp"
 *
 * e.g: mex tCGGPsumrnd.cpp -I"C:/Program Files/boost/boost_1_58_0" COMPFLAGS="$COMPFLAGS /openmp" LINKFLAGS="$LINKFLAGS /openmp"
 *
 *========================================================*/

#include "mex.h"
#include <cmath>
#include <ctime>
#include <limits>
#include <omp.h>

typedef double double_;

#if __cplusplus != 199711L /* C++11 std */
        
#include <random> 
typedef std::default_random_engine generator; 
typedef std::uniform_real_distribution<double_> uniform_dist; 
typedef std::gamma_distribution<double_> gamma_dist; 

#else /* Boost */

#include <boost/math/special_functions/gamma.hpp>
using boost::math::lgamma;

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/gamma_distribution.hpp>
typedef boost::random::mt19937 generator;
typedef boost::random::uniform_real_distribution<double_> uniform_dist;
typedef boost::random::gamma_distribution<double_> gamma_dist;

#endif

void tCGGPsumrnd( const unsigned p,  const double logalpha, const double sigma,
        const double tau, const double * a, const double * b, const double * gamma, 
        double_ T, const double_ Tmax, const unsigned seed,
        double_ * S_w );

void tCGGPsumrnd_tau0( const unsigned p,  const double logalpha, const double sigma,
        const double * a, const double * b, const double * gamma, 
        double_ T, const double_ Tmax, const unsigned seed,
        double_ * S_w );

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    /* check for proper number of arguments */
    if(nrhs<8) {
        mexErrMsgTxt("At least 8 inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgTxt("One output required.");
    }
    
    /* get the value of the input  */
    const unsigned p = mxGetScalar( prhs[0] );
    const double logalpha = mxGetScalar( prhs[1] );
    const double sigma = mxGetScalar( prhs[2] );
    const double tau = mxGetScalar( prhs[3] );
    double * a = mxGetPr( prhs[4] );
    const unsigned na = mxGetNumberOfElements( prhs[4] );
    double * b = mxGetPr( prhs[5] );
    const unsigned nb = mxGetNumberOfElements( prhs[5] );
    double * gamma = mxGetPr( prhs[6] );
    const unsigned ng = mxGetNumberOfElements( prhs[6] );
    const double T = mxGetScalar( prhs[7] );
    double seed = time(NULL);
    if ( nrhs >= 9 )
        seed = mxGetScalar( prhs[8] );
    // number of concurent jobs
    unsigned nthreads = omp_get_max_threads();
    if ( nrhs >= 10 )
        nthreads = unsigned(mxGetScalar( prhs[9] ));
    
    if ( na < p && na==1 ) {
        // duplicate value of a
        double a_val = *a;
        a = new double[p];
        for ( unsigned k = 0; k < p; ++k )
            a[k] = a_val;
    }
    if ( nb < p && nb==1 ) {
        // duplicate value of b
        double b_val = *a;
        b = new double[p];
        for ( unsigned k = 0; k < p; ++k )
            b[k] = b_val;
    }
    if ( ng < p && ng==1 ) {
        // duplicate value of gamma
        double g_val = *a;
        gamma = new double[p];
        for ( unsigned k = 0; k < p; ++k )
            gamma[k] = g_val;
    }
    
    // expected log number of jumps
    const double_ cst = logalpha - lgamma( 1.0 - sigma ) - log(sigma);
    const double_ log_njumps = cst - sigma * log(T);
    
    /* cut into nthreads intervals [Tcut(i), Tcut(i+1)] such that the expected
     * number of jumps is njumps/nthreads in each interval */
    double_ * Tcut = new double_[nthreads+1];
    Tcut[0] = T;
    for (unsigned i = 1; i < nthreads; ++i) {
        // TODO: find a better bound of the expected number of jumps
        Tcut[i] = pow( pow(Tcut[i-1], -sigma) - exp( log_njumps - log(double(nthreads)) - cst ), -1.0 / sigma);
    }
    Tcut[nthreads] = std::numeric_limits<double_>::infinity();
            
    double_ ** S_w = new double_*[nthreads];
    for (unsigned i = 0; i < nthreads; ++i)
        S_w[i] = new double_[p];
    
    /* call the computational routine */
    if ( tau < 1e-8 ) {
        // case tau==0
        #pragma omp parallel for
        for (int i = 0; i < nthreads; ++i) {
            tCGGPsumrnd_tau0( p, logalpha, sigma, a, b, gamma, Tcut[i], Tcut[i+1], unsigned(seed)+i,
                S_w[i] );
        }
    } else {
        // case tau>0
        #pragma omp parallel for
        for (int i = 0; i < nthreads; ++i) {
            tCGGPsumrnd( p, logalpha, sigma, tau, a, b, gamma, Tcut[i], Tcut[i+1], unsigned(seed)+i,
                S_w[i] );
        }
    }
    
    /* reduction */
    for (unsigned i = 1; i < nthreads; ++i) {
        for (unsigned k = 0; k < p; ++k)
            S_w[0][k] += S_w[i][k];
    }
    
    /* create the output */
    plhs[0] = mxCreateDoubleMatrix( 1, p, mxREAL );
    for (unsigned k = 0; k < p; ++k)
        mxGetPr( plhs[0] )[k] = S_w[0][k];
    
    /* deallocate memory */
    delete[] Tcut;
    for (unsigned i = 0; i < nthreads; ++i)
        delete[] S_w[i];
    delete[] S_w;
    if ( na < p && na==1 )
        delete[] a;
    if ( nb < p && nb==1 )
        delete[] b;
    if ( ng < p && ng==1 )
        delete[] gamma;
}

/* The computational routine */

void tCGGPsumrnd( const unsigned p,  const double logalpha, const double sigma,
        const double tau, const double * a, const double * b, const double * gamma, 
        double_ T, const double_ Tmax, const unsigned seed,
        double_ * S_w )
{
    const double_ sigmap1 = sigma + 1.0;
    const double_ tauinv = 1.0 / tau;
    const double_ log_cst = logalpha - lgamma( 1.0 - sigma ) - log(tau);
    double_ * gamma_b = new double_[p];
    for ( unsigned k = 0; k < p; ++k )
        gamma_b[k] = gamma[k]/b[k];
    double_ log_r;    
    double_ log_mgf;       
    double_ t_new;    
    double_ log_mgf_new; 
    double_ log_G;
    generator gen(seed);
    uniform_dist unif( 0.0, 1.0 );
    gamma_dist gam;

    for ( unsigned k = 0; k < p; ++k )
        S_w[k] = 0.0;
       
    while (true) {
        log_r = log( -log( unif(gen) ) );
        log_mgf = 0.0;
        for ( unsigned k = 0; k < p; ++k )
            log_mgf -= a[k] * log1p( T * gamma_b[k] );
        log_G = log_mgf + log_cst - sigmap1 * log(T) - tau * T;
        if ( log_r > log_G )
            break;
        t_new = T - tauinv * log1p( -exp( log_r - log_G ) );
        if ( t_new > Tmax )
            break;
        log_mgf_new = 0.0;
        for ( unsigned k = 0; k < p; ++k )
            log_mgf_new -= a[k] * log1p( t_new * gamma_b[k] );
        if ( log( unif(gen) ) < sigmap1 * ( log(T) - log(t_new) ) + log_mgf_new - log_mgf ) {
            for ( unsigned k = 0; k < p; ++k ) {
                gam = gamma_dist( a[k], 1 / ( b[k] + t_new * gamma[k] ) );
                S_w[k] += exp( log(t_new) + log( gam(gen) ) );
                //gam = gamma_dist( a[k], 1 / ( b[k] / t_new + gamma[k] ) );
                //S_w[k] += gam(gen);
            }
        }
        T = t_new;
    }
    
    delete[] gamma_b;
}


void tCGGPsumrnd_tau0( const unsigned p,  const double logalpha, const double sigma,
        const double * a, const double * b, const double * gamma, 
        double_ T, const double_ Tmax, const unsigned seed,
        double_ * S_w )
{
    const double_ log_cst = logalpha - lgamma( 1.0 - sigma ) - log(sigma);
    const double_ msigma = -sigma;
    const double_ msigmainv = 1.0 / msigma;
    double_ * gamma_b = new double_[p];
    for ( unsigned k = 0; k < p; ++k )
        gamma_b[k] = gamma[k]/b[k];
    double_ log_r;    
    double_ log_mgf;       
    double_ t_new;    
    double_ log_mgf_new; 
    double_ log_temp; 
    double_ log_tmsigma;
    generator gen(seed);
    uniform_dist unif( 0.0, 1.0 );
    gamma_dist gam;

    for ( unsigned k = 0; k < p; ++k )
        S_w[k] = 0.0;
       
    while (true) {
        log_r = log( -log( unif(gen) ) );
        log_mgf = 0.0;
        for ( unsigned k = 0; k < p; ++k )
            log_mgf -= a[k] * log1p( T * gamma_b[k] );
        log_temp = log_r - log_mgf - log_cst;
        log_tmsigma = msigma * log(T);
        if ( log_temp > log_tmsigma )
            break;
        t_new = pow( exp(log_tmsigma) - exp(log_temp) , msigmainv );
        if ( t_new > Tmax )
            break;
        log_mgf_new = 0.0;
        for ( unsigned k = 0; k < p; ++k )
            log_mgf_new -= a[k] * log1p( t_new * gamma_b[k] );
        if ( log( unif(gen) ) < log_mgf_new - log_mgf ) {
            for ( unsigned k = 0; k < p; ++k ) {
                gam = gamma_dist( a[k], 1 / ( b[k] / t_new + gamma[k] ) );
                S_w[k] += gam(gen);
            }
        }
        T = t_new;
    }
    
    delete[] gamma_b;
}
