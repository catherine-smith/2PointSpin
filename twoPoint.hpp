/*
 *  file   twoPoint.hpp
 *
 *  Header file for twoPointSpin.cpp
 *  
 *  Includes Adaptable Binary Function Object for calculating
 *  matrix elements for two-point spin correlation values.
 *
 * !! Only collection of GroundState implemented !!
 *
 *  Catherine J. Stevenson 
 *  Jan 23 2008 - Modified to enable threading
 *  Mar 12 2008 - Symmetrized 2-Pt Matrix Element (ME)
 *              - 2-Pt units are now h-bar^2, not (h_bar/2)^2
 *                (incorporated the 1/4 into the ME)
 */

#ifndef GUARD_twoPoint_hpp
#define GUARD_twoPoint_hpp

#include <utility>
#include <complex>
#include <string>
#include <cmath>
#include <cstdlib>
#include "DiagonHeaders.h"
#include <boost/math/special_functions/laguerre.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/constants/constants.hpp>

double Pi();

double Root_Two_Pi();

/// Contains input file information
struct Input {
  int numFiles, input_Sz, numThreads;
  double input_omega0, input_B, r0, theta0, ETol, CTol, r_start, r_end, 
    r_step, theta_start, theta_end, arcLength;
  std::string output;
  std::vector<std::string> datFiles;
};

/// Reads information from input file
std::istream& Read_Input(std::istream&, Input&);

/// Write information from input file to stdout
void Write_Stdout(Input&);

/// Contains file information about eigenstate
struct Header {
  int N, Lz, tS, tSz, MPsize, numDegenStates, numEigs, numOmega0, numB;
  double Ez, omega0, B, Eg, ETol, CTol, sumOfSquares;
  std::string file;
};

typedef AntiSymmState<FDState, FDQNumLess> antiSymmState;

/// Stores Header info, energy and explicit eigenstate (as LinCombState object)
struct EigenstateInfo {
  Header h;
  double energy;
  LinCombState<antiSymmState, double> linCombState;
};

/// Finds eigenstates with specific omega0 and B values:
std::istream& Find_States(std::istream&, Input&, Header&,
			  std::vector<EigenstateInfo>&);

std::istream& next_omega0(std::istream&, int);

std::istream& next_B(std::istream&, int, int);

std::istream& antiSymmetrize(std::istream&, Header&, 
			     LinCombState<antiSymmState, double>&);

/// Sorts a vector of QMPairs by increasing (absolute value) coefficients:
bool compare(const QMPair<antiSymmState, double>&, 
	     const QMPair<antiSymmState, double>&);

/// Removes states from LinCombState whose coefficients are smaller than
/// a user-input tolerance (CTol)
void remove_States(double&, LinCombState<antiSymmState, double>&, double&);

/// Extracts degenerate eigenstates from input .dat files
std::istream& get_Eigenstates(std::istream&, Header&, std::vector<EigenstateInfo>&);

/// Fock-Darwin orbital wave function in real space. Calculates wave function
/// at user-supplied coordinates. Constructor takes FDState quantum numbers n
/// and m, and effective length scale lo (in units of length). Function takes
/// polar coordinates r (same units as lo) and theta. Output units are inverse
/// of the supplied units of lo.
class FDWaveFunction {
  int _n;
  int _m;
  double _lo; // units of length
public:
  FDWaveFunction(int n, int m, double lo): _n(n), _m(m), _lo(lo) { }
  // In function below, r must have same units as lo:
  std::complex<double> operator()(double& r, double& theta) { 

    using namespace boost::math;
    unsigned nPrime = _n + _m;
    int mPrime = _n - _m;
    unsigned absMPrime = mPrime < 0 ? -mPrime : mPrime;

    // Angular part:
    std::complex<double> I(0, 1.0);
    std::complex<double> mPrimeComplex(mPrime, 0);
    std::complex<double> thetaComplex(theta, 0);
    std::complex<double> AngularComponent;
    std::complex<double> Root2Pi(Root_Two_Pi(), 0);
    AngularComponent = exp(I*mPrimeComplex*thetaComplex)/Root2Pi;
    
    // Radial part:
    unsigned nr = (nPrime - absMPrime)/2;
    double Root = sqrt(tgamma(nr+1)/tgamma(nr+absMPrime+1));
    //    double Root = sqrt(factorial(nr)/factorial(nr+absMPrime));
    double arg = pow(r,2)/(2*pow(_lo,2));
    double LaguerrePoly = laguerre(nr, absMPrime, arg);
    double Radial = (Root*exp(-arg/2)*pow(sqrt(arg), fabs(mPrime))
		     *LaguerrePoly)/(_lo);
    std::complex<double> RadialComponent(Radial, 0);
    
    return pow(-1,nr)*AngularComponent*RadialComponent;
  } 
};

/// S+z(r0)S+/-z(r1): Symmetrized Two-Point Correlation Matrix Element
/// This is an Adaptable Binary Function Object.
/// Input are two FDstates, two sets of polar coordinates (one, r0,
/// held constant) and an effective length scale (lo). Template is
/// instantiated with int of value 1 for calculating <S+z(r0)S+z(r1)> 
/// or value -1 for calculating <S+z(r0)S-z(r1)>. 
/// Resulting Matrix Element is in units of h-bar^2*(inverse-
/// length scale units)^4.
template <int Sz>
class SzSzME: public Arity4Function<FDState, FDState, FDState, FDState, 
				    std::complex<double> > {
  double _r0;
  double _theta0;
  double _r1;
  double _theta1;
  double _lo;
public:
  // Constructor takes fixed coordinates r0 and theta0, and "probe" coord
  // r1 and theta1, and effective length lo:
  SzSzME(double r0, double theta0, double r1, double theta1, double lo):
    _r0(r0), _theta0(theta0), _r1(r1), _theta1(theta1), _lo(lo) {
    if (std::abs(Sz) != 1) {
      std::cout << "Spin projection incompatible with spin matrix element. Aborting.\n";
      abort();      
    }
  }
  
  // Public access to r0, theta0, r1, theta1 used to construct SzSzME:
  double r0() const { return _r0; }
  double theta0() const { return _theta0; }
  double r1() const { return _r1; }
  double theta1() const { return _theta1; }

  // Define <SzSzME> matrix element:
  std::complex<double> operator()(const FDState& bra1, const FDState& bra2,
				  const FDState& ket1, const FDState& ket2) {
    // This operator conserves spin:
    if (bra1.spin() != ket1.spin() || bra2.spin() != ket2.spin())
      return 0;
//     // We are calculating S+z(r0)S+z/S-z(r1), therefore both bra1.spin()
//     // and ket1.spin() must be up. Also, bra/ket2.spin() must be Sz:
//     if ( (ket1.spin() != 1) || (ket2.spin() != Sz) )
//       return 0;

    // Initialize wavefunctions from bra and ket quantum numbers:
    FDWaveFunction PhiA(bra1.n(), bra1.m(), _lo);
    FDWaveFunction PhiB(bra2.n(), bra2.m(), _lo);
    FDWaveFunction PhiC(ket1.n(), ket1.m(), _lo);
    FDWaveFunction PhiD(ket2.n(), ket2.m(), _lo);

    // Symmetric ME has two terms:
    std::complex<double> ME1 = 0;
    std::complex<double> ME2 = 0;

    //Evaluate each term at appropriate r and theta:
    if ( (ket1.spin() == 1) && (ket2.spin() == Sz) ) 
      {
	std::complex<double> phi1 = PhiA(_r0, _theta0);
	std::complex<double> phi2 = PhiB(_r1, _theta1);
	std::complex<double> phi3 = PhiC(_r0, _theta0);
	std::complex<double> phi4 = PhiD(_r1, _theta1);
	
	ME1 = 0.25*conj(phi1)*conj(phi2)*phi3*phi4;
      } 
    
    if ( (ket1.spin() == Sz) && (ket2.spin() == 1) )
      {
	std::complex<double> phi5 = PhiB(_r0, _theta0);
	std::complex<double> phi6 = PhiA(_r1, _theta1);
	std::complex<double> phi7 = PhiD(_r0, _theta0);
	std::complex<double> phi8 = PhiC(_r1, _theta1);
	
	ME2 = 0.25*conj(phi5)*conj(phi6)*phi7*phi8;
      }
    
    return (ME1 + ME2);      
  }
};

/// Calculates two-point correlation values (<S+z(r0)S+z(r1)> or <S+z(r0)S-z(r1)>)
/// at a specific spatial location (r, theta). 
class twoPointCalc {
  std::complex<double> _plusCorrelation;
  std::complex<double> _minusCorrelation;
  std::vector<EigenstateInfo>* _eigs_p; // pointer to vec of EigenstateInfo objects
  SzSzME<1> _MEPlus;
  SzSzME<-1> _MEMinus;
public:
  // Constructor takes address of vector of EigenstateInfo objects, two matrix 
  // elements, SzSzME<1> and SzSzME<-1>, instantiated at the same coordinates, and
  // input Sz spin for desired calculation (<S+z(r0)S+z(r1)> or <S+z(r0)S-z(r1)>):
  twoPointCalc(std::vector<EigenstateInfo>* const eigs_p, SzSzME<1> MEPlus,
	   SzSzME<-1> MEMinus): _plusCorrelation(0), _minusCorrelation(0), 
				_eigs_p(eigs_p), _MEPlus(MEPlus),
				_MEMinus(MEMinus) {  
    // Check that both MEPlus and MEMinus are instantiated at the same coords:
    if (_MEPlus.r0() != _MEMinus.r0() || _MEPlus.theta0() != _MEMinus.theta0() ||
	_MEPlus.r1() != _MEMinus.r1() || _MEPlus.theta1() != _MEMinus.theta1()) {
      std::cout << "Matrix element instantiation corrupted. Aborting. \n";
      abort();
    }
  }

  //Public access to correlation, coordinates: Since the coords for each _MEPlus
  // and _MEMinus are equal, it doesn't matter which ones we return:
  double plusCorrelation_real() const { return _plusCorrelation.real(); }
  double plusCorrelation_imag() const { return _plusCorrelation.imag(); }
  double minusCorrelation_real() const { return _minusCorrelation.real(); }
  double minusCorrelation_imag() const { return _minusCorrelation.imag(); }
  double r0() const { return _MEPlus.r0(); }
  double theta0() const { return _MEPlus.theta0(); }
  double r1() const { return _MEPlus.r1(); }
  double theta1() const { return _MEPlus.theta1(); }

  void operator() () {
    for (size_t i = 0; i != (*_eigs_p).size(); ++i) {
      _plusCorrelation += twoBodySOp((*_eigs_p)[i].linCombState,
				 (*_eigs_p)[i].linCombState, _MEPlus);
      _minusCorrelation += twoBodySOp((*_eigs_p)[i].linCombState,
				 (*_eigs_p)[i].linCombState, _MEMinus);
    }
    
    // The final correlation value is the average of the correlations of each
    // denegerate eigenstate:
    _plusCorrelation.real() = _plusCorrelation.real()/(*_eigs_p).size();
    _plusCorrelation.imag() = _plusCorrelation.imag()/(*_eigs_p).size();
    _minusCorrelation.real() = _minusCorrelation.real()/(*_eigs_p).size();
    _minusCorrelation.imag() = _minusCorrelation.imag()/(*_eigs_p).size();
    
    return;
  }

    // Writing to stdout, back-up file can be coded later. Will require mutexes.
};


#endif
