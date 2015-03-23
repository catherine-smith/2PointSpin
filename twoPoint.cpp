/*
 *  file   twoPoint.cpp
 *  
 *  See header file for description
 *
 *  Catherine J. Stevenson 
 *  Jan 23, 2008 - Modified to enable threading
 *
 */

#include "twoPoint.hpp"
#include "DiagonHeaders.h"
#include "CloseAbs.hpp"
#include <boost/math/special_functions/laguerre.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/constants/constants.hpp>
#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <complex>
#include <ios>
#include <iomanip>

using std::fabs;    using std::istream;  using std::vector;
using std::complex; using std::string;   using std::cout;
using std::endl;    using std::ostream;  using std::streamsize;
using std::setprecision;

double Pi() {
  return boost::math::constants::pi<double>();
}

double Root_Two_Pi() {
  return boost::math::constants::root_two_pi<double>();
}

// FUNCTION TO READ INPUT FILE:
istream& Read_Input(istream& in, Input& input) {
  in >> input.numFiles;
  string junk;
  getline(in, junk); // skip remaining text in line
  // Get .dat file names:
  for (int i = 0; i != input.numFiles; ++i) {
    string file;
    getline(in, file);
    input.datFiles.push_back(file);
  }
  // Get the number of threads for calculation:
  in >> input.numThreads;
  getline(in, junk);
  // Get omega0 and B (Only states with these values will be examined):
  in >> input.input_omega0 >> input.input_B;
  getline(in, junk);
  // Get reference coordinates for two-point correlation,
  // convert theta from degrees into rad:
  double theta0;
  in >> input.r0 >> theta0;
  getline(in, junk);
  input.theta0 = theta0*Pi()/180;
  // Get energy tolerance for comparing .dat file GS energies:
  in >> input.ETol;
  getline(in, junk);
  // Get coefficient tolerance for calculation:
  in >> input.CTol;
  getline(in, junk);
  // Get starting radius, end radius, and radius step size:
  in >> input.r_start >> input.r_end >> input.r_step;
  getline(in, junk);
  // Get starting angle, end angle:
  double theta_start, theta_end;
  in >> theta_start >> theta_end;
  getline(in, junk);
  // convert into radian and store in Input structure:
  input.theta_start = theta_start*Pi()/180;
  input.theta_end = theta_end*Pi()/180;
  // Get arc length: theta step size (for iterating through theta) will be
  // (arc length)/r:
  in >> input.arcLength; // in units of Pi
  getline(in, junk);
  // Get name of output file:
  in >> input.output;
  getline(in, junk);
  
  return in;
}

// FUNCTION TO WRITE INPUT DATA TO STDOUT:
void Write_Stdout(Input& input) {

  // Machine's default precision (for resetting precision later)
  streamsize prec = cout.precision();

  // Number of files used:
  cout << "Files read (total " << input.numFiles << "): " << endl;
  // Name of files:
  for (size_t i = 0; i != input.datFiles.size(); ++i) {
    cout << input.datFiles[i] << endl;
  }
  // Reference coordinate:
  cout << "Reference coordinate r0 = (" << input.r0 << ", "
       << setprecision(4) << input.theta0*180/Pi() << ")" << endl;
  // Omega, B:
  cout << '\n' << "State examined has omega0 = " << input.input_omega0
      << " meV, B = " << input.input_B << " T" << endl;
  // Energy Tolerance, Coefficient tolerance:
  cout << "Energy equality tolerance: " << input.ETol << " (meV)" << endl;
  cout << "Coefficient tolerance: " << input.CTol << endl;
  // Area calculated:
  cout << '\n' << "r_start: " << input.r_start << '\n' << "r_end: "
       << input.r_end << '\n' << "r_step: " << input.r_step << endl;
  cout << "theta_start: " << setprecision(4) << input.theta_start*180/Pi()
       << '\n' << "theta_end: " << input.theta_end*180/Pi()
       << setprecision(prec) << endl;
  // Name of output file:
  cout << "Output file name: " << input.output << '\n' << endl;
 
  return;
}


// FUNCTION TO FIND EIGENSTATES WITH SPECIFIC OMEGA0 AND B VALUES:
istream& Find_States(istream& in, Input& input, Header& h, 
		     vector<EigenstateInfo>& ES)  {
  string skip;
  // Gather header information:
  in >> h.N >> h.Lz >> h.tS >> h.tSz >> h.Ez >> h.numOmega0;
  getline(in, skip);
  
  // We want to collect info only for specified omega0 and B in input file
  bool wrongOmega = true;
  bool wrongB = true;
  int numOmegaRead = 0;
  int numBRead = 0;
  while (wrongOmega) {
    // Read current omega in file:
    double fileOmega;
    in >> fileOmega >> h.numB;
    getline(in, skip); // (Skip rest of line)
    numOmegaRead += 1;
    // If it is not the correct omega,
    if (fileOmega != h.omega0) {
      // And if it is the last omega in the file,
      if (numOmegaRead == h.numOmega0) {
	// Then this file does not contain our requested omega/B combo:
	cout << "Requested omega0 not found in file." << endl;
	break;	    
      }
      // Otherwise skip to next omega in .dat file:
      // Skip data for each B field value:
      next_omega0(in, h.numB);
    }
    // If it *is* the correct omega, check B:
    if (fileOmega == h.omega0) {
      wrongOmega = false;
      while (wrongB) {
	// Read in current B in file:
	double fileB;
	int MPsize, numEigs;
	// Get B field value, MPbasis size and # eigenvalues in .dat file
	in >> fileB >> MPsize >> numEigs;
	getline(in, skip); // (Skip rest of line)
	numBRead += 1;
	// If it is not the correct B,
	if (fileB != h.B) {
	  // And if it is the last B in the file,
	  if (numBRead == h.numB) {
	    // Then this file does not contain our requested omega/B combo:
	    cout << "Requested B not found in file." << endl;
	    break;
	  }
	  // Otherwise skip to next B
	  next_B(in, MPsize, numEigs);
	}
	// If it *is* the correct B,
	if (fileB == h.B) {
	  wrongB = false;
	  cout << "Successfully found requested omega and B in file." << endl;
	  h.MPsize = MPsize;
	  h.numEigs = numEigs;
	  // Pass omega0 and B to the FD Parameter set:
// 	      // Don't know if I need this (calculates omega+, omega-, gMuB)
// 	      FDParam fdParam = genFDParamGaAs(h.omega0, h.B);
	  
	  // For each degenerate state, store state and other info as 
	  // EigenstateInfo objects in ES:
	  get_Eigenstates(in, h, ES);   
	}
      } // End loop over B
    }
  } // End loop over omega
  return in;
}

// FUNCTION TO SKIP OVER DATA WITH UNWANTED OMEGA0 VALUE:
istream& next_omega0(istream& in, int numB)
{
  string skip;
  // For each B field value, skip data:
  for (int i = 0; i != numB; ++i) {
    double fileB;
    int MPsize, numEigs;
    in >> fileB >> MPsize >> numEigs;
    getline(in, skip); // (Skip rest of line)
    // Skip basis states:
    for (int j = 0; j != MPsize; ++j) {
      getline(in, skip);
    }
    // Skip energies and coefficients:
    for (int j = 0; j != numEigs; ++j) {
      string skipEnergy, skipCoeffs;
      getline(in, skipEnergy);
      for (int k = 0; k != MPsize; ++k) {
	getline(in, skipCoeffs);
      }
    }
  }
  return in;
}


// FUNCTION TO SKIP OVER DATA WITH UNWANTED B VALUE:
istream& next_B(istream& in, int MPsize, int numEigs)
{
  string skip;
  // Skip basis states:
  for (int i = 0; i != MPsize; ++i) {
    getline(in, skip);
  }
  // Skip the coefficients for each eigenvalue:
  for (int i = 0; i != numEigs; ++i) {
    string skipEnergy, skipCoeffs;
    getline(in, skipEnergy);
    for (int j = 0; j != MPsize; ++j) {
      getline(in, skipCoeffs); // Skip the coefficients
    }
  }
  return in;
}


// FUNCTION TO READ AND STORE ANTISYMM BASIS STATES IN .dat FILE:
istream& antiSymmetrize(istream& in, Header& h, 
			LinCombState<antiSymmState, double>& basisStates)
{
  // Get rid of any previous content in the LinCombState:
  // basisStates.clear();
  if (basisStates.size() != 0)
    std::cerr << __FILE__ << ": " << __LINE__ << ": warning: "
	      << "Received non-empty LinCombState.  "
	      << "Continuing ... bravely ...\n";

  // Read each sp state: 
  for (int i = 0; i != h.MPsize; ++i) { // For each basis state in .dat file
    // Create an empty AntiSymmState:
    FDQNumLess fdLess;
    antiSymmState basis(fdLess);

    // Create an initial coefficient = +1 for now.  Later, when FDStates will
    // be added to the AntiSymmState, it will be multiplied by the appropriate
    // phase.
    double coeff = 1;

    for (int j = 0; j != h.N; ++j) { // For each sp state in each basis state
      // Create an FDstate:
      int n, m, s;
      in >> n >> m >> s;
      FDState fdState(n, m, s);
      
      // Place FDstate into the AntiSymmState:
      coeff *= basis.create(fdState);

      // The AntiSymmState should be well-defined
      if (CloseAbs<double>()(coeff, 0)) {
	std::cerr << __FILE__ << ": " << __LINE__ << ": critical error: "
		  << "Received poorly formed AntiSymmState.  Aborting.\n";
	abort();
      }

    } // End loop over sp states
     
    // Store the AntiSymmState in the antiSymmState vector:
    basisStates.add(basis, coeff);

    // Clear any end-of-line garbage:
    string junk;
    getline(in, junk);

  } // End loop over basis states

  return in;
}


// FUNCTION TO READ COEFFICIENTS OF BASIS STATES IN EACH EIGENSTATE:
istream& read_coeffs(istream& in, Header& h, vector<double>& coeffs)
{
  // Get rid of any previous content in the coeffs vector:
  coeffs.clear();
  double coefficient;
  string junk;
  for (int i = 0; i != h.MPsize; ++i) {
    in >> coefficient;
    // Skip any end-of-line junk:
    getline(in, junk);
    coeffs.push_back(coefficient);
  } // End loop over coefficients
 
  return in;
}

// Sort a vector of QMPair by increasing absolute coefficients:
bool compare(const QMPair<antiSymmState, double>& l, 
	     const QMPair<antiSymmState, double>& r)
{
  return fabs(l.coeff()) < fabs(r.coeff());
}

// FUNCTION TO REMOVE STATES FROM LINCOMBSTATE WHOSE COEFFS ARE LESS THAN CTol:
void remove_States(double& CTol, LinCombState<antiSymmState, double>& lcState, 
		   double& sumOfSquares)
{
  // Copy LinCombState to a vector of QMPairs:
  typedef QMPair<antiSymmState, double> qmPair;
  vector<qmPair> vecPairs;
  for (LinCombState<antiSymmState, double>::const_iterator iter = lcState.begin();
       iter != lcState.end(); ++iter) {
    qmPair qmpair(*iter);
    vecPairs.push_back(qmpair);
  }
  
  // Sort the vector of QMPairs in order of increasing abs(coeffs)
  sort(vecPairs.begin(), vecPairs.end(), compare);

  // Add the squares of the ordered coefficients until they
  // are less than or equal to CTol:
  vector<antiSymmState> statesToRemove;
  // Initialize sumOfSquares to zero:
  sumOfSquares = 0;
  int index = 0;
  bool LessThanCTol = true;
  while (LessThanCTol) {
    // Add the square of the coeff to sumOfSquares:
    sumOfSquares += pow(vecPairs[index].coeff(), 2);
    // Check if sumOfSquares is less than or equal to CTol:
    LessThanCTol = (sumOfSquares <= CTol);
    // If it is, then append the corresponding state to statesToRemove:
    if (LessThanCTol) {
      statesToRemove.push_back(vecPairs[index].state());
    }
    ++index;
  }

  // If the addition of the square of the last coeff makes sumOfSquares
  // *exceed* CTol, we must remove this last square from sumOfSquares
  // to keep the coeffs contributing to sumOfSquares consistent with
  // the states in statesToRemove:
  sumOfSquares -= pow(vecPairs[index-1].coeff(), 2);
  
  // Remove these states from the LinCombState:
  for (size_t i = 0; i != statesToRemove.size(); ++i) {
    lcState.remove(statesToRemove[i]);
  }
  // Sort the remaining states...
  vector<qmPair> remainingStates;
  for (LinCombState<antiSymmState, double>::const_iterator iter = lcState.begin();
       iter != lcState.end(); ++iter) {
    qmPair qmpair(*iter);
    remainingStates.push_back(qmpair);
  }
  sort(remainingStates.begin(), remainingStates.end(), compare);
  // ... and print states to standard out
  for (size_t i = 0; i != remainingStates.size(); ++i) {
    cout << remainingStates[i].state() << "  " << remainingStates[i].coeff() << endl;
  }
  cout << endl;
  return;
}

// FUNCTION TO COLLECT EIGENSTATE INFORMATION FROM .dat FILE:
istream& get_Eigenstates(istream& in, Header& h, vector<EigenstateInfo>& eStates)
{
  // Number of degenerate states in the .dat file:
  h.numDegenStates = 0;
  
  // Number of energies read so far in the .dat file:
  int numEnergiesChecked = 0;

  // Define storage for Antisymmetrized Basis States:
  // vector<antiSymmState> basisStates;
  LinCombState<antiSymmState, double> basisStates;
  // Read and store the AntiSymmStates:
  antiSymmetrize(in, h, basisStates);

  // Initialize energy tolerance for determining energy equality:
  CloseAbs<double> close(h.ETol);

  // Check energies of states to determine if they are part of groundstate set:
  bool goodEnergy = true;
  while (goodEnergy) {
    
    // Read in energy of state:
    double fileE;
    in >> fileE;
    string skip;
    getline(in, skip); // (skip energy header)
    
    // Keep track of energies read so far:
    numEnergiesChecked += 1;
    
    bool EEqual = close(fileE, h.Eg);
    bool ELessThan = (fileE < h.Eg) && !EEqual; 
    goodEnergy = (ELessThan || EEqual);

    // If fileE is greater than Eg, then this state and any remaining in the .dat
    // file are NOT in the groundstate manifold:
    if (!goodEnergy)
      break;

    // We require EEqual != ELessThan
    if (EEqual == ELessThan) {
      std::cerr << "Woah! Logic error.  Aborting.\n";
      abort();
    }
    
    // If fileE is equal to the previously stored energy, then we have 
    // found another degenerate state:
    if (EEqual) {
      // Add this degenerate state to the total counted so far:
      ++(h.numDegenStates); // h.numDegenStates += 1;
      // Create EigenstateInfo object for storing state details:
      EigenstateInfo eigenstateInfo;
      eigenstateInfo.h = h;
      eigenstateInfo.energy = fileE; // Store state energy
      cout << "State energy: " << fileE << " meV." << endl; // Write out energy
      // Read and store basis state coefficients:
      vector<double> coeffs;
      read_coeffs(in, h, coeffs);
      // Copy the full eigenstate into a LinCombState
      // ... first the states and the phases
      LinCombState<antiSymmState, double> lcState = basisStates;
      // ... now add in the coefficients
      for (size_t m = 0; m != coeffs.size(); ++m) {
	lcState[m].coeff() *= coeffs[m];
      }
      // Remove basisStates from LinCombState with coefficients whose squares
      // sum to smaller than CTol:
      remove_States(h.CTol, lcState, h.sumOfSquares);
      // Store linCombState in eigenstateInfo:
      eigenstateInfo.linCombState = lcState;
      // Update sumOfSquares in eigenstateInfo:
      eigenstateInfo.h.sumOfSquares = h.sumOfSquares;
      // Store eigenstateInfo object:
      eStates.push_back(eigenstateInfo);
    }
    
    // If fileE is less than Eg, any previously stored states are not truly in the
    // (ground) state manifold:
    if (ELessThan) {
      // Empty vector storing degenerate (ground) states and their .dat files:
      eStates.clear();
      // Replace the old (ground) state energy this state energy:
      h.Eg = fileE;
      // Keep track of the number of denerate states in this .dat file:
      ++(h.numDegenStates); // h.numDegenStates += 1;
      // Create EigenstateInfo object for storing state details:
      EigenstateInfo eigenstateInfo;
      eigenstateInfo.h = h;
      eigenstateInfo.energy = fileE; // Store state energy
      cout << "State energy: " << fileE << " meV." << endl; // Write out energy
      // Read and store basis state coefficients:
      vector<double> coeffs;
      read_coeffs(in, h, coeffs);
      // Copy the full eigenstate into a LinCombState
      // ... first the states and the phases
      LinCombState<antiSymmState, double> lcState = basisStates;
      // ... now add in the coefficients
      for (size_t m = 0; m != coeffs.size(); ++m) {
	lcState[m].coeff() *= coeffs[m];
      }
      // Remove basisStates from LinCombState with coefficients whose squares
      // sum to smaller than CTol:
      remove_States(h.CTol, lcState, h.sumOfSquares);
      // Store linCombState in eigenstateInfo:
      eigenstateInfo.linCombState = lcState;
      // Update sumOfSquares in eigenstateInfo:
      eigenstateInfo.h.sumOfSquares = h.sumOfSquares;
      // Store eigenstateInfo object:
      eStates.push_back(eigenstateInfo);
    }

    // If the last energy we read was the last one in the file, then we are now at 
    // EOF, i.e. there are no more energies to read:
    if (numEnergiesChecked == h.numEigs)
      break;
  }    
  return in;
}
