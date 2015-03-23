/*
 *  file   twoPointSpin.cpp
 *  
 *  Calculates two-point spin correlations in real space 
 * for the global ground state from a series of .dat files.
 *
 *  !!! Assumes .dat file is structured in order of INCREASING 
 *  variable values !!!
 *
 *  Catherine J. Stevenson 
 *  Jan 23, 2008 - Modified to enable threading
 *  February 6 - Modified to output both <S+zS+z> and <S+zS-z> ME concurrently.
 *             - # iterations through theta depends on arc length now.
 *
 */

#include "twoPoint.hpp"
#include "DiagonHeaders.h"
#include "ThreadManager.hpp"
#include <boost/function.hpp>
#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <ios>
#include <limits>

using std::cout;        using std::cin;           using std::endl;
using std::vector;      using std::string;        using std::ofstream;
using std::streamsize;  using std::setprecision;  using std::ifstream;
using std::complex;

int main(int argc, char** argv) {
  // Program requires an input file:
  if (argc != 2) {
    cout << "prog: usage: prog infile\n";
    abort();
  }

  // Collect input information:
  ifstream paramFile(argv[1]);
  Input input;
  Read_Input(paramFile, input);  
  // Write input data to stdout:
  Write_Stdout(input);

//////////////////////////////////////////
//  GET EIGENSTATES FROM EACH .dat FILE //
//////////////////////////////////////////
  // Data storage for degenerate Eigenstates:
  vector<EigenstateInfo> eigenStates;
  // Header info storage:
  Header h;

  // Write a dummy energy (to be replaced with true (ground) state
  // energy after reading data), other input data to header:
  h.Eg = std::numeric_limits<double>::max();
  h.ETol = input.ETol;
  h.CTol = input.CTol;
  h.omega0 = input.input_omega0;
  h.B = input.input_B;

  // Check each .dat file for degenerate (ground) states:
  for (size_t i = 0; i != input.datFiles.size(); ++i) {    
    cout << "Now reading file " << input.datFiles[i] << endl;
    // Write .dat file name to header:
    h.file = input.datFiles[i];
    
    // Open .dat file:
    ifstream in(input.datFiles[i].c_str());

    // Search for states with the omega0/B values specified in input:
    if (in) {
      Find_States(in, input, h, eigenStates);
    } else {
      cout << "Unable to read .dat file " << input.datFiles[i] << endl;
    }

    // Close .dat file:
    in.close();
//     cout << '\n' << endl;

  } // End loop over datFiles

  // Write data to stdout:
  cout << "States used in calculation:" << endl;
  for (size_t i = 0; i != eigenStates.size(); ++i) {
    double newNorm = 1 - eigenStates[i].h.sumOfSquares;
    cout << "File: " << eigenStates[i].h.file << '\t' << "Energy: " 
	 << eigenStates[i].energy << " meV" << endl;
    cout << "Normalization: " << newNorm << '\t' << "Size: "
	 << eigenStates[i].linCombState.size() << '\n' << endl;
  }
  cout << endl;
  
/////////////////////////////////////
// GET TWO-POINT SPIN CORRELATIONS //
/////////////////////////////////////
  // Open an output file:
  ofstream out(input.output.c_str());

  // Machine's default precision (for resetting precision later)
  streamsize prec = cout.precision();

 // Write to file the eigenstate info:
  out << "# " << h.N << "  " << h.omega0 << "  " << h.B << "  " << h.Eg 
      << "  ! N omega0 B Energy" << endl;

  // Define probe coordinates: Must be POLAR COORDINATES: r must have same
  // units as lo below (currently nanometers), theta is in radian (theta is
  // input in degrees, then stored in the Input structure as radians).
  double r, theta;
  // Write to output file the caculation done:
  out << "#  r0 = " << input.r0 << ", theta0 = "
      << setprecision(4) << input.theta0*180/Pi() << setprecision(prec) 
      << endl; 
  
  // Define effective length scale (lo) in GaAs (in units of nm):
  double lo = sqrt(pow(GaAs::a0,2)*GaAs::Ry/sqrt(pow(h.omega0,2)
	      + 0.25*pow(h.B*GaAs::meV_per_Tesla_orbital,2)));

  // Instantiate all twoPointCalc objects for each (r,theta) coordinate:
  vector<twoPointCalc> correlations;
  for (r = input.r_start; r <= input.r_end; r += input.r_step) {
    if (r == 0) {
      theta = 0;
      // Initiate the matrix elements:
      SzSzME<1> SzSzPlus(input.r0, input.theta0, r, theta, lo);
      SzSzME<-1> SzSzMinus(input.r0, input.theta0, r, theta, lo);      
      // Instantiate a twoPointCalc object for this (r,theta) coordinate,
      twoPointCalc corrCalc(&eigenStates, SzSzPlus, SzSzMinus);
      // and store for later calculation:
      correlations.push_back(corrCalc);
    } else { // for all non-zero r, increment theta via a constant arc length:
      double arc = (input.arcLength)*Pi();
      // Go 3 steps past theta_end for gnuplot:
      for (theta = input.theta_start; theta <= input.theta_end + 2*arc/r;
	   theta += arc/r) {
      // Initiate the matrix elements:
      SzSzME<1> SzSzPlus(input.r0, input.theta0, r, theta, lo);
      SzSzME<-1> SzSzMinus(input.r0, input.theta0, r, theta, lo);
      
      // Instantiate a twoPointCalc object for this (r,theta) coordinate,
      twoPointCalc corrCalc(&eigenStates, SzSzPlus, SzSzMinus);
      // and store for later calculation:
      correlations.push_back(corrCalc);
      } // End loop over theta values for larger r
    }
  } // end r loop

  ThreadManager threadMan(input.numThreads);
  
  // Calculate the correlations at every (r,theta) coordinate:
  for (vector<twoPointCalc>::iterator i = correlations.begin(); i != correlations.end();
       ++ i) {
    threadMan.add_job(boost::ref(*i));
  }
  
  threadMan.wait_all(); // Wait until all jobs (calculations) are done.

  // Write coordinates, correlation data to out file:
  out << "# r \t theta \t <S+zS+z>.real \t \t \t <S+zS+z>.imag \t \t \t <S+zS-z>.real"
      << "\t \t \t" << "<S+zS-z>.imag" << endl;
  for (size_t i = 0; i != correlations.size(); ++i) {
    out << correlations[i].r1() << '\t' << setprecision(5)
	<< (correlations[i].theta1())*180/Pi() << '\t'
	<< setprecision(12)
	<< correlations[i].plusCorrelation_real() << '\t' << '\t' 
	<< correlations[i].plusCorrelation_imag() << '\t' << '\t' 
	<< correlations[i].minusCorrelation_real() << '\t' << '\t'
	<< correlations[i].minusCorrelation_imag() 
	<< setprecision(prec) << endl;
    // If we're not at the last set of values, and the next set of values are for a new
    // r, print a blank line for gnuplot first:
    if (i+1 != correlations.size() && correlations[i].r1() != correlations[i+1].r1()) {
      // Block line for gnuplot:
      out << endl;
    }
  }

  cout << "Done." << endl;
  
  // Close file:
  out.close();
  
  return 0;
}
