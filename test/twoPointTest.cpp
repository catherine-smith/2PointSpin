/*
 *  twoPointTest.cpp
 *
 *  Note the bulk of the testing was done for 1PointSpin.cpp.
 *  Here we only test the twoPoint matrix element and the Find_States
 *  function.
 *
 *  Catherine J. Stevenson
 *  Mar 13 2008 - Modified test for SYMMETRIC 2-Pt ME
 *              - Symmetric TwoBodySOp is used.
 *
 */

#include "twoPoint.hpp"
#include "DiagonHeaders.h"
#include "CloseAbs.hpp"
#include <boost/test/unit_test.hpp>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>

using boost::unit_test::test_suite;
using std::complex;      using std::fabs;    using std::cin;
using std::endl;         using std::string;  using std::vector;
using std::ifstream;     using std::cout;    using std::sqrt;
using std::setprecision; using std::streamsize;

void twoPointTest() {

  // Instantiate CloseAbs:
  CloseAbs<double> closer(1.e-23);
  CloseAbs<double> close(1.e-10);

  // Test parameters:
  double w_1 = 1.0; //hbar*omega0 = 1 meV
  double B_0 = 0; // B = 0 T
//   double r_0 = 0;
//   double r_1 = 1;
//   double r_3 = 3;
//   double r_11 = 11;
//   double r_15 = 15;
//   double r_30 = 30;
//   double t_0 = 0;
  double pi_5 = Pi()/5;
  double pi_4 = Pi()/4;
  double pi_2 = Pi()/2;
  double twoPi_3 = 2*Pi()/3;
  // Characteristic length for B = 0 T and h-bar*w = 1.0 meV:
  double lo_0 = sqrt(pow(GaAs::a0,2)*GaAs::Ry/sqrt(pow(w_1,2)
		   +0.25*pow(B_0*GaAs::meV_per_Tesla_orbital,2)));
  BOOST_CHECK(closer(lo_0, 23.84646922721940509121778));

  // Some single particle states:
  FDState fd00u(0, 0, 1), fd00d(0, 0, -1), fd01u(0, 1, 1), fd01d(0, 1, -1),
    fd02u(0, 2, 1), fd02d(0, 2, -1), fd10u(1, 0, 1), fd10d(1, 0, -1), 
    fd11u(1, 1, 1), fd11d(1, 1, -1), fd12u(1, 2, 1), fd12d(1, 2, -1),
    fd13u(1, 3, 1), fd13d(1, 3, -1), fd14u(1, 4, 1), fd14d(1, 4, -1),
    fd20u(2, 0, 1), fd20d(2, 0, -1), fd21u(2, 1, 1), fd21d(2, 1, -1),
    fd26u(2, 6, 1), fd26d(2, 6, -1), fd31u(3, 1, 1), fd31d(3, 1, -1);


  // ##### CHECK 1 ##### Does S+z(r0)S+/-z(r1) ME behave properly?
  // A. S+z(r0)S+z(r1) ME (B = 0 T, omega0 = 1meV)   
  // Instantiate 4 MEs at 4 different locations:
  // r_0 = (0, 0), r_1 = (1, 0):  
  SzSzME<1> Plus01(0, 0, 1, 0, lo_0);
  BOOST_CHECK_EQUAL(Plus01.r0(), 0);
  BOOST_CHECK_EQUAL(Plus01.theta0(), 0);
  BOOST_CHECK_EQUAL(Plus01.r1(), 1);
  BOOST_CHECK_EQUAL(Plus01.theta1(), 0);
  // r_0 = (1, 0), r_1 = (3, pi/4):
  SzSzME<1> Plus02(1, 0, 3, pi_4, lo_0);
  BOOST_CHECK_EQUAL(Plus02.r0(), 1);
  BOOST_CHECK_EQUAL(Plus02.theta0(), 0);
  BOOST_CHECK_EQUAL(Plus02.r1(), 3);
  BOOST_CHECK_EQUAL(Plus02.theta1(), pi_4);
  // r_0 = (11, pi/4), r_1 = (15, 2pi/3):  
  SzSzME<1> Plus03(11, pi_4, 15, twoPi_3, lo_0);
  BOOST_CHECK_EQUAL(Plus03.r0(), 11);
  BOOST_CHECK_EQUAL(Plus03.theta0(), pi_4);
  BOOST_CHECK_EQUAL(Plus03.r1(), 15);
  BOOST_CHECK_EQUAL(Plus03.theta1(), twoPi_3);
  // r_0 = (30, pi_5), r_1 = (30, pi_2):  
  SzSzME<1> Plus04(30, pi_5, 30, pi_2, lo_0);
  BOOST_CHECK_EQUAL(Plus04.r0(), 30);
  BOOST_CHECK_EQUAL(Plus04.theta0(), pi_5);
  BOOST_CHECK_EQUAL(Plus04.r1(), 30);
  BOOST_CHECK_EQUAL(Plus04.theta1(), pi_2);

  // No. 1: (11^, 21^|S+z(r_0)S+z(r_1)|11^, 21v) = 0
  // Spin not conserved.
  complex<double> MEPlus03_0 = Plus03(fd11u, fd21u, fd11u, fd21d);
  BOOST_CHECK_EQUAL(MEPlus03_0.real(), 0);
  BOOST_CHECK_EQUAL(MEPlus03_0.imag(), 0);

  // No. 2: (26^, 31v|S+z(r_0)S+z(r_1)|31v, 26^) = 0  
  // Spin not conserved.
  complex<double> MEPlus04_0 = Plus04(fd26u, fd31d, fd31d, fd26u);
  BOOST_CHECK_EQUAL(MEPlus04_0.real(), 0);
  BOOST_CHECK_EQUAL(MEPlus04_0.imag(), 0);

  // No. 3: (20^, 14v|S+z(r_0)S+z(r_1)|20^, 14v) = 0  
  // All spins must be up for S+z(r_0)S+z(r_1)
  complex<double> MEPlus02_0 = Plus02(fd20u, fd14d, fd20u, fd14d);
  BOOST_CHECK_EQUAL(MEPlus02_0.real(), 0);
  BOOST_CHECK_EQUAL(MEPlus02_0.imag(), 0);

  // No. 4: (00^, 01^|S+z(r_0)S+z(r_1)|00^, 01^)
  complex<double> MEPlus01_1 = Plus01(fd00u, fd01u, fd00u, fd01u);
  BOOST_CHECK(closer(MEPlus01_1.real(), 1.72037808884308130207676694417679e-11));
  BOOST_CHECK(closer(MEPlus01_1.imag(), 0));
  complex<double> MEPlus02_1 = Plus02(fd00u, fd01u, fd00u, fd01u);
  BOOST_CHECK(closer(MEPlus02_1.real(), 1.7068177356555718971228535065414e-10));
  BOOST_CHECK(closer(MEPlus02_1.imag(), 0));
  complex<double> MEPlus03_1 = Plus03(fd00u, fd01u, fd00u, fd01u);
  BOOST_CHECK(closer(MEPlus03_1.real(), 4.394988082562185397561530832692e-9));
  BOOST_CHECK(closer(MEPlus03_1.imag(), 0));
  complex<double> MEPlus04_1 = Plus04(fd00u, fd01u, fd00u, fd01u);
  BOOST_CHECK(closer(MEPlus04_1.real(), 6.36688698286884420579369125654e-9));
  BOOST_CHECK(closer(MEPlus04_1.imag(), 0));

  // No. 5: (00^, 01^|S+z(r_0)S+z(r_1)|01^, 00^)
  complex<double> MEPlus01_2 = Plus01(fd00u, fd01u, fd01u, fd00u);
  BOOST_CHECK_EQUAL(MEPlus01_2.real(), 0);
  BOOST_CHECK_EQUAL(MEPlus01_2.imag(), 0);
  complex<double> MEPlus02_2 = Plus02(fd00u, fd01u, fd01u, fd00u);
  BOOST_CHECK(closer(MEPlus02_2.real(), 7.24141437078913795247163495140e-11));
  BOOST_CHECK(closer(MEPlus02_2.imag(), 0));
  complex<double> MEPlus03_2 = Plus03(fd00u, fd01u, fd01u, fd00u);
  BOOST_CHECK(closer(MEPlus03_2.real(), 1.08490515662624900904040165728e-9));
  BOOST_CHECK(closer(MEPlus03_2.imag(), 0));
  complex<double> MEPlus04_2 = Plus04(fd00u, fd01u, fd01u, fd00u);
  BOOST_CHECK(closer(MEPlus04_2.real(), 3.74236227154322663339433668327e-9));
  BOOST_CHECK(closer(MEPlus04_2.imag(), 0));

  // No. 6: (01^, 11^|S+z(r_0)S+z(r_1)|11^, 13^)
  complex<double> MEPlus01_3 = Plus01(fd01u, fd11u, fd11u, fd13u);
  BOOST_CHECK_EQUAL(MEPlus01_3.real(), 0);
  BOOST_CHECK_EQUAL(MEPlus01_3.imag(), 0);
  complex<double> MEPlus02_3 = Plus02(fd01u, fd11u, fd11u, fd13u);
  BOOST_CHECK(close(MEPlus02_3.real(), -1.30298731752361938343696897886e-12));
  BOOST_CHECK(close(MEPlus02_3.imag(), 4.21215404003235331371864627774e-12));
  complex<double> MEPlus03_3 = Plus03(fd01u, fd11u, fd11u, fd13u);
  BOOST_CHECK(close(MEPlus03_3.real(), 2.37066375064275534708238614435e-10));
  BOOST_CHECK(close(MEPlus03_3.imag(), -4.87296590301574820271430885702e-10));
  complex<double> MEPlus04_3 = Plus04(fd01u, fd11u, fd11u, fd13u);
  BOOST_CHECK(close(MEPlus04_3.real(), -1.57910292847137520632904382995e-11));
  BOOST_CHECK(close(MEPlus04_3.imag(), 3.09916399669186933557420205225e-11));

  // B. S+z(r0)S-z(r1) ME (B = 0 T, omega0 = 1meV)   
  // Instantiate 4 MEs at 4 different locations:
  // r_0 = (0, 0), r_1 = (1, 0):  
  SzSzME<-1> Minus01(0, 0, 1, 0, lo_0);
  BOOST_CHECK_EQUAL(Minus01.r0(), 0);
  BOOST_CHECK_EQUAL(Minus01.theta0(), 0);
  BOOST_CHECK_EQUAL(Minus01.r1(), 1);
  BOOST_CHECK_EQUAL(Minus01.theta1(), 0);
  // r_0 = (1, 0), r_1 = (3, pi/4):
  SzSzME<-1> Minus02(1, 0, 3, pi_4, lo_0);
  BOOST_CHECK_EQUAL(Minus02.r0(), 1);
  BOOST_CHECK_EQUAL(Minus02.theta0(), 0);
  BOOST_CHECK_EQUAL(Minus02.r1(), 3);
  BOOST_CHECK_EQUAL(Minus02.theta1(), pi_4);
  // r_0 = (11, pi/4), r_1 = (15, 2pi/3):  
  SzSzME<-1> Minus03(11, pi_4, 15, twoPi_3, lo_0);
  BOOST_CHECK_EQUAL(Minus03.r0(), 11);
  BOOST_CHECK_EQUAL(Minus03.theta0(), pi_4);
  BOOST_CHECK_EQUAL(Minus03.r1(), 15);
  BOOST_CHECK_EQUAL(Minus03.theta1(), twoPi_3);
  // r_0 = (30, pi_5), r_1 = (30, pi_2):  
  SzSzME<-1> Minus04(30, pi_5, 30, pi_2, lo_0);
  BOOST_CHECK_EQUAL(Minus04.r0(), 30);
  BOOST_CHECK_EQUAL(Minus04.theta0(), pi_5);
  BOOST_CHECK_EQUAL(Minus04.r1(), 30);
  BOOST_CHECK_EQUAL(Minus04.theta1(), pi_2);

  // No. 7: (11^, 21^|S+z(r_0)S+z(r_1)|11^, 21v) = 0
  // Spin not conserved.
  complex<double> MEMinus03_0 = Minus03(fd11u, fd21u, fd11u, fd21d);
  BOOST_CHECK_EQUAL(MEMinus03_0.real(), 0);
  BOOST_CHECK_EQUAL(MEMinus03_0.imag(), 0);

  // No. 8: (26^, 31v|S+z(r_0)S+z(r_1)|31v, 26^) = 0  
  // Spin not conserved.
  complex<double> MEMinus04_0 = Minus04(fd26u, fd31d, fd31d, fd26u);
  BOOST_CHECK_EQUAL(MEMinus04_0.real(), 0);
  BOOST_CHECK_EQUAL(MEMinus04_0.imag(), 0);

  // No. 9: (20^, 14v|S+z(r_0)S+z(r_1)|20^, 14v) = 0  
  // All spins must be up for S+z(r_0)S+z(r_1)
  complex<double> MEMinus02_0 = Minus02(fd20u, fd14u, fd20u, fd14u);
  BOOST_CHECK_EQUAL(MEMinus02_0.real(), 0);
  BOOST_CHECK_EQUAL(MEMinus02_0.imag(), 0);

  // No. 10: (01^, 26v|S+z(r_0)S+z(r_1)|02^, 14v)
  complex<double> MEMinus01_1 = Minus01(fd01u, fd26d, fd02u, fd14d);
  BOOST_CHECK_EQUAL(MEMinus01_1.real(), 0);
  BOOST_CHECK_EQUAL(MEMinus01_1.imag(), 0);
  complex<double> MEMinus02_1 = Minus02(fd01u, fd26d, fd02u, fd14d);
  BOOST_CHECK(closer(MEMinus02_1.real(), -7.16393497492944441671801461062e-21));
  BOOST_CHECK(closer(MEMinus02_1.imag(), -7.16393497492944441671801461062e-21));
  complex<double> MEMinus03_1 = Minus03(fd01u, fd26d, fd02u, fd14d);
  BOOST_CHECK(closer(MEMinus03_1.real(), -1.78790115364491660629880654510e-13));
  BOOST_CHECK(closer(MEMinus03_1.imag(), -6.67253794431383838136528141442e-13));
  complex<double> MEMinus04_1 = Minus04(fd01u, fd26d, fd02u, fd14d);
  BOOST_CHECK(closer(MEMinus04_1.real(), -1.89232885956177748801555574112e-10));
  BOOST_CHECK(closer(MEMinus04_1.imag(), -2.60456722988666494177896410128e-10));

  // No. 11: (01v, 26^|S+z(r_0)S+z(r_1)|02v, 14^)
  complex<double> MEMinus01_2 = Minus01(fd01d, fd26u, fd02d, fd14u);
  BOOST_CHECK_EQUAL(MEMinus01_2.real(), 0);
  BOOST_CHECK_EQUAL(MEMinus01_2.imag(), 0);
  complex<double> MEMinus02_2 = Minus02(fd01d, fd26u, fd02d, fd14u);
  BOOST_CHECK(closer(MEMinus02_2.real(), -8.88493782660097640913436630035e-23));
  BOOST_CHECK(closer(MEMinus02_2.imag(), 8.88493782660097640913436630035e-23));
  complex<double> MEMinus03_2 = Minus03(fd01d, fd26u, fd02d, fd14u);
  BOOST_CHECK(closer(MEMinus03_2.real(), -5.49975843546383833053828493530e-14));
  BOOST_CHECK(closer(MEMinus03_2.imag(), 2.05253779105065629546243500521e-13));
  complex<double> MEMinus04_2 = Minus04(fd01d, fd26u, fd02d, fd14u);
  BOOST_CHECK(closer(MEMinus04_2.real(), -1.89232885956177748801555574113e-10));
  BOOST_CHECK(closer(MEMinus04_2.imag(), 2.60456722988666494177896410128e-10));


  // ##### CHECK 2 ##### Check MEs using twoBodySOp:
  // A. <S+z(r0)S+z(r1)>  (B = 0 T, omega0 = 1meV)
  // No. 12: <00^, 01^|S+z(r_0)S+z(r_1)|00^, 01^>:
  // Antisymm state:
  FDQNumLess fdLess;
  antiSymmState aSymm01(fdLess);
  BOOST_CHECK_EQUAL(aSymm01.create(fd00u), 1);
  BOOST_CHECK_EQUAL(aSymm01.create(fd01u), -1);
  // twoBodySOp results:
  complex<double> corrPlus01_1 = twoBodySOp(aSymm01, aSymm01, Plus01);
  BOOST_CHECK(closer(corrPlus01_1.real(), 1.72037808884308130207676694418e-11));
  BOOST_CHECK(closer(corrPlus01_1.imag(), 0));
  complex<double> corrPlus02_1 = twoBodySOp(aSymm01, aSymm01, Plus02);
  BOOST_CHECK(closer(corrPlus02_1.real(), 9.82676298576658101875690011410e-11));
  BOOST_CHECK(closer(corrPlus02_1.imag(), 0));
  complex<double> corrPlus03_1 = twoBodySOp(aSymm01, aSymm01, Plus03);
  BOOST_CHECK(closer(corrPlus03_1.real(), 3.31008292593593638852112917544e-9));
  BOOST_CHECK(closer(corrPlus03_1.imag(), 0));
  complex<double> corrPlus04_1 = twoBodySOp(aSymm01, aSymm01, Plus04);
  BOOST_CHECK(closer(corrPlus04_1.real(), 2.62452471132561757239935457331e-9));
  BOOST_CHECK(closer(corrPlus04_1.imag(), 0));

  // No. 13: <01^, 11^, 14v|S+z(r_0)S+z(r_1)|11^, 13^, 21v>:
  // Antisymm state:
  antiSymmState aSymm02(fdLess);
  BOOST_CHECK_EQUAL(aSymm02.create(fd01u), 1);
  BOOST_CHECK_EQUAL(aSymm02.create(fd11u), -1);
  BOOST_CHECK_EQUAL(aSymm02.create(fd14d), 1);
  antiSymmState aSymm03(fdLess);
  BOOST_CHECK_EQUAL(aSymm03.create(fd11u), 1);
  BOOST_CHECK_EQUAL(aSymm03.create(fd13u), -1);
  BOOST_CHECK_EQUAL(aSymm03.create(fd21d), 1);
  // twoBodySOp results:
  complex<double> corrPlus01_2 = twoBodySOp(aSymm02, aSymm03, Plus01);
  BOOST_CHECK_EQUAL(corrPlus01_2.real(), 0);
  BOOST_CHECK_EQUAL(corrPlus01_2.imag(), 0);
  complex<double> corrPlus03_2 = twoBodySOp(aSymm02, aSymm03, Plus03);
  BOOST_CHECK_EQUAL(corrPlus03_2.real(), 0);
  BOOST_CHECK_EQUAL(corrPlus03_2.imag(), 0);

  // No. 14: <01^, 11^, 14^|S+z(r_0)S+z(r_1)|11^, 13^, 21^>:
  // Antisymm state:
  antiSymmState aSymm04(fdLess);
  BOOST_CHECK_EQUAL(aSymm04.create(fd01u), 1);
  BOOST_CHECK_EQUAL(aSymm04.create(fd11u), -1);
  BOOST_CHECK_EQUAL(aSymm04.create(fd14u), 1);
  antiSymmState aSymm05(fdLess);
  BOOST_CHECK_EQUAL(aSymm05.create(fd11u), 1);
  BOOST_CHECK_EQUAL(aSymm05.create(fd13u), -1);
  BOOST_CHECK_EQUAL(aSymm05.create(fd21u), 1);
  // twoBodySOp results:
  complex<double> corrPlus01_3 = twoBodySOp(aSymm04, aSymm05, Plus01);
  BOOST_CHECK(close(corrPlus01_3.real(), 0));
  BOOST_CHECK(close(corrPlus01_3.imag(), 0));
  complex<double> corrPlus02_3 = twoBodySOp(aSymm04, aSymm05, Plus02);
  BOOST_CHECK(close(corrPlus02_3.real(), -1.28615386507685174735821995210e-16));
  BOOST_CHECK(close(corrPlus02_3.imag(), -1.10089730593635271714491851356e-16));
  complex<double> corrPlus03_3 = twoBodySOp(aSymm04, aSymm05, Plus03);
  BOOST_CHECK(close(corrPlus03_3.real(), 3.89716350642589266807065173862e-11));
  BOOST_CHECK(close(corrPlus03_3.imag(), 6.54127445209308363485852935080e-11));
  complex<double> corrPlus04_3 = twoBodySOp(aSymm04, aSymm05, Plus04);
  BOOST_CHECK(close(corrPlus04_3.real(), 2.82568940096126463512501708336e-9));
  BOOST_CHECK(close(corrPlus04_3.imag(), 4.47545234302224126653862685693e-10));

  // B. <S+z(r0)S-z(r1)>  (B = 0 T, omega0 = 1meV)
  // No. 15 <11v, 13^|S+z(r_0)S+z(r_1)|00^, 01v>:
  // Antisymm state:
  antiSymmState aSymm06(fdLess);
  BOOST_CHECK_EQUAL(aSymm06.create(fd11d), 1);
  BOOST_CHECK_EQUAL(aSymm06.create(fd13u), -1);
  antiSymmState aSymm07(fdLess);
  BOOST_CHECK_EQUAL(aSymm07.create(fd00u), 1);
  BOOST_CHECK_EQUAL(aSymm07.create(fd01d), -1);
  // twoBodySOp results:
  complex<double> corrMinus01_1 = twoBodySOp(aSymm06, aSymm07, Minus01);
  BOOST_CHECK(closer(corrMinus01_1.real(), 0));
  BOOST_CHECK(closer(corrMinus01_1.imag(), 0));
  complex<double> corrMinus02_1 = twoBodySOp(aSymm06, aSymm07, Minus02);
  BOOST_CHECK(closer(corrMinus02_1.real(), -1.30413400231379953500955138756e-12));
  BOOST_CHECK(closer(corrMinus02_1.imag(), 1.30413400231379953500955138756e-12));
  complex<double> corrMinus03_1 = twoBodySOp(aSymm06, aSymm07, Minus03);
  BOOST_CHECK(closer(corrMinus03_1.real(), -5.61018682786428318907203165495e-10));
  BOOST_CHECK(closer(corrMinus03_1.imag(), 3.23904287527153656097583630618e-10));
  complex<double> corrMinus04_1 = twoBodySOp(aSymm06, aSymm07, Minus04);
  BOOST_CHECK(closer(corrMinus04_1.real(), -5.06725784221743660881042542488e-10));
  BOOST_CHECK(closer(corrMinus04_1.imag(), 1.64645187882710187654934584030e-10));

  // No. 16 <11^, 14v, 26^|S+z(r_0)S+z(r_1)|02v, 11^, 31^>:
  // Antisymm state:
  antiSymmState aSymm08(fdLess);
  BOOST_CHECK_EQUAL(aSymm08.create(fd11u), 1);
  BOOST_CHECK_EQUAL(aSymm08.create(fd14d), -1);
  BOOST_CHECK_EQUAL(aSymm08.create(fd26u), 1);  
  antiSymmState aSymm09(fdLess);
  BOOST_CHECK_EQUAL(aSymm09.create(fd02d), 1);
  BOOST_CHECK_EQUAL(aSymm09.create(fd11u), -1);
  BOOST_CHECK_EQUAL(aSymm09.create(fd31u), 1);
  // twoBodySOp results:
  complex<double> corrMinus01_2 = twoBodySOp(aSymm08, aSymm09, Minus01);
  BOOST_CHECK(closer(corrMinus01_2.real(), 0));
  BOOST_CHECK(closer(corrMinus01_2.imag(), 0));
  complex<double> corrMinus02_2 = twoBodySOp(aSymm08, aSymm09, Minus02);
  BOOST_CHECK(closer(corrMinus02_2.real(), -2.89808806954259157494186659360e-23));
  BOOST_CHECK(closer(corrMinus02_2.imag(), -2.89808806954259157494186659360e-23));
  complex<double> corrMinus03_2 = twoBodySOp(aSymm08, aSymm09, Minus03);
  BOOST_CHECK(closer(corrMinus03_2.real(), -1.28757200017254525821341778240e-13));
  BOOST_CHECK(closer(corrMinus03_2.imag(), -7.43380040900643876137483956182e-14));
  complex<double> corrMinus04_2 = twoBodySOp(aSymm08, aSymm09, Minus04);
  BOOST_CHECK(closer(corrMinus04_2.real(), -1.51785954620069089121779984236e-10));
  BOOST_CHECK(closer(corrMinus04_2.imag(), 2.08915443720521088706583551244e-10));


  // ##### CHECK 3 ##### Does the twoPointCalc class work properly?
  // twoPointCalc calculates the "total" correlation at a specific point
  // by averaging the correlations of all degenerate states.
  // Assume Psi_1 and Psi_2 are degenerate.
  // Make Psi_1 LinCombState:
  antiSymmState aSymm1(fdLess);
  BOOST_CHECK_EQUAL(aSymm1.create(fd00d), 1);
  BOOST_CHECK_EQUAL(aSymm1.create(fd00u), -1);
  BOOST_CHECK_EQUAL(aSymm1.create(fd01u), 1);
  antiSymmState aSymm2(fdLess);
  BOOST_CHECK_EQUAL(aSymm2.create(fd01d), 1);
  BOOST_CHECK_EQUAL(aSymm2.create(fd01u), -1);
  BOOST_CHECK_EQUAL(aSymm2.create(fd10u), 1);
  antiSymmState aSymm3(fdLess);
  BOOST_CHECK_EQUAL(aSymm3.create(fd00u), 1);
  BOOST_CHECK_EQUAL(aSymm3.create(fd11d), -1);
  BOOST_CHECK_EQUAL(aSymm3.create(fd12u), 1);
  LinCombState<antiSymmState, double> Psi_1;
  Psi_1.add(aSymm1, 1/(sqrt(2)));
  Psi_1.add(aSymm2, 0.5);
  Psi_1.add(aSymm3, 0.5);
  // Make corresponding EigenstateInfo object:
  Header head;
  head.N = 3;
  head.Lz = 1;
  head.tS = 1;
  head.tSz = 1;
  head.MPsize = 3;
  head.numDegenStates = 2;
  head.numEigs = 1;
  head.numOmega0 = 1;
  head.numB = 1;
  head.Ez = 0.012734447;
  head.omega0 = 1;
  head.B = 0;
  head.Eg = 10.8407;
  head.ETol = 1e-5;
  head.CTol = 1e-20;
  head.sumOfSquares = 0;
  head.file = "dummyfile";
  EigenstateInfo eigenstateInfo_1;
  eigenstateInfo_1.h = head;
  eigenstateInfo_1.energy = head.Eg;
  eigenstateInfo_1.linCombState = Psi_1;
  // Store eigenstateInfo_1 in a vector to be used in twoPointCalc:
  vector<EigenstateInfo> eigenStates;
  eigenStates.push_back(eigenstateInfo_1);
  //Repeat for Psi_2:
  // Make Psi_2 LinCombState:
  antiSymmState aSymm4(fdLess);
  BOOST_CHECK_EQUAL(aSymm4.create(fd01d), 1);
  BOOST_CHECK_EQUAL(aSymm4.create(fd10u), -1);
  BOOST_CHECK_EQUAL(aSymm4.create(fd12u), 1);
  LinCombState<antiSymmState, double> Psi_2;
  Psi_2.add(aSymm1, sqrt(3)/(sqrt(8)));
  Psi_2.add(aSymm4, sqrt(3)/(sqrt(8)));
  Psi_2.add(aSymm3, 0.5);
  // Make corresponding EigenstateInfo object:
  EigenstateInfo eigenstateInfo_2;
  eigenstateInfo_2.h = head;
  eigenstateInfo_2.energy = head.Eg;
  eigenstateInfo_2.linCombState = Psi_2;
  // Store eigenstateInfo_2 in a vector to be used in twoPointCalc:
  eigenStates.push_back(eigenstateInfo_2);

  // Instantiate two Matrix Elements for each possible correlation
  // calculation (<S+z(r0)S+z(r1)> or <S+z(r0)S-z(r1)>).
  SzSzME<1> MEPlus(11, pi_4, 15, twoPi_3, lo_0);
  SzSzME<-1> MEMinus(11, pi_4, 15, twoPi_3, lo_0);
  // Instantiate the twoPointCalc object:
  twoPointCalc cc(&eigenStates, MEPlus, MEMinus);
  // Do the correlation calculation:
  cc();
  // Check the results for the "averaged" correlation:
  BOOST_CHECK(closer(cc.plusCorrelation_imag(), 0));
  BOOST_CHECK(closer(cc.plusCorrelation_real(), 3.33969613895424354460235915032e-9));
  BOOST_CHECK(closer(cc.minusCorrelation_imag(), 0));
  BOOST_CHECK(closer(cc.minusCorrelation_real(), 8.84753157007167253164345953098e-9));
  BOOST_CHECK(closer(cc.r0(), 11));
  BOOST_CHECK(closer(cc.theta0(), pi_4));
  BOOST_CHECK(closer(cc.r1(), 15));
  BOOST_CHECK(closer(cc.theta1(), twoPi_3));

 // ##### CHECK 4 ##### Does the Find_States function work properly?
  // Read input:
  Input input;
  Read_Input(cin, input);
  //  Check files in list:
  size_t two = 2;
  size_t four = 4;
  size_t five = 5;
  BOOST_CHECK_EQUAL(input.numFiles, 5);  
  BOOST_CHECK_EQUAL(input.datFiles.size(), five);
  BOOST_CHECK_EQUAL(input.datFiles[0], "test1.dat");
  BOOST_CHECK_EQUAL(input.datFiles[1], "test2.dat");
  BOOST_CHECK_EQUAL(input.datFiles[2], "test3.dat");
  BOOST_CHECK_EQUAL(input.datFiles[3], "test4.dat");
  BOOST_CHECK_EQUAL(input.datFiles[4], "test5.dat");
  BOOST_CHECK_EQUAL(input.numThreads, 2);
  BOOST_CHECK_EQUAL(input.input_omega0, 0.5);
  BOOST_CHECK_EQUAL(input.input_B, 0.1);
  BOOST_CHECK_EQUAL(input.r0, 0);
  BOOST_CHECK_EQUAL(input.theta0, 0);
  BOOST_CHECK_EQUAL(input.ETol, 1e-5);
  BOOST_CHECK_EQUAL(input.CTol, 1e-3);
  BOOST_CHECK_EQUAL(input.r_start, 0);
  BOOST_CHECK_EQUAL(input.r_end, 6);
  BOOST_CHECK_EQUAL(input.r_step, 2);
  // Note input theta in degrees is converted into radians:
  BOOST_CHECK_EQUAL(input.theta_start, 0);
  BOOST_CHECK(close(input.theta_end, 0.174532925199433)); // 10 degrees
  BOOST_CHECK_EQUAL(input.arcLength, 2);
  BOOST_CHECK_EQUAL(input.output, "ReadTest.out");

  // Data storage for degenerate Eigenstates:
  vector<EigenstateInfo> eStates;
  // Header info storage:
  Header h;
  size_t zero = 0;
  BOOST_CHECK_EQUAL(eStates.size(), zero);

  // Write a dummy energy (to be replaced with true (ground) state
  // energy after reading data), other input data to header:
  h.Eg = std::numeric_limits<double>::max();
  BOOST_CHECK(h.Eg > 1e+308);
  h.ETol = input.ETol;
  BOOST_CHECK_EQUAL(h.ETol, 1e-5);
  h.CTol = input.CTol;
  BOOST_CHECK_EQUAL(h.CTol, 1e-3);
  h.omega0 = input.input_omega0;
  BOOST_CHECK_EQUAL(h.omega0, 0.5);
  h.B = input.input_B;
  BOOST_CHECK_EQUAL(h.B, 0.1);
  
  // Check each .dat file for degenerate (ground) states:
  for (size_t i = 0; i != input.datFiles.size(); ++i) {
    // Write .dat file name to header:
    h.file = input.datFiles[i];
    // Open .dat file:
    ifstream in(input.datFiles[i].c_str());
    // Search for states with the omega0/B values specified in input:
    if (in) {
      Find_States(in, input, h, eStates);
    }
    // Close .dat file:
    in.close();
  } // End loop over .dat files
  
  // The test ground state energy should be smallest energy found in test2.dat:
  BOOST_CHECK_EQUAL(h.Eg, eStates[0].energy);
  // Check that each eigenstate has the correct information (that hasn't just been
  // checked) stored with it:
  BOOST_CHECK_EQUAL(eStates.size(), five);
  BOOST_CHECK_EQUAL(h.Eg, eStates[0].energy); // eState 1 has lowest energy
  // First eigenstate:
  BOOST_CHECK_EQUAL(eStates[0].h.N, 3);
  BOOST_CHECK_EQUAL(eStates[0].h.Lz, 1);
  BOOST_CHECK_EQUAL(eStates[0].h.tS, 1);
  BOOST_CHECK_EQUAL(eStates[0].h.tSz, 1);
  BOOST_CHECK_EQUAL(eStates[0].h.Ez, 0.012734447);
  BOOST_CHECK_EQUAL(eStates[0].h.numOmega0, 1);
  BOOST_CHECK_EQUAL(eStates[0].h.omega0, 0.5);
  BOOST_CHECK_EQUAL(eStates[0].h.numB, 1);
  BOOST_CHECK_EQUAL(eStates[0].h.B, 0.1);
  BOOST_CHECK_EQUAL(eStates[0].h.MPsize, 5);
  BOOST_CHECK_EQUAL(eStates[0].h.numEigs, 13);
  BOOST_CHECK_EQUAL(eStates[0].h.ETol, 1e-5);
  BOOST_CHECK_EQUAL(eStates[0].h.CTol, 1e-3);
  BOOST_CHECK_EQUAL(eStates[0].h.file, "test2.dat");
  BOOST_CHECK_EQUAL(eStates[0].h.numDegenStates, 1); // first degen state in file
  BOOST_CHECK(closer(eStates[0].h.Eg, 3.09999999999111)); // global groundstate
  BOOST_CHECK_EQUAL(h.Eg, eStates[0].energy); // global groundstate
  BOOST_CHECK_EQUAL(eStates[0].energy, 3.09999999999111); // actual state energy
  BOOST_CHECK_EQUAL(eStates[0].linCombState.size(), five);
  BOOST_CHECK_EQUAL(eStates[0].h.sumOfSquares, 0);
  // Second eigenstate:
  BOOST_CHECK_EQUAL(eStates[1].h.N, 3);
  BOOST_CHECK_EQUAL(eStates[1].h.Lz, 1);
  BOOST_CHECK_EQUAL(eStates[1].h.tS, 1);
  BOOST_CHECK_EQUAL(eStates[1].h.tSz, 1);
  BOOST_CHECK_EQUAL(eStates[1].h.Ez, 0.012734447);
  BOOST_CHECK_EQUAL(eStates[1].h.numOmega0, 1);
  BOOST_CHECK_EQUAL(eStates[1].h.omega0, 0.5);
  BOOST_CHECK_EQUAL(eStates[1].h.numB, 1);
  BOOST_CHECK_EQUAL(eStates[1].h.B, 0.1);
  BOOST_CHECK_EQUAL(eStates[1].h.MPsize, 5);
  BOOST_CHECK_EQUAL(eStates[1].h.numEigs, 13);
  BOOST_CHECK_EQUAL(eStates[1].h.ETol, 1e-5);
  BOOST_CHECK_EQUAL(eStates[1].h.CTol, 1e-3);
  BOOST_CHECK_EQUAL(eStates[1].h.numDegenStates, 2); // 2nd degen state in file
  BOOST_CHECK(closer(eStates[1].h.Eg, 3.09999999999111)); // global grounstate
  BOOST_CHECK_EQUAL(eStates[1].h.file, "test2.dat");
  BOOST_CHECK_EQUAL(eStates[1].energy, 3.10); // actual state energy
  BOOST_CHECK_EQUAL(eStates[1].linCombState.size(), five);
  BOOST_CHECK_EQUAL(eStates[1].h.sumOfSquares, 0);
  // Third eigenstate:
  BOOST_CHECK_EQUAL(eStates[2].h.N, 3);
  BOOST_CHECK_EQUAL(eStates[2].h.Lz, 1);
  BOOST_CHECK_EQUAL(eStates[2].h.tS, 1);
  BOOST_CHECK_EQUAL(eStates[2].h.tSz, 1);
  BOOST_CHECK_EQUAL(eStates[2].h.Ez, 0.012734447);
  BOOST_CHECK_EQUAL(eStates[2].h.numOmega0, 1);
  BOOST_CHECK_EQUAL(eStates[2].h.omega0, 0.5);
  BOOST_CHECK_EQUAL(eStates[2].h.numB, 1);
  BOOST_CHECK_EQUAL(eStates[2].h.B, 0.1);
  BOOST_CHECK_EQUAL(eStates[2].h.MPsize, 5);
  BOOST_CHECK_EQUAL(eStates[2].h.numEigs, 13);
  BOOST_CHECK_EQUAL(eStates[2].h.ETol, 1e-5);
  BOOST_CHECK_EQUAL(eStates[2].h.CTol, 1e-3);
  BOOST_CHECK_EQUAL(eStates[2].h.numDegenStates, 3); // 3rd degen state in file
  BOOST_CHECK(closer(eStates[2].h.Eg, 3.09999999999111)); // global grounstate
  BOOST_CHECK_EQUAL(eStates[2].h.file, "test2.dat");
  BOOST_CHECK_EQUAL(eStates[2].energy, 3.100000456); // actual state energy
  BOOST_CHECK_EQUAL(eStates[2].linCombState.size(), four);
  BOOST_CHECK(close(eStates[2].h.sumOfSquares, 0.0009));
  // Fourth eigenstate:
  BOOST_CHECK_EQUAL(eStates[3].h.N, 3);
  BOOST_CHECK_EQUAL(eStates[3].h.Lz, 1);
  BOOST_CHECK_EQUAL(eStates[3].h.tS, 1);
  BOOST_CHECK_EQUAL(eStates[3].h.tSz, 1);
  BOOST_CHECK_EQUAL(eStates[3].h.Ez, 0.012734447);
  BOOST_CHECK_EQUAL(eStates[3].h.numOmega0, 3);
  BOOST_CHECK_EQUAL(eStates[3].h.omega0, 0.5);
  BOOST_CHECK_EQUAL(eStates[3].h.numB, 3);
  BOOST_CHECK_EQUAL(eStates[3].h.B, 0.1);
  BOOST_CHECK_EQUAL(eStates[3].h.MPsize, 2);
  BOOST_CHECK_EQUAL(eStates[3].h.numEigs, 3);
  BOOST_CHECK_EQUAL(eStates[3].h.ETol, 1e-5);
  BOOST_CHECK_EQUAL(eStates[3].h.CTol, 1e-3);
  BOOST_CHECK_EQUAL(eStates[3].h.numDegenStates, 1); // first degen state in new file
  BOOST_CHECK(closer(eStates[3].h.Eg, 3.09999999999111)); // global grounstate
  BOOST_CHECK_EQUAL(eStates[3].h.file, "test3.dat");
  BOOST_CHECK_EQUAL(eStates[3].energy, 3.1000000023); // actual state energy
  BOOST_CHECK_EQUAL(eStates[3].linCombState.size(), two);
  BOOST_CHECK(closer(eStates[3].h.sumOfSquares, 0));
  // Fifth eigenstate:
  BOOST_CHECK_EQUAL(eStates[4].h.N, 3);
  BOOST_CHECK_EQUAL(eStates[4].h.Lz, 1);
  BOOST_CHECK_EQUAL(eStates[4].h.tS, 1);
  BOOST_CHECK_EQUAL(eStates[4].h.tSz, 1);
  BOOST_CHECK_EQUAL(eStates[4].h.Ez, 0.012734447);
  BOOST_CHECK_EQUAL(eStates[4].h.numOmega0, 3);
  BOOST_CHECK_EQUAL(eStates[4].h.omega0, 0.5);
  BOOST_CHECK_EQUAL(eStates[4].h.numB, 3);
  BOOST_CHECK_EQUAL(eStates[4].h.B, 0.1);
  BOOST_CHECK_EQUAL(eStates[4].h.MPsize, 2);
  BOOST_CHECK_EQUAL(eStates[4].h.numEigs, 3);
  BOOST_CHECK_EQUAL(eStates[4].h.ETol, 1e-5);
  BOOST_CHECK_EQUAL(eStates[4].h.CTol, 1e-3);
  BOOST_CHECK_EQUAL(eStates[4].h.numDegenStates, 2); // 2nd degen state in new file
  BOOST_CHECK(closer(eStates[4].h.Eg, 3.09999999999111)); // global grounstate
  BOOST_CHECK_EQUAL(eStates[4].h.file, "test3.dat");
  BOOST_CHECK_EQUAL(eStates[4].energy, 3.1000045); // actual state energy
  BOOST_CHECK_EQUAL(eStates[4].linCombState.size(), two);
  BOOST_CHECK(closer(eStates[4].h.sumOfSquares, 0));

  return;
}

test_suite* init_unit_test_suite(int, char**)
{
  test_suite* test = BOOST_TEST_SUITE("twoPointSpin Test");
  test->add(BOOST_TEST_CASE(&twoPointTest));
  return test;
}
