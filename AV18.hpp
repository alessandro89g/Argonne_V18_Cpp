#ifndef AV18_HPP
#define AV18_HPP

#include "Row.hpp"
#include <cmath>
#define ushort unsigned short

/**********************   Argonne v18      ********************************
*   Argonne v18 and vn' and Super-Soft Core (C) potential package
*
*   prepared 1 Sep 1994 by R. B. Wiringa, Physics Division,
*   Argonne National Laboratory, Argonne, IL 60439
*   e-mail: wiringa@theory.phy.anl.gov
*
*   reference:
*    "Accurate nucleon-nucleon potential with charge-independence breaking"
*     R. B. Wiringa, V. G. J. Stoks, and R. Schiavilla,
*     Physical Review C51, 38 (1995) - WSS95.
*
*   option for v8' reprojection of v18 potential added 10 Jan 2001
*
*   reference:
*    "Quantum Monte Carlo calculations of nuclei with A<=7"
*     B. S. Pudliner, V. R. Pandharipande, J. Carlson, Steven C. Pieper,
*     and R. B. Wiringa, Physical Review C56, 1720 (1997) - PPCPW97
*
*   option for v6', v4', vx', v2', v1' potentials added 16 Jul 2002
*
*   reference:
*    "Evolution of Nuclear Spectra with Nuclear Forces"
*     R. B. Wiringa and Steven C. Pieper,
*     Physical Review Letters 89, 182501 (2002) - WP02
*
*   option for Super-Soft Core (C) added 14 Feb 2007
*
*   reference:
*    "Construction d'un potentiel nucleon-nucleon a coeur tres mou (SSC)"
*     R. de Tourreil and D. W. L. Sprung,
*     Nuclear Physics A201, 193 (1973) - TS73
*
*   option for modfied Argonne v8' and modified SSC   v8' added 14 Feb 2007
*     modifications selected by Steve Pieper
*     correction to modified SSC   v8' on 4 Apr 2007
*
*   this file contains 4 subroutines:
*     subroutine av18pw(lpot,l,s,j,t,t1z,t2z,r,vpw)
*     subroutine av18op(lpot,r,vnn)
*     subroutine empot(lpot,r,vem)
*     subroutine consts(lpot,hc,mpi0,mpic,mp,mn,alpha,mup,mun)
*
*   the variable lpot selects between v18, v8' and other options
*
*   av18pw gives the full potential in a particular partial wave
*   av18op gives the strong interaction part in operator format
*   empot  gives the electromagneti   part in operator format
*   consts gives values of fundamental constants and masses used
*
*   notes:
*   1) av18pw includes full EM interaction for lpot=1;
*      for lpot>1 it includes only C1(pp), i.e.,
*      Coulomb with a form factor for pp channels.
*
*   2) empot does not include the energy-dependence of the Coulomb
*      interaction used in eq.(4) of WSS95, i.e., it uses alpha,
*      not alpha'.
*
*   3) the vacuum polarization in empot is a short-range approximation
*      to eq.(7) suitable for bound states, but not for scattering.
*      it is taken from eq.(3.13) of Rev. Mod. Phys. 44, 48 (1972)
*
*      8/28/97 error in this formula detected and corrected:
*      should be -(gamma+5/6) instead of printed (-gamma+5/6)
*
*   4) these subroutines should be compiled with a compiler option
*      that forces all floating point constants to be evaluated at
*      real*8 significance, e.g., on an IBM RS6000 the xlf compiler
*      option qdpc=e should be used; on SGI machines, the -r8 option
*      should be used; on a Cray no action is needed.
*      if such an option is not available and the default precision is
*      real*4 (32 bits), then all constants should be explicitly
*      converted to double precision by appending a D0.
*
*   5) consts now (14 Feb 2007) depend upon potential:
*      need to call to generate appropriate hbar**2/M
*/

//  C++ version thanks to Alessandro Grassi

/*  *id* av18pw **********************************************************
*   subroutine for partial-wave projection of argonne v18 potential
*   or super-soft core (C) v14 potential
*   or reprojected vn' potential
*   calls subroutines av18op, empot
*   ----------------------------------------------------------------------
*   arguments for av18pw
*   lpot: switch for potential choice
*       -----------------------------------------------
*           Argonne                Super-Soft Core (C)
*         = 1 : av18              = 101 : ssc   v14
*         = 2 : av8'              = 102 : ssc   v8'
*         = 3 : av6'
*         = 4 : av4'
*         = 5 : avx'
*         = 6 : av2'
*         = 7 : av1'
*         = 8 : modified av8'     = 108 : modified ssc   v8'
*       -----------------------------------------------
*   l:    orbital angular momentum of pair (0,1,2,...)
*   s:    total spin of pair (0 or 1)
*   j:    total angular momentum of pair (0,1,2,...)
*   t:    total isospin of pair (0 or 1)
*   t1z:  isospin of particle 1 (1 for p, -1 for n)
*   t2z:     "    "     "     2 (1 for p, -1 for n)
*   r:    separation in fm
*   v:    returned potential in MeV (2x2 array)
*         (includes all strong and em terms)
*   ----------------------------------------------------------------------
*   order of terms in v(l,m):
*        single channel                 coupled channel (l=j-1,s=1)
*        v(1,1) = v(l,s,j,t,t1z,t2z)    v(1,1) = v(l,s,j,t,t1z,t2z)
*        v(2,1) = 0                     v(2,1) = v(l<->l+2)
*        v(1,2) = 0                     v(1,2) = v(l<->l+2)
*        v(2,2) = 0                     v(2,2) = v(l+2,s,j,t,t1z,t2z)
*   ----------------------------------------------------------------------
*/
void av18pw(const ushort lpot, const ushort l, const ushort s, const ushort  j,
            const ushort t, const short t1z, const short t2z, double r, Matrix<double,2>& vpw);

/* *id* av18op **********************************************************
* subroutine for strong interaction part of argonne v18 potential
* or super-soft core (C) v14 potential
* or reprojected vn' potential in operator format
* calls subroutine consts
* ----------------------------------------------------------------------
* arguments for av18pot
* lpot: switch for potential choice
*     -----------------------------------------------
*         Argonne                Super-Soft Core (C)
*       = 1 : av18              = 101 : sscc v14
*       = 2 : av8'              = 102 : sscc v8'
*       = 3 : av6'
*       = 4 : av4'
*       = 5 : avx'
*       = 6 : av2'
*       = 7 : av1'
*       = 8 : modified av8'     = 108 : modified sscc v8'
*     -----------------------------------------------
* r:    separation in fm
* vnn:  output potential in MeV (18 component array)
* ----------------------------------------------------------------------
* order of operators l in vnn(l):
* l:    1=1                              2=t1.t2
*       3=s1.s2                          4=(s1.s2)(t1.t2)
*       5=S12 [=3(s1.r)(s2.r)-s1.s2]     6=S12(t1.t2)
*       7=L.S                            8=L.S(t1.t2)
*       9=L**2                          10=L**2(t1.t2)
*      11=L**2(s1.s2)                   12=L**2(s1.s2)(t1.t2)

*      13=(L.S)**2                      14=(L.S)**2(t1.t2)
*      15=T12 [=3*t1z*t2z-t1.t2]        16=(s1.s2)T12
*      17=S12*T12                       18=t1z+t2z
* where s1=sigma_1, t1=tau_1, t1z=tau_1(z), etc.
* --------------------------------------------------------------------*/
void av18op(const ushort lpot, double& r, Row<double,18>& vnn);

/**id* empot ***********************************************************
* subroutine for electromagnetic part of Argonne v18 potential
* for avn' models returns pp Coulomb only
* calls subroutine consts
* ----------------------------------------------------------------------
* arguments for empot
* lpot: switch for potential choice
*       = 1 : full EM potential
*       > 1 : C1(pp) only
* r:    input separation in fm
* vem:  output potential in MeV (14 component array)
* ----------------------------------------------------------------------
* order of operators in vem(l)
* l:    1=C1    (pp)          2=DF    (pp)          3=C2      (pp)
*       4=VP    (pp)                                5=C1      (np)
*       6=s1.s2 (pp)          7=s1.s2 (nn)          8=s1.s2   (np)
*       9=S12   (pp)         10=S12   (nn)         11=S12     (np)
*      12=L.S   (pp)         13=L.S   (nn)         14=L.S     (np)
* C1 = one-photon-exchange Coulomb with form factor
* C2 = two-photon-exchange Coulomb
* DF = Darwin-Foldy
* VP = vacuum polarization (short-range approximation)
* all other terms from magnetic moment (MM) interactions
* --------------------------------------------------------------------*/
void empot(const ushort lpot, const double r, Row<double,14>& vem);

/**id* consts **********************************************************
* subroutine for constants in av18 and sscc potentials
* ----------------------------------------------------------------------
* arguments for consts
* lpot:  input potential
* hc:    output value for hbar*c (MeV-fm)
* mpi0:    "      "    "  neutral pion mass (MeV)
* mpic:    "      "    "  charged pion mass (MeV)
* mp:      "      "    "  proton mass (MeV)
* mn:      "      "    "  neutron mass (MeV)
* alpha:   "      "    "  electromagnetic constant alpha
* mup:     "      "    "  proton magnetic moment (nm)
* mun:     "      "    "  neutron magnetic moment (nm)
* --------------------------------------------------------------------*/
void consts(ushort lpot, double& hc, double& mpi0, double& mpic, double& mp,
            double& mn, double& alpha, double& mup, double& mun);

#endif // AV18_HPP
