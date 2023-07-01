#include "AV18.hpp"
#include <functional>


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
void av18pw(const ushort lpot, const ushort l, const ushort s, const ushort  j, const ushort t,
            const short t1z, const short t2z, double r, double** vpw){

    int s1ds2,t1dt2,t12;
    double vnn[18];
    double vem[14];
/* ------------------------
   strong interaction terms
-------------------------*/
    av18op(lpot,r,vnn);
    s1ds2=4*s-3;
    t1dt2=4*t-3;
    t12=3*t1z*t2z-t1dt2;
    double vc=vnn[0]+t1dt2*vnn[1]+s1ds2*vnn[2]+s1ds2*t1dt2*vnn[3]
         +t12*vnn[14]+s1ds2*t12*vnn[15]+(t1z+t2z)*vnn[17];
    double vt=vnn[4]+t1dt2*vnn[5]+t12*vnn[16];
    double vls=vnn[6]+t1dt2*vnn[7];
    double vl2=vnn[8]+t1dt2*vnn[9]+s1ds2*vnn[10]+s1ds2*t1dt2*vnn[11];
    double vls2=vnn[12]+t1dt2*vnn[13];
/*  ---------------------
 *  electromagnetic terms
 *  -------------------*/
    empot(lpot,r,vem);

    if (t1z+t2z<0)
        goto l10;
    else if (t1z+t2z==0)
        goto l20;
    else
        goto l30;
l10:
    vc=vc+s1ds2*vem[6];
    vt=vt+vem[9];
    goto l40;
l20:
    vc=vc+vem[4]+s1ds2*vem[7];
    vt=vt+vem[10];
    vls=vls+vem[13];
    goto l40;
l30:
    vc=vc+vem[0]+vem[1]+vem[2]+vem[3]+s1ds2*vem[5];
    vt=vt+vem[8];
    vls=vls+vem[11];
l40:
// ---------------------
    int ncc, ls, lsm, lsp;
    double s12, s12m, s12p;

    ncc=1;
    if (s==1 && j > l)
        ncc=2;
    if (ncc == 1) {
        s12=0.;
        if (s == 1 && l == j)
            s12=2.;
        if (l == (j+1)) s12=-2.*(j+2.)/(2.*j+1.);
        ls=(j*(j+1)-l*(l+1)-s*(s+1))/2;
        vpw[0][0]=vc+s12*vt+ls*vls+l*(l+1)*vl2+ls*ls*vls2;
        vpw[1][0]=0;
        vpw[0][1]=0;
        vpw[1][1]=0;
    } else if (ncc == 2) {
        s12m=-2.*(j-1.)/(2.*j+1.);
        s12=std::sqrt(36.*j*(j+1))/(2.*j+1.);
        s12p=-2.*(j+2.)/(2.*j+1.);
        lsm=j-1;
        lsp=-(j+2);
        vpw[0][0] = vc+s12m*vt+lsm*vls+l*(l+1)*vl2+lsm*lsm*vls2;
        vpw[1][0] = s12*vt;
        vpw[0][1] = s12*vt;
        vpw[1][1] = vc+s12p*vt+lsp*vls+(l+2)*(l+3)*vl2+lsp*lsp*vls2;
    }
}




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
void av18op(const ushort lpot, double& r, double* vnn) {
    double mpi0,mpic,mp,mn,mup,mun;
    double mpi,mu0,muc,mu;

    const double small = 1.E-4;
    const double vsmall = 1.E-10;

    double fsq, cpi, rws, aiws, x, x0, xc;
    double tpi, ypi0, tpi0, ypic, tpic;
    double tpi2, ws, ws0, wsp, wsx, wsx2;
    double dypi00, dypic0, ypi0p, ypicp, ypi;
    double p11pp, p11np, p11nn, pt1pp, pt1np, pt1nn, pls1, pl211, pls21,
        p10, pt0,pls0, pl210, pls20, p01pp, p01np, p01nn, pl201, p00, pl200,
        p11, p11cd, p11cs, pt1, pt1cd, pt1cs, p01, p01cd, p01cs;
    double rr4, rc4, rc6, hr, pq0, alpha, rcut, hc;
/*-------------------
* statement functions
* -----------------*/
    std::function<double(double,double)> yc =
        [](const double t, const double x) {return exp(-t)/x;};

    std::function<double(double,double)> yt =
        [](const double t, const double x) {return (1+3/t+3/std::pow(t,2))*exp(-t)/x;};

    std::function<double(double,double)> yls =
        [](const double t, const double x) {return -(1+t)*exp(-t)/std::pow(x,3);};

    std::function<double(double,double)> yl2 =
        [](const double t, const double x) {return (1+2/t)*exp(-t)/std::pow(x,3);};
//-------------------
    for (int l=0; l<18; l++)
        vnn[l] = 0;
/* ---------------------------------
* argonne potential and derivatives
*  -------------------------------*/
    if (lpot < 100) {
        consts(lpot,hc,mpi0,mpic,mp,mn,alpha,mup,mun);
        mpi=(mpi0+2.*mpic)/3.;
        mu0=mpi0/hc;
        muc=mpic/hc;
        mu=mpi/hc;
        fsq=.075;
        cpi=2.1;
        rws=.5;
        aiws=5.;
        x=mu*r;
        x0=mu0*r;
        xc=muc*r;
        if (r <= small) {
            tpi=3*std::pow(cpi,2)*r/std::pow(mu,3);
            ypi0=std::pow(mpi0/mpic,2)*(mpi0/3)*cpi*r/mu0;
            tpi0=3*cpi*ypi0/std::pow(mu0,2);
            ypic=(mpic/3)*cpi*r/muc;
            tpic=3*cpi*ypic/std::pow(muc,2);
        } else{
            rcut=1-exp(-cpi*r*r);
            ypi=exp(-x)*rcut/x;
            tpi=(1+(3+3/x)/x)*ypi*rcut;
            ypi0=std::pow(mpi0/mpic,2)*(mpi0/3)*exp(-x0)*rcut/x0;
            tpi0=(1+(3+3/x0)/x0)*ypi0*rcut;
            ypic=(mpic/3)*exp(-xc)*rcut/xc;
            tpic=(1+(3+3/xc)/xc)*ypic*rcut;
        }
        ypi0=fsq*ypi0;
        ypic=fsq*ypic;
        tpi0=fsq*tpi0;
        tpic=fsq*tpic;
        tpi2=tpi*tpi;
        ws=1/(1+exp((r-rws)*aiws));
        ws0=1/(1+exp(-rws*aiws));
        wsp=ws*(1+aiws*exp(-rws*aiws)*ws0*r);
        wsx=ws*x;
        wsx2=wsx*x;
        dypi00=std::pow(mpi0/mpic,2)*(mpi0/3)*cpi/mu0;
        dypic0=(mpic/3)*cpi/muc;
        ypi0p=ypi0-fsq*dypi00*ws*r/ws0;
        ypicp=ypic-fsq*dypic0*ws*r/ws0;
        ypi=(ypi0+2*ypic)/3;
        tpi=(tpi0+2*tpic)/3;
        p11pp=  -7.62701*tpi2+1815.4920*wsp+1847.8059*wsx2+ypi0p;
        p11np=  -7.62701*tpi2+1813.5315*wsp+1847.8059*wsx2-ypi0p+2*ypicp;
        p11nn=  -7.62701*tpi2+1811.5710*wsp+1847.8059*wsx2+ypi0p;
        pt1pp=   1.07985*tpi2 -190.0949*wsx -811.2040*wsx2+tpi0;
        pt1np=   1.07985*tpi2 -190.0949*wsx -811.2040*wsx2-tpi0+2*tpic;
        pt1nn=   1.07985*tpi2 -190.0949*wsx -811.2040*wsx2+tpi0;
        pls1=    -.62697*tpi2 -570.5571*wsp +819.1222*wsx2;
        pl211=    .06709*tpi2 +342.0669*wsp -615.2339*wsx2;
        pls21=    .74129*tpi2   +9.3418*wsp -376.4384*wsx2;
        p10=    -8.62770*tpi2+2605.2682*wsp +441.9733*wsx2-ypi0p-2*ypicp;
        pt0=    1.485601*tpi2-1126.8359*wsx +370.1324*wsx2-tpi0-2*tpic;
        pls0=     .10180*tpi2  +86.0658*wsp -356.5175*wsx2;
        pl210=   -.13201*tpi2 +253.4350*wsp   -1.0076*wsx2;
        pls20=    .07357*tpi2 -217.5791*wsp  +18.3935*wsx2;
        p01pp= -11.27028*tpi2+3346.6874*wsp-3*ypi0p;
        p01np= -10.66788*tpi2+3126.5542*wsp-3*(-ypi0p+2*ypicp);
        p01nn= -11.27028*tpi2+3342.7664*wsp-3*ypi0p;
        pl201=    .12472*tpi2  +16.7780*wsp;
        p00=    -2.09971*tpi2+1204.4301*wsp-3*(-ypi0p-2*ypicp);
        pl200=   -.31452*tpi2 +217.4559*wsp;
        p11=(p11pp+p11nn+p11np)/3;
        p11cd=(.5*(p11pp+p11nn)-p11np)/6;
        p11cs=(p11pp-p11nn)/4;
        pt1=(pt1pp+pt1nn+pt1np)/3;
        pt1cd=(.5*(pt1pp+pt1nn)-pt1np)/6;
        pt1cs=(pt1pp-pt1nn)/4;
        p01=(p01pp+p01nn+p01np)/3;
        p01cd=(.5*(p01pp+p01nn)-p01np)/6;
        p01cs=(p01pp-p01nn)/4;
/*------------------------
* option for v8' reduction
* ----------------------*/
        if (lpot >= 2) {
            p00=p00+2*pl200;
            p11=p11+2*pl211+4*pls21/3;
            pt1=pt1-5*pls21/12;
            pls1=pls1-.5*pls21;
            pls0=pls0-2*pl210-3*pls20;
        }
/*------------------------
* option for v6' redutcion
* ----------------------*/
        if (lpot >= 3  &&  lpot <= 7) p10=p10-.3*pls0;
/*------------------------
* option for v4' reduction
* ----------------------*/
        if (lpot >= 4  &&  lpot <= 7) p10=p10+.8735*pt0;
/*------------------------
* option for vx' reduction
* ----------------------*/
        if (lpot == 5) {
            vnn[0]=.0625*(9*p11+3*p10+3*p01+p00);
            vnn[1]=.0125*(9*p11-5*p10-5*p01+p00);
            vnn[2]=vnn[1];
            vnn[3]=vnn[1];
            return;
        }
/*------------------------
* option for v2' reduction
* ----------------------*/
        else if (lpot == 6) {
            vnn[0]=.25*(3*p01+p10);
            vnn[1]=.25*(  p01-p10);
            return;
        }
/*------------------------
* option for v1' reduction
* ----------------------*/
        else if (lpot == 7) {
            vnn[0]=.5*(p01+p10);
            return;
        }
/*------------------------
* option for modified v8'
* ----------------------*/
        if (lpot == 8) p11=p11-.37*tpi2;
//------------------------
        vnn[0]=.0625*(9*p11+3*p10+3*p01+p00);
        vnn[1]=.0625*(3*p11-3*p10  +p01-p00);
        vnn[2]=.0625*(3*p11  +p10-3*p01-p00);
        vnn[3]=.0625*(  p11  -p10  -p01+p00);
//------------------------
        if (lpot == 4) return;
//------------------------
        vnn[4]=.25*(3*pt1+pt0);
        vnn[5]=.25*(  pt1-pt0);
//------------------------
        if (lpot == 3) return;
//------------------------
        vnn[6]=.25*(3*pls1+pls0);
        vnn[7]=.25*(  pls1-pls0);
//------------------------
        if (lpot == 2  ||  lpot == 8) return;
//------------------------
        vnn[8]=.0625*(9*pl211+3*pl210+3*pl201+pl200);
        vnn[9]=.0625*(3*pl211-3*pl210+  pl201-pl200);
        vnn[10]=.0625*(3*pl211+  pl210-3*pl201-pl200);
        vnn[11]=.0625*(  pl211-  pl210-  pl201+pl200);
        vnn[12]=.25*(3*pls21+pls20);
        vnn[13]=.25*(  pls21-pls20);
        vnn[14]=.25*(3*p11cd+p01cd);
        vnn[15]=.25*(  p11cd-p01cd);
        vnn[16]=pt1cd;
        vnn[17]=p01cs;
        return;
    }
/*---------------------------------------------
* super-soft core (C) potential and derivatives
* -------------------------------------------*/
    else if (lpot > 100) {
        if (r <= vsmall) r=vsmall;
        x=.7*r;
        rr4=std::pow(r,4);
        rc4=1-exp(-rr4);
        rc6=1-exp(-std::pow(r,6));
        hr=10.463;
        p11=144.83*exp(-rr4/std::pow(0.88787,2))+(-241.34*yc(3.3788*x,x)+(hr/3)*yc(x,x))*rc4;
        p10=215.32*exp(-rr4/std::pow(0.85807,2))+(-883.6*yc(3.5042*x,x)-hr*yc(x,x))*rc4;
        p01=375.*exp(-rr4/std::pow(0.47552,2))+(-1001.6*yc(3.6071*x,x)-hr*yc(x,x))*rc4;
        p00=75.653*exp(-rr4/std::pow(3.0000,2))+(-286.26*yc(2.0254*x,x)+3*hr*yc(x,x))*rc4;
        pt1=36.*exp(-rr4/std::pow(1.0805,2))+(-110.*yt(3.9529*x,x)+(hr/3)*yt(x,x))*rc6;
        pt0=-58.951*exp(-rr4/std::pow(1.3171,2))+(395.18*yt(4.3098*x,x)-hr*yt(x,x))*rc6;
        pls1=(520.*yls(5.661*x,x)-54.85*yls(4.0141*x,x))*rc6;
        pls0=(-40.466*yls(5.768*x,x)-40.408*yls(4.0676*x,x))*rc6;
        pl211=(6.65*yl2(1.965*x,x)-.959*yl2(x,x))*rc6;
        pl210=(17.626*yl2(2.6463*x,x)-.35261*yl2(x,x))*rc6;
        pl201=(14.*yl2(2.5*x,x)-.35*yl2(x,x))*rc6;
        pl200=(15.633*yl2(2.01*x,x)+.72581*yl2(x,x))*rc6;
        pq0=-3.9904*yl2(2.4583*x,x)*rc6;
/*------------------------
* option for v8' reduction
* ----------------------*/
        if (lpot >= 102) {
            p00=p00+2*pl200;
            pls0=pls0-2*pl210-10*pq0;
            p11=p11+2*pl211;
        }
/*------------------------
* option for modified v8'
* ----------------------*/
        if (lpot == 108) p11=p11-111*yc(3.3788*x,x)*rc4;
//------------------------
        vnn[0]=.0625*(9*p11+3*p10+3*p01+p00);
        vnn[1]=.0625*(3*p11-3*p10+  p01-p00);
        vnn[2]=.0625*(3*p11+  p10-3*p01-p00);
        vnn[3]=.0625*(  p11-  p10-  p01+p00);
        vnn[4]=.25*(3*pt1+pt0);
        vnn[5]=.25*(  pt1-pt0);
        vnn[6]=.25*(3*pls1+pls0)+.75*pq0;
        vnn[7]=.25*(  pls1-pls0)-.75*pq0;
//------------------------
        if (lpot >= 102) return;
//------------------------
        vnn[8]= .0625*(9*pl211+3*pl210+3*pl201+pl200)-.75*pq0;
        vnn[9]=.0625*(3*pl211-3*pl210+  pl201-pl200)+.75*pq0;
        vnn[10]=.0625*(3*pl211+  pl210-3*pl201-pl200)-.25*pq0;
        vnn[11]=.0625*(  pl211-  pl210-  pl201+pl200)+.25*pq0;
        vnn[12]=1.5*pq0;
        vnn[13]=-1.5*pq0;

    }
}




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
void empot(const ushort lpot, const double r, double* vem) {
    double mpi0,mpic,mp,mn,mup,mun;
    double kr,me,mr;
    double alpha, b, br, pi, gamma, beta, fcoulr, ftr3, flsr3, hc,
        fivp, fdelta, fnpr;

    const double small = 1.E-5;

    for (int l=0; l<14; l++) {
        vem[l] = 0;;
    }

    consts(lpot,hc,mpi0,mpic,mp,mn,alpha,mup,mun);
    b=4.27;
    br=b*r;
    pi=acos(-1.);
    me=0.510999;
    mr=mp*mn/(mp+mn);
    gamma=0.577216;
    beta=.0189;
    if (r < small) {
        fcoulr=5*b/16;
        ftr3=std::pow(b,3)*std::pow(br,2)/720;
        flsr3=std::pow(b,3)/48;
        kr=me*small/hc;
    }
    else {
        fcoulr=(1-(1+11*br/16+3*std::pow(br,2)/16+std::pow(br,3)/48)*exp(-br))/r;
        ftr3=(1-(1+br+std::pow(br,2)/2+std::pow(br,3)/6+std::pow(br,4)/24+std::pow(br,5)/144)*exp(-br))/std::pow(r,3);
        flsr3=(1-(1+br+std::pow(br,2)/2+7*std::pow(br,3)/48+std::pow(br,4)/48)*exp(-br))/std::pow(r,3);
        kr=me*r/hc;
    }
    fivp=-gamma-5./6.+abs(log(kr))+6*pi*kr/8;
    fdelta=std::pow(b,3)*(1+br+std::pow(br,2)/3)*exp(-br)/16;
    fnpr=std::pow(b,3)*(15+15*br+6*std::pow(br,2)+std::pow(br,3))*exp(-br)/384;
    vem[0]=alpha*hc*fcoulr;
    //------------------------
    if (lpot >= 2) return;
    //------------------------
    vem[1]=-alpha*std::pow(hc,3)*fdelta/(4*std::pow(mp,2));
    vem[2]=-std::pow(vem[0],2)/mp;
    vem[3]=2*alpha*vem[0]*fivp/(3*pi);
    vem[4]=alpha*hc*beta*fnpr;
    vem[5]=-alpha*std::pow(hc,3)*std::pow(mup,2)*fdelta/(6*std::pow(mp,2));
    vem[6]=-alpha*std::pow(hc,3)*std::pow(mun,2)*fdelta/(6*std::pow(mn,2));
    vem[7]=-alpha*std::pow(hc,3)*mup*mun*fdelta/(6*mn*mp);
    vem[8]=-alpha*std::pow(hc,3)*std::pow(mup,2)*ftr3/(4*std::pow(mp,2));
    vem[9]=-alpha*std::pow(hc,3)*std::pow(mun,2)*ftr3/(4*std::pow(mn,2));
    vem[10]=-alpha*std::pow(hc,3)*mup*mun*ftr3/(4*mp*mn);
    vem[11]=-alpha*std::pow(hc,3)*(4*mup-1)*flsr3/(2*std::pow(mp,2));
    vem[12]=0;
    vem[13]=-alpha*std::pow(hc,3)*mun*flsr3/(2*mn*mr);
}





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
void consts(const ushort lpot,double& hc,double& mpi0, double& mpic, double& mp,
            double& mn,double& alpha, double& mup, double& mun) {
    hc=197.327053;
    if (lpot < 100) {
        mpi0=134.9739;
        mpic=139.5675;
        mp=938.27231;
        mn=939.56563;
    } else if (lpot > 100) {
        mpi0=.7*hc;
        mpic=.7*hc;
        mp=std::pow(hc,2)/41.47;
        mn=std::pow(hc,2)/41.47;
    }
    alpha=1./137.035989;
    mup= 2.7928474;
    mun=-1.9130427;
}
