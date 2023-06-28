//
//  main.cpp
//  BDEI
//
//  Created by Frédéric Hecht on 10/04/2020.
//  mail: Frederic.Hecht@sorbonne-universite.fr 
//  Copyright © 2020 Frédéric Hecht. All rights reserved.

#include "BDEI.hpp"
//
// khi-2  95%  : 3.84, 5.99, 7.81 , 9.49
double khi_2_95[] = {0., 3.84, 5.99, 7.81, 9.49};// dat wiki https://fr.wikipedia.org/wiki/Loi_du_χ²


// following python logging values
const int debugVal = 3;
const int infoVal = 2;
const int warningVal = 1;
const int errorVal = 0;

int debug = errorVal;

int size_pool = std::thread::hardware_concurrency();
std::mt19937_64 rand_gen;

const char *algo_name[] = {"NEWUOA_BOUND", "MMA", "LBFGS", "DIRECT", 0};
const char *BDEIParams[] = {"mu", "lambda", "psi", "p"};

/*
 ODE
 Ue_t = mu( Ue - Ui)
 Ui_t = (la+psi) Ui - la Ui Ue -psi(1-p)

 equation 1 origine
 Ue_t = mu( Ue - Ui)
 Ui_t = c1 Ui - la Ui Ue -c0

 equation retrograde
 Ue_t = -mu( Ue - Ui)
 Ui_t = -c1 Ui + la Ui Ue +c0
 on part de 1,1

 de plus donc dans l'equation retrograde
 1 >= Ue > Ui. >0

 equation retrograde
 donc si Ue_t <0 => Ue decroissant


 comme Uc est la seule solution positive de -c1 Uc + la UcUc +c0 =0
 donc  -c1 Ui + la Ui Ui +c0 <0
 -c1 Ue + la Ue Ue +c0 <0


 donc Ue  est decroissant et pour Ui je ne sais pas


 // J'aimerais avoir etude a l'infini de Ve = Ue - Uc, Vi = Ui - Uc
 // je suppose que Ve->0 et Vi ->0 =>   Ve Vi est negligable

 Ve = Ue - Uc, Vi = Ui - Uc
 Ve_t = -mu (Ve -Vi)
 Vi_t = -c1 (Vi+Uc) + la (Vi+Uc) (Ve+Uc) -c0
 => Vi_t = -c1 Vi + la Uc ( Vi+Ve) + la Vi Ve

 donc ODE Lineaire est:
 a l'infini en negligant  Vi Ve
 We_t = -mu (We -Wi)
 Wi_t = -c1 Wi + la Uc ( Wi+We)

 W_t =  M W   M = [ [-mu, mu ], [la Uc, -c1 + la Uc ]] ;
 */

template<class R>
struct Diff {
    R val, dval;

    Diff(R x, R dx) : val(x), dval(dx) {}

    Diff(R x) : val(x), dval() {}

    Diff() : val(), dval() {}
};


template<class R>
ostream &operator<<(ostream &f, const Diff<R> &a) {
    f << a.val << " ( d = " << a.dval << ") ";
    return f;
}

template<class R>
Diff<R> &operator+=(Diff<R> &a, const Diff<R> &b) {
    a.val += b.val;
    a.dval += b.dval;
    return a;
}

template<class R>
Diff<R> operator+(const Diff<R> &a, const Diff<R> &b) { return Diff<R>(a.val + b.val, a.dval + b.dval); }

template<class R>
Diff<R> operator+(const R &a, const Diff<R> &b) { return Diff<R>(a + b.val, b.dval); }

template<class R>
Diff<R> operator*(const Diff<R> &a, const Diff<R> &b) {
    return Diff<R>(a.val * b.val, a.dval * b.val + a.val * b.dval);
}

template<class R>
Diff<R> operator*(double a, const Diff<R> &b) { return Diff<R>(a * b.val, a * b.dval); }

template<class R>
Diff<R> operator*(const Diff<R> &a, const double b) { return Diff<R>(a.val * b, a.dval * b); }

template<class R>
Diff<R> operator/(const Diff<R> &a, const Diff<R> &b) {
    return Diff<R>(a.val / b.val, (a.dval * b.val - a.val * b.dval) / (b.val * b.val));
}

template<class R>
Diff<R> operator/(const R &a, const Diff<R> &b) { return Diff<R>(a / b.val, (-a * b.dval) / (b.val * b.val)); }

template<class R>
Diff<R> operator/(const Diff<R> &a, const double &b) { return Diff<R>(a.val / b, a.dval / b); }

template<class R>
Diff<R> operator-(const Diff<R> &a, const Diff<R> &b) { return Diff<R>(a.val - b.val, a.dval - b.dval); }

template<class R>
Diff<R> operator-(const Diff<R> &a) { return Diff<R>(-a.val, -a.dval); }

template<class R>
Diff<R> sin(const Diff<R> &a) { return Diff<R>(sin(a.val), a.dval * cos(a.val)); }

template<class R>
Diff<R> cos(const Diff<R> &a) { return Diff<R>(cos(a.val), -a.dval * sin(a.val)); }

template<class R>
Diff<R> exp(const Diff<R> &a) { return Diff<R>(exp(a.val), a.dval * exp(a.val)); }

template<class R>
Diff<R> log(const Diff<R> &a) { return Diff<R>(log(a.val), a.dval / (a.val)); }

template<class R>
Diff<R> fabs(const Diff<R> &a) { return (a.val > R()) ? a : -a; }

template<class R>
Diff<R> max(const Diff<R> &a, const Diff<R> &b) { return (a.val > b.val) ? a : b; }

template<class R>
Diff<R> min(const Diff<R> &a, const Diff<R> &b) { return (a.val < b.val) ? a : b; }

template<class R>
Diff<R> abs(const Diff<R> &a) { return (a.val > R()) ? a : -a; }

template<class R>
bool operator<(const Diff<R> &a, const Diff<R> &b) { return a.val < b.val; }

template<class R>
bool operator<(const Diff<R> &a, double b) { return a.val < b; }

template<class R>
bool operator>(const Diff<R> &a, const Diff<R> &b) { return a.val > b.val; }

template<class R>
bool operator>(const Diff<R> &a, const R &b) { return a.val > b; }

template<class R>
bool operator<=(const Diff<R> &a, const Diff<R> &b) { return a.val <= b.val; }

template<class R>
bool operator<=(const Diff<R> &a, double b) { return a.val <= b; }

template<class R>
bool operator>=(const Diff<R> &a, const Diff<R> &b) { return a.val >= b.val; }

template<class R>
bool operator>=(const Diff<R> &a, const R &b) { return a.val >= b; }

template<class R>
Diff<R> sqrt(Diff<R> x) { return Diff<R>(sqrt(x.val), 0.5 * x.dval / sqrt(x.val)); }

typedef double R;
typedef double RR;
//typedef Diff<double> RR;// pout la diff Auto

R val(R x) { return x; }

R val(Diff<double> x) { return x.val; }

template<class R> R R2x2[2][2];

template<class R>
R inv(R P[][2], R P1[][2]) {
    R det = P[0][0] * P[1][1] - P[1][0] * P[0][1];
    P1[0][0] = P[1][1] / det;
    P1[0][1] = -P[0][1] / det;
    P1[1][0] = -P[1][0] / det;
    P1[1][1] = P[0][0] / det;
    return det;
}

template<class R>
void Mult(R A[][2], R x[], R Ax[]) {
    Ax[0] = A[0][0] * x[0] + A[0][1] * x[1];
    Ax[1] = A[1][0] * x[0] + A[1][1] * x[1];
}

template<class R>
void Mult(R A[][2], R B[][2], R AB[][2]) {
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            AB[i][j] = 0.;
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            for (int k = 0; k < 2; ++k) {
                AB[i][k] += A[i][j] * B[j][k];
            }
}

template<typename R>
ostream &operator<<(ostream &cout, const vector<R> &v) {
    size_t n = v.size() - 1;;
    for (int i = 0; i < n; ++i)
        cout << v[i] << " ; ";
    cout << v[n];
    return cout;
}

static const string emptystr = "";

template<class RR>
class BDEI_pb {
public:
    //    order 4 data:   mu,la,psi,p;
    typedef double R;
    const RR eps;
    const RR mu, la, psi, p;
    const RR c1, c0;// data modele
    RR a0, a1;
    RR U[2], Uc[2];
    RR uc, lp, lm;
    RR P[2][2], P1[2][2];
    RR piI, piE;
    int n;
    vector<RR> y;

    BDEI_pb(const RR *a, const R pie)
            : eps(1e-7),
              mu(max(a[0], eps)), la(max(a[1], eps)), psi(max(a[2], eps)), p(min(max(a[3], eps), RR(0.99999999))),
              c1(la + psi), c0(psi * (RR(1.) - p)), n(100000), y(3 * (n + 1)), piE(pie) { init(); }

    void init() {
        if (val(piE) < 0.) {
            const RR mu_plus_psi = mu + psi;
            if (val(la) == val(psi)) {
                piI = mu / mu_plus_psi;
            } else {
                const RR two_la_minus_psi = (val(la) == val(psi)) ? RR(1.): (RR(2.) * (la - psi));
                const RR det = sqrt(mu_plus_psi * mu_plus_psi + RR(2.) * mu * two_la_minus_psi);
                const RR piPlus = (det - mu_plus_psi) / two_la_minus_psi;
                piI = (0. <= val(piPlus) <= 1.) ? piPlus : ((-det - mu_plus_psi) / two_la_minus_psi);
            }
            piE = RR(1.) - piI;
        } else {
            piI = RR(1.) - piE;
        }

//        if (debug >= debugVal)
//            cout << setprecision(10) << " ** mu = " << mu << " , lambda =  " << la
//                 << " , psi = " << psi << " , p = " << p
//                 << " , freqs = " << piE << " , " << piI << endl;
//        assert(piI >= 0);
//        assert(piI <= 1.);
        /*
         F[0] = -mu*(y[0]-y[1]);
         F[1]= -c1*y[1] + la*y[0]*y[1] +c0;
         */
        //  la uc^2 - c1 uc + c0 =0
        RR delta = (c1) * (c1) - 4 * la * c0;
        // assert(delta>=R(0.));
        uc = (c1 - sqrt(delta)) / (2 * la);
        Uc[0] = Uc[1] = uc;
        // [ [-mu, mu ], [la Uc,-c1 +  la Uc ]]
        // M =[[ a,b] [c,d]]
        // lambda^± =  1/2(a+d ± sqrt(delta))
        // delta = (a-d)^2 + 4 bc
        // a = -mu, d =  -c1+la Uc
        // b = mu , c =  la Uc
        //  delta = (-mu - la Uc)^2 + 4 mu (-c1+la Uc)
        // 1/2 ( -mu+la Uc  ± sqrt( delta)
        RR a = -mu, d = -c1 + la * uc;
        RR b = mu, c = la * uc;
        RR delta1 = (a - d) * (a - d) + 4. * b * c;
        lm = 0.5 * (a + d - sqrt(delta1));// lm
        lp = 0.5 * (a + d + sqrt(delta1));//lp
        //  P = [ [b,b],[a-lp, a- lm] ]
        P[0][0] = b;
        P[0][1] = b;
        P[1][0] = a - lp;
        P[1][1] = a - lm;

        RR y0[] = {1., 1.};
        int nn = Run2(y0, n, &y[0], 1e-6);
        n = nn;
        y.resize((n + 1) * 3);
        // [Wi,We] =  P^-1 [ [ exp( lp t) ,0] [0, exp(lm t)]] P
    }

    int findinterval(RR t) {//
        int i = 0, j = n, k;
        if (t < y[3 * i]) return 0;//
        else if (t > y[j * 3])
            return -1;
        while (i < j - 1)
            if (y[3 * (k = (i + j) / 2)] > t) j = k;
            else i = k;
        return i;
    }

    RR fU(R tt, int i) {
        RR t = tt;
        int j = findinterval(t);
        if (j == -1) return uc;
        RR *y0 = &y[3 * j], *y1 = &y[3 * (j + 1)], dy0[2], dy1[2];
        RR dt = (*y1 - *y0);
        RR theta = (t - *y0) / dt;
        if (0) {
            RR cb = theta;
            RR ca = RR(1.) - cb;
            RR f = ca * y0[1 + i] + cb * y1[1 + i];
            return f;
        }
        F(y0 + 1, dy0);
        F(y1 + 1, dy1);
        RR cb = theta * theta * (RR(3.) - RR(2.) * theta);
        RR ca = RR(1.) - cb;
        RR cta = ((RR(1) - theta) * (RR(1) - theta) * theta) * dt;
        RR ctb = ((theta - RR(1)) * theta * theta) * dt;
        RR f = ca * y0[1 + i] + cb * y1[1 + i] + cta * dy0[i] + ctb * dy1[i];
        return f;
    }

    void F(RR *y, RR *F) {
        /**
        Unobserved probabilities ODEs:
        c1 = la + mu; c0 = psi (1 - p)
        **/
        F[0] = -mu * (y[0] - y[1]);
        F[1] = -c1 * y[1] + la * y[0] * y[1] + c0;
    }

    void D2F(RR *y, RR *dF, RR *ddF, RR *dddF = 0) {
        /**
        Observed probabilities ODEs:
        c1 = la + mu; c0 = psi (1 - p)
        **/
        F(y, dF);
        ddF[0] = -mu * (dF[0] - dF[1]);
        ddF[1] = -c1 * dF[1] + la * dF[0] * y[1] + la * y[0] * dF[1];
        if (dddF) {
            dddF[0] = -mu * (ddF[0] - ddF[1]);
            dddF[1] = -c1 * ddF[1] + la * ddF[0] * y[1] + la * y[0] * ddF[1] +
                      2 * la * dF[0] * dF[1];
        }
    }


    int F(RR *y, R h, RR *b) {// resolution y-h F(y) = b
        /*
         //     y[0] - h mu ( -y[0] + y[1] ) = b0
         //     -b0 + y[0] ( 1+ h mu) =  h mu  y[1]
         //     y[1] = -b/(mu h) + y[0] ( ( 1+ h mu)/ (mu h)
         //  put in 2 equations
         y1= a0 + a1 y0
         */
        RR hmu = h * mu, a0 = -b[0] / hmu, a1 = (RR(1.) + hmu) / hmu;
        //  y1 = a0 + a1 y0
        //  y1- h* (-c1*y[1] + la*y[0]*y[1] +c0) - b1 =0 ;

        //  eq aa y0^2 + bb y0 + cc
        // h = 0 => y = b
        RR aa = -h * la * a1, aa2 = 2 * aa;//  term en y0^2
        RR bb = a1 - h * (-c1 * a1 + la * a0);// term y0
        RR cc = a0 - h * (-c1 * a0 + c0) - b[1]; // term ctn
        RR delta = bb * bb - 4 * aa * cc;
        int nF = 0, k = 0;
        if (delta < 0) return 0;
        RR delta2 = sqrt(delta);

        RR y0p = (-bb + delta2) / aa2;
        RR y0m = (-bb - delta2) / aa2;
        // chosir la bonne racine !!!
        RR y1p = a0 + a1 * y0p;
        RR y1m = a0 + a1 * y0m;
        if (val(y0p) >= 0 && val(y1p) >= 0 && val(y0p) <= 1. && val(y1p) <= 1.) {
            y[k++] = y0p;
            y[k++] = y1p;
            ++nF;
        }
        if (val(y0m) >= 0 && val(y1m) >= 0 && val(y0m) <= 1. && val(y1m) <= 1.) {
            y[k++] = y0m;
            y[k++] = y1m;
            ++nF;
        }
        // verif
#ifndef NDEBUG
        if (nF > 0) {
            // assert(nF>0);
            RR f[2], e[2];
            F(y, f);
            e[0] = y[0] - h * f[0] - b[0];
            e[1] = y[1] - h * f[1] - b[1];
            if (abs(e[0]) + abs(e[1]) > 1e-1) {
                cout << " bizarre:   " << abs(e[0]) + abs(e[1]) << " nF " << " " << delta << " " << aa << " " << bb
                     << " " << cc << endl;
            }
        }
#endif
        return nF;
    }

    void J(RR *y, RR *J) {
        J[0] = mu;
        J[1] = -mu;
        J[2] = 2 * la * y[1];
        J[3] = c1 + la * y[0];
    }

    int Run0(RR *y0, int nx, RR *y, R eps = 1e-4) { // Euler explicite en : temps exp(-t*cc)
        // size of y = 3*(n+1)
        // exp(min(lm,lp)*T) = eps
        R T = log(eps) / min(val(lm), val(lp));

        if (debug >= debugVal) cout << " err " << exp(min(lm, lp) * T) << endl;
        RR f[2], df[2];
        RR *yp = y;
        R t = 0, tp = 0;//lt=exp(-t*cc);

        *yp++ = t;
        *yp++ = y0[0];
        *yp++ = y0[1];
        if (debug >= debugVal) cout << t << " " << y0[0] << " " << y0[1] << endl;
        //int i=0;
        for (int i = 0; i < nx; ++i) {
            tp = t;
            RR *yp = y + 3 * i;
            RR *yy = y + 3 * (i + 1);
            //y = yp + h *F(yp)

            D2F(yp + 1, f, df);
            // calcul de dt
            R n2 = max(abs(val(df[0])), abs(val(df[1])));
            // err = n2 * dt*dt/8
            R dt = sqrt(8 * eps / n2);
            //dt = min(dt,T-tp);
            t += dt;
            yy[0] = t;
            yy[1] = yp[1] + dt * f[0];
            yy[2] = yp[2] + dt * f[1];
//            if (debug >= debugVal)
//                cout << t << " " << yy[1] << " " << yy[2] << " " << max(abs(yy[1] - uc), abs(yy[2] - uc)) << endl;
            if (max(abs(yy[1] - uc), abs(yy[2] - uc)) < eps) return i + 1;
        }
        return nx;
    }

    int Run2(RR *y0, int n, RR *y, R eps = 1e-6) { //  Crank Nicolson implicite  en : temps exp(-t*cc)
        //  (yn1 - yn)/dt = (F(yn1)+F(un))/2
        //  yn1 - dt/2 F(yn1) = yn + dt F(yn)/2
        // size of y = 3*(n+1)
        //R cc = max(val(-lm),val(-lp));
        RR f[2], Y[4], df[2], ddf[2], dddf[3];
        if (0) {
            // petit verif
            RR f[2], y[2] = {0.5, 0.4};
            F(y, f);
            //  y - h F = b
            R h = 0.1;
            RR b[] = {y[0] - h * f[0], y[1] - h * f[1]};
            int ns = F(Y, h, b);
//            if (debug >= debugVal) cout << ns << " : " << Y[0] << " " << Y[1] << ", " << Y[2] << " " << Y[3] << endl;
        }
        RR *yp = y;
        R t = 0;
        *yp++ = t;
        *yp++ = y0[0];
        *yp++ = y0[1];
//        if (debug >= debugVal) cout << t << " " << y0[0] << " " << y0[1] << endl;
        int ns = 0;
        for (int i = 0; i < n; ++i) {
            RR *yp = y + 3 * i;
            RR *yy = y + 3 * (i + 1);
            //  yn1 - dt/2 F(yn1) = yn + dt/2 F(yn)
            D2F(yp + 1, df, ddf, dddf);
            RR n3 = max(abs(dddf[0]), abs(dddf[1]));
            //  err = dt^3*n3/8;
            R dt = pow(val(8 * eps / n3), 1. / 3.);
            // int F(R *y,R h,R *b)
            R h = dt / 2;
            F(yp + 1, f);
            RR b[] = {yp[1] + f[0] * h, yp[2] + f[1] * h};
            ns = F(Y, h, b);
            //ns =0;
            if (ns == 0) { //  go to Euler explicite
                //F(yp+1,f);
                R n2 = max(abs(val(ddf[0])), abs(val(ddf[1])));
                R dt = -sqrt(8 * eps / n2);

                Y[0] = yp[1] + dt * f[0];
                Y[1] = yp[2] + dt * f[1];
            }
            t += dt;
            int ii = 0;
            yy[0] = t;
            yy[1] = Y[ii];
            yy[2] = Y[ii + 1];

//            if (debug >= debugVal) cout << t << " " << yy[1] << " " << yy[2] << endl;
            if (max(abs(yy[1] - uc), abs(yy[2] - uc)) < eps)
                return i + 1;
        }
        return n;
    }

    /*
     Branche ODE sur une branche de temps l_b  ici Ue et Ui sont connue,
     t in [ t_b  , t_b + l_b]  // retrograde  (t decroit)
     Pe_t = mu*Pe(t) -mu *Pi(t)
     Pi_t = (c1 Pi(t) - la (Pi(t) Ue(T-t) + Pe(i) Ui(T-t)
     CI: AU TEMPS FINAL t_b-l_b
     Pe_t(t_b+l_b) = 0
     Pi_t(t_b+l_b) = (si term) psi p
     ( sinon) la*(Pi^l(t_b+l_b)Pe^r(t_b)+ Pi^r(t_b)Pe^r(t_b);
     */
    int ODE(R T, R delta, RR *P, R eps, int id = 0, const string &shift = emptystr) {
        return ODE(T, delta, P, 1, eps, id, shift);
    }

    int ODE(R T, R delta, RR *P, int np, R eps, int id = 0, const string &shift = emptystr) {// attention dt < 0 ...
        // P[0] == Pe, P[1] == Pi
        // et P donnee initiale en tb,
        R t = 0;// temps final

        int n = 0; // nb iteration
        // R t0=T;
        while (1) {
            n++;
            RR Pt[2], Ptt[2];
            RR U[] = {fU(T + t, 0), fU(T + t, 1)}, dU[2], ddU[2], dddU[3];
            D2F(U, dU, ddU, dddU);// derive de U pour calcul le pas de temps
            R n2 = 0;
            for (int k = 0; k < np; ++k) {
                int i0 = 2 * k, i1 = i0 + 1;
                Pt[0] = mu * P[i0] - mu * P[i1];
                Pt[1] = c1 * P[i1] - la * (P[i1] * U[0] + P[i0] * U[1]);
                Ptt[0] = mu * Pt[0] - mu * Pt[1];
                Ptt[1] = c1 * Pt[1]
                         - la * (Pt[1] * U[0] + Pt[0] * U[1])
                         - la * (P[1] * dU[0] + P[0] * dU[1]);
                n2 = max(n2, max(abs(val(Ptt[0])), abs(val(Ptt[1]))));
            }
            // err = dt^2*n2/8
            // Euler implicite pour commencer.
            R dt = sqrt((8. * eps) / n2);
            R trestant = delta - t;
            bool fini = false;
            if (dt * 1.2 > trestant) {
                dt = trestant;
                fini = true;
            }
            t += dt;
            //  **** Attention a signe de dt>0 ****
            //  (P^n1 - P^n)/  = dt A P^n1,
            //  M = Id - dt A
            // A =  [[    mu ,-mu ], [-la*Ui , c1-la*Ue ]]

            RR Ut[] = {fU(T + t, 0), fU(T + t, 1)};
            RR M1[2][2], M[2][2] = {
                    {RR(1) + mu * dt, -mu * dt},
                    {(-la * Ut[1]) * dt, RR(1.) + (c1 - la * Ut[0]) * dt}
            };
            RR det = inv(M, M1);
            for (int k = 0; k < np; ++k) {
                int k2 = k + k, i0 = k2, i1 = i0 + 1;
                RR b[2] = {P[i0], P[i1]};
//                if (debug >= debugVal)
//                    cout << shift << "  " << T + t << " " << P[i0] << " " << P[i1]
//                         << " " << Ut[0] << " " << Ut[1] << " " << id << " " << dt << "  TT " << endl;
                Mult(M1, b, P + k2);
            }
            if (fini) break;
        }
//        if (debug >= debugVal)
//            cout << shift << " ODEe " << id << " [ " << T << " .. " << T + t << "] , " << P[0] << " " << P[1] << " np="
//                 << np << endl;

        return n;
    }

private: // no copy
    BDEI_pb(const BDEI_pb &);

    void operator=(const BDEI_pb &);
};

struct TreeBranch {
    double t;
    double value;
    double T;
    TreeBranch *b[2];// 2 pointeur on 2 branch if existe
    int id;// id[1] == 0 => just a branch //
    TreeBranch() : t(0), value(0), T(0), b{0, 0}, id{-1} {}

    bool internal() const {
        assert((b[0] == 0 && b[1] == 0) == (id != -1));
        return b[0];
    }
};


Solution::Solution(R l,
    R mu_, R mu_min_, R mu_max_,
    R la_, R la_min_, R la_max_,
    R psi_, R psi_min_, R psi_max_,
    R p_, R p_min_, R p_max_,
    R cpu_time_, int nb_iter_)
        : likelihood(l), mu(mu_), la(la_), psi(psi_), p(p_),
        mu_max(mu_max_), mu_min(mu_min_),
        la_max(la_max_), la_min(la_min_),
        psi_max(psi_max_), psi_min(psi_min_),
        p_max(p_max_), p_min(p_min_),
        cpu_time(cpu_time_), nb_iter(nb_iter_) {}

Solution::Solution(R l, R mu_, R la_, R psi_, R p_, R cpu_time_, int nb_iter_)
        : likelihood(l), mu(mu_), la(la_), psi(psi_), p(p_),
        la_max(la_), la_min(la_),
        mu_max(mu_), mu_min(mu_),
        psi_max(psi_), psi_min(psi_),
        p_max(p_), p_min(p_),
        cpu_time(cpu_time_), nb_iter(nb_iter_) {}

struct Forest {
    vector<TreeBranch *> f;
    int nt, ni;// nb of tips ( leaves of tree), nb of internal node
    R T;

    Forest(string fn);

    long size() const { return f.size(); }

    TreeBranch *operator[](int i) const { return f[i]; }
};

R SetTime(TreeBranch *b, R t, int &nt, int &ni) {
    /* Sets the time at the beginning of the branch.
    The time is tree-specific, i.e. the time at the start of the root branch is 0 */

    if (!b) return t;
    b->t = t;
    R T = t + b->value;// T branchement
    R Tb = T;
    if (!b->b[0])
    {
        ++nt;// nb of tips (leaves)
    }
    else
    {
        ++ni;
        for (int i = 0; i < 2; ++i)
        {
            T = max(T, SetTime(b->b[i], Tb, nt, ni));
        }
    }
    return T;
}


TreeBranch *Read(ifstream &f, int i) {
    TreeBranch *p = 0;
    int c = f.peek();
    if (c == EOF) return 0;
    p = new TreeBranch;
    if (c == '(') {
        char cass[] = ",)";
        c = f.get();
        for (int step = 0; step < 2; ++step) {
            p->b[step] = Read(f, 2 * i + step);
            c = f.get();
            if (c != cass[step]) {
                throw std::invalid_argument("Invalid newick format of the input tree(s)");
            }
        }
        c = f.get();
        while ((c != ':') && (c != ';') && (c != EOF)) {
            // read internal node's name and ignore it
            c = f.get();
        }
        if (c == EOF) {
            throw std::invalid_argument("Invalid newick format of the input tree(s)");
        }
    } else {
        c = f.get();
        while ((c != ':') && (c != ';') && (c != EOF)) {
            // ignore the id in the newick file
            c = f.get();
        }
        if (c == EOF) {
            throw std::invalid_argument("Invalid newick format of the input tree(s)");
        }
        // set tip id
        p->id = i;
    }

    if (c == ':') {
        // read branch length into p->value
        f >> p->value;
        c = f.peek();
        // we expect the time T for this branch to be provided after the branch length in square brackets
        if (c == '[') {
            c = f.get();
            f >> p->T;
            c = f.get();
            if (c != ']') {
                throw std::invalid_argument("Invalid newick format of the input tree(s)");
            }
        }
    } else {
        if (c != ';') {
            throw std::invalid_argument("Invalid newick format of the input tree(s)");
        }
        // if the root branch length is not specified, set it to zero
        p->value = 0;
    }
    return p;
}

Forest::Forest(string fn)
        : T(0), nt(0), ni(0) {
    ifstream ff(fn);
    assert(ff);
    TreeBranch *p = 0;
    do {
        p = Read(ff, 1);

        if (p) {
            f.push_back(p);
            SetTime(p, 0, nt, ni);
            if (p->value) {
                // have not yet read the ;
                if (ff.get() != ';') {
                    throw std::invalid_argument("Invalid newick format of the input tree(s)");
                }
            }
            // skip any whitespaces, newlines etc. till the next tree start or the EOF
            int c = ff.peek();
            while((c != '(') && (c != EOF) && (c != ':')) {
                ff.get();
                c = ff.peek();
            }
        }
    } while (p);
    if (debug >= infoVal)
    {
        cout << "Observed forest contains "  << f.size() << " tree(s) with "
        << nt << " tips and " << ni << " internal nodes" << endl;
    }
}

typedef struct {
    double a, b;
} my_constraint_data;

template<class RR>
void ODE2(RR *v, BDEI_pb<RR> &pb, R t0, R t1, R eps = 1e-6) {
    RR p[] = {1., 0, 0., 1.};
    pb.ODE(t0, t1 - t0, p, 2, eps);
    v[0] = p[0];
    v[2] = p[1];
    v[1] = p[2];
    v[3] = p[3];
}

struct DataOde {
    R T, l1, l2;
    int iT, iT1, iT2; // for optimisation

    DataOde(R TT, R ll1, R ll2 = 0) : T(TT), l1(ll1), l2(ll2),
                                      iT(std::numeric_limits<int>::max()),
                                      iT1(std::numeric_limits<int>::max()),
                                      iT2(std::numeric_limits<int>::max()) { assert(T >= 0); }

    bool internal() const { return l2 > 0.; }

    template<class RR>
    RR operator()(BDEI_pb<RR> &pb, R eps = 1e-6) const {
        if (l2 == 0) { // root
            RR p[] = {0., 1.};
            pb.ODE(T - l1, l1, p, eps);
            return log(pb.piE * p[0] + pb.piI * p[1]);
        } else { // internal node
            RR pp[] = {0., 1., 0., 1.};
            pb.ODE(T - l1, l1, pp + 0, eps);
            pb.ODE(T - l2, l2, pp + 2, eps);
            return log(pp[0] * pp[3] + pp[1] * pp[2]);
        }
    }
};

inline ostream &operator<<(ostream &cout, const DataOde &d) {
    if (d.l2) cout << " i " << d.T << " " << d.l1 << " " << d.l2;
    else cout << " r " << d.T << " " << d.l1;
    if (d.iT1 >= 0)
        cout << " :  " << d.iT << " " << d.iT1 << " " << d.iT2;
    else if (d.iT >= 0)
        cout << " : " << d.iT << " " << d.iT1;
    return cout;
}

void SetDataOde(TreeBranch *b, R T, vector<DataOde> &vdo, int lvl) {
    if (lvl == 0)
    { //racine ...
        vdo.push_back(DataOde(T, b->value));
    }
    if (b->internal()) {
        R t = T - (b->t + b->value); // de
        vdo.push_back(DataOde(t, b->b[0]->value, b->b[1]->value));
        SetDataOde(b->b[0], T, vdo, lvl + 1);
        SetDataOde(b->b[1], T, vdo, lvl + 1);
    }
}


void SetDataOde(Forest &f, vector<DataOde> &vdo) {
    for (int i = 0; i < f.size(); ++i) {
        SetDataOde(f[i], f[i]->T, vdo, 0);
    }
}

int oneeee = 0;

double duratinit = 0, duratopt = 0, duratbuild = 0, duraterr = 0;


template<class RR>
RR JCout(RR *x, R pie, R u, R ut, const vector<DataOde> &vdo, const vector<tuple<double, int, int >> &vs, R eps = 1e-5);

class J_vdo {
public:
    const vector<DataOde> &vdo;
    vector<tuple<double, int, int >> *vs;
    const R *p0;
    R u;
    R ut;
    R eps;
    R pie;
    long count;
    int debug;
    vector<int> num;
    int nn;

    double Jm;
    vector<double> xopt;
    double Jopt;//  value
    int dd, dsens;
    vector<double> dir; //
    int nnewton;

    J_vdo(const vector<DataOde> &vdoo, const R *pp, R piee, R uu, R uut, R epss = 1e-6, int dd = 0)
            : vdo(vdoo), vs(0), p0(pp), pie(piee), u(uu), ut(uut), eps(epss), count(0), debug(dd), num(4), Jm(nan("")), xopt(), dd(0),
              dsens(0), dir(), nnewton(0) { init(); }

    J_vdo(const vector<DataOde> &vdoo, vector<tuple<double, int, int >> &vss, const R *pp, R piee, R uu, R uut, R epss = 1e-6,
          int dd = 0)
            : vdo(vdoo), vs(&vss), p0(pp), pie(piee), u(uu), ut(uut), eps(epss), count(0), debug(dd), num(4), Jm(nan("")), xopt(), dd(0),
              dsens(0), dir(4), nnewton(0) { init(); }

    void set(int ddd, int dds) {
        dd = ddd;
        dsens = dds;
        assert(abs(dds) == 1);
        assert(dd >= 0 && dd < nn);
    }

    void set(vector<double> &dd) { // direction ...
        dir = dd;
        assert(dd.size() == 4);
    }

    J_vdo(const vector<DataOde> &vdoo, vector<tuple<double, int, int >> &vss, const R *pp, R piee, R uu, R uut, R epss, int dd,
          double JJm, vector<double> xxopt, double JJopt)
            : vdo(vdoo), vs(&vss), p0(pp), pie(piee), u(uu), ut(uut), eps(epss), count(0),
              debug(dd), num(4), Jm(JJm), xopt(xxopt), Jopt(JJopt),
              dd(-1), dsens(-1), nnewton(0) { init(); }

    void init() {
        nn = 0;
        for (int i = 0; i < 4; ++i)
            if (p0[i] < 0) num[nn++] = i;
        num.resize(nn);
        vector<R> pp(4);
        copy(p0, p0 + 4, pp.begin());
        nnewton = 0;
        if (debug >= debugVal) cout << "J_vdo:  num = " << num << " :: " << pp << endl;
    }

    R operator()(const std::vector<double> &xx, std::vector<double> &grad) {
        vector<double> x(4);
        copy(p0, p0 + 4, x.begin());
        for (int i = 0; i < nn; ++i)
            x[num[i]] = xx[i];
        return J(x, grad);
    }

    R J(std::vector<double> &x, std::vector<double> &grad) {
        count++;
        R Cout = 0;
        assert(x.size() == 4);
        for (int i = 0; i < 4; ++i)
            x[i] = max(1e-6, x[i]);
        // sampling probability
        x[3] = min(x[3], 0.99999999);
        assert(x.size() == 4);

        int n = (int) grad.size();
        // CIs
        if (n == 1) { // grad in direction dir ..
            assert(dir.size() == 4); // dir is def ..
            typedef Diff<R> DR;
            vector<DR> dx(4);
            for (int i = 0; i < 4; ++i)
                dx[i] = DR(x[i], dir[i]);
            DR dCout = JCout(&dx[0], pie, u, ut, vdo, *vs, eps);
            grad[0] = dCout.dval;
            Cout = dCout.val;
        } // optimum search
        else {
            assert(n == nn || n == 4);
            typedef Diff<R> DR;
            vector<DR> dx(4);
            R ccout = 0;
            for (int j = 0; j < n; ++j) {
                if (n == 4) // all derivative
                    for (int i = 0; i < 4; ++i)
                        dx[i] = DR(x[i], R(i == j));
                else
                    for (int i = 0; i < 4; ++i)
                        dx[i] = DR(x[i], R(i == num[j]));
                DR Cout = JCout(&dx[0], pie, u, ut, vdo, *vs, eps);
                grad[j] = Cout.dval;
                if (ccout == 0) ccout = Cout.val;
            }
            Cout = ccout;
        }
        if (debug >= debugVal) {
            cout << " J " << x << "  = " << Cout << " " << count << " / " << nn << endl;
        }
        return Cout;
    }

    void DJ(const std::vector<double> &xx, std::vector<double> &grad) {
        R h = 1e-3;

        for (int i = 0; i < grad.size(); ++i) {
            std::vector<double> x = xx, dx;
            x[i] += h;
            R cp = operator()(x, dx);
            x[i] -= 2 * h;
            R cm = operator()(x, dx);
            grad[i] = (cp - cm) / (2 * h);
        }
    }

    double J1(const std::vector<double> &xx, std::vector<double> &grad) {
        double J = (xx[dd] - xopt[num[dd]]) * dsens;
        if (!grad.empty()) {
            std::fill(grad.begin(), grad.end(), 0.);
            grad[dd] = dsens;
        }
        return J;
    }

    void rhominmax(R &rhomin, R &rhomax) const {
        rhomin = 0;// rho is positif ...
        rhomax = 1000;
        for (int i = 0; i < 4; ++i) {
            R xmin = 1e-6, di = dir[i], xi = xopt[i];
            //  xi +di*rho > xmin  so:
            // di>0 : rho > (xmin-xi)/di
            // di<0 : rho < (xmin-xi)/di

            if (di > 1e-6) rhomin = max(rhomin, (xmin - xi) / di);
            else if (di < -1e-6) rhomax = min(rhomax, (xmin - xi) / di);
            if (i == 3) {
                R xmax = 0.9999999;
                //  xi +di*rho < xmax  so:
                // d>0 : rho < (xmax-xi)/di
                // d<0 : rho > (xmax-xi)/d

                if (di > 1e-6) rhomax = min(rhomax, (xmax - xi) / di);
                else if (di < -1e-6) rhomin = max(rhomin, (xmax - xi) / di);
            }
        }

    }

    double J2(R rho, R &dg) { //
        vector<double> x(xopt), g1(1);
        for (int i = 0; i < 4; ++i)
            x[i] += rho * dir[i];

        double jj = val(J(x, g1));
        dg = g1[0];
        return jj;
    }

    double J2(R rho, vector<double> &dg) { //
        vector<double> x(xopt);
        for (int i = 0; i < 4; ++i)
            x[i] += rho * dir[i];
        return J(x, dg);
    }

    double Newton(double rho0, int niter = 100) {
        const int ndj = 10;
        double ddd[ndj];
        assert(rho0 > 0);
        double rmin, rmax; //  borne on rho ..
        rhominmax(rmin, rmax);
        double rdd = min((rmax - rmin) * 0.01, 0.01);
        double rmy = (min(rmax, 1.) + rmin) * 0.5;
        rho0 = max(rmin + rdd, rho0);
        rho0 = min(rmax - rdd, rho0);
        double lim = 0.1;
        for (int step = 0; step < 3; ++step) {

            assert(rho0);
            // double sig = copysign(1.0, rho0);
            double eps = 1e-1;
            double rho = rho0, jr, djr, rhop = rho0;
            double krho = 1;

            for (int i = 0; i < niter; ++i) {
                nnewton++;
                rhop = rho;
                jr = J2(rho, djr) - Jm;
                if (abs(djr) < 1e-4) break;
                if (abs(jr) < eps)
                    return rho;

                double dd = max(min(lim, jr / djr), -lim);//  limiteur
                // fonction croissante en rho ???
                // jr < 0 => dd <0
                // jr > 0 => dd > 0

                if (jr < 0 && dd > 0) dd = -0.1;///
                if (jr > 0 && dd < 0) dd = 0.1;///
                ddd[i % ndj] = dd;
                if (abs(dd) < 1e-4 && abs(jr) > 100 && rho > 1e-2) {
                    rho /= 3.;
                    continue;
                }

                rho -= dd;
                //rho = sig*max(rho*sig,0.0001);//  pour etre sur de ne pas changer de signe ...
                rho = max(rmin, min(rmax, rho));
                if (abs(rho - rhop) < 1e-6) {
                    krho += 1.;
                    if (rho < rmy) rho = rho0 * krho;// bof Bof F. Kecht ..
                    else rho = rho0 / krho;
                    rho = max(rmin, min(rmax, rho));
                }
            }
            if (debug >= debugVal)
                cout << "Pb::  Newton 1d  do not converge !!!!!!" << rho0 << " dir = "
                     << dir[0] << " " << dir[1] << " " << dir[2] << " " << dir[3] << " step = " << step << " retry "
                     << endl;
            {
                rho0 /= 2;
                lim /= 2.;
            }
        }
        if (debug >= debugVal)
            cout << " Pb:: Newton Err: 1d do converge (return Nan) data:"
                 << dir[0] << " " << dir[1] << " " << dir[2] << " " << dir[3] << endl;
        return nan("");
    }

    double C1(const std::vector<double> &xx, std::vector<double> &grad) {
        double C = operator()(xx, grad) - Jm;
        return C;
    }

    static double wrap(const std::vector<double> &x, std::vector<double> &grad, void *data) {
        return (*reinterpret_cast<J_vdo *>(data))(x, grad);
    }

    static double wrapJ(const std::vector<double> &x, std::vector<double> &grad, void *data) {
        return (reinterpret_cast<J_vdo *>(data))->J1(x, grad);
    }

    static double wrapC(const std::vector<double> &x, std::vector<double> &grad, void *data) {
        return (reinterpret_cast<J_vdo *>(data))->C1(x, grad);
    }

};

template<class RR>
vector<RR> Build(const vector<tuple<double, int, int >> &vs, BDEI_pb<RR> &pb, R eps = 1e-6) {
    int ivend = (int) vs.size();
    int iend = get<2>(vs[ivend - 1]);
    vector<RR> m22(iend * 4);
    R t0, t1 = get<0>(vs[0]);
    int k0 = 0, k1 = get<2>(vs[0]), kk = 0;
    {
        thread_pool tp(size_pool);

        for (int i = 0; i < ivend; ++i) {
            k0 = k1;
            t0 = t1;
            t1 = get<0>(vs[i]);
            k1 = get<2>(vs[i]);
            if (k0 < k1)// true interval
            {
                assert(k1 - 1 == kk / 4);//  ok ..
                RR *pm = &m22[kk];
                tp.enqueue_work([pm, &pb, t0, t1, eps]() {
                    ODE2(pm, pb, t0, t1, eps);
                });
                kk += 4;
            }
        }
    }
    return m22;
}

template<class RR>
void mult(RR *m, RR *p) {
    RR x[] = {p[0], p[1]};
    p[0] = m[0] * x[0] + m[1] * x[1];
    p[1] = m[2] * x[0] + m[3] * x[1];
}

template<class RR>
RR JCout(RR *x, R pie, R u, R ut, const vector<DataOde> &vdo, const vector<tuple<double, int, int >> &vs, R eps) {
    RR Cout = 0;
    auto start = high_resolution_clock::now();
    BDEI_pb<RR> pb(x, pie);
    auto afterpb = high_resolution_clock::now();
    int ni = 0;// number of internal node ..
    int n = (int) vdo.size();
    vector<RR> m22 = Build(vs, pb, eps);
    auto afterbuild = high_resolution_clock::now();
    vector<RR> pp(4 * vdo.size());
    RR *pm22 = &m22[0];
    RR *ppp = &pp[0];
    RR pi[] = {pb.piE, pb.piI};
    RR *ppi = &pi[0];
    {
        thread_pool tp(size_pool);

        // calculate the likelihood contribution of every node
        for (int i = 0; i < vdo.size(); ++i) {
            const DataOde &vdoi = vdo[i];
            int kp = i * 4; //  init
            RR *p4 = ppp + kp;
            if (vdoi.internal()) ni++;
            tp.enqueue_work([vdoi, p4, pm22, ppi]() {

                p4[0] = 0; // PEI (left)
                p4[1] = 1; // PII (left)
                p4[2] = 0; // PEI (right)
                p4[3] = 1; // PII (right)

                int iT = vdoi.iT;
                for (int j = vdoi.iT1; j < iT; ++j) {
                    mult(pm22 + 4 * j, p4);
                }
                for (int j = vdoi.iT2; j < iT; ++j)// IT2 => max int  if unset => no loop
                    mult(pm22 + 4 * j, p4 + 2);

                if (vdoi.internal())
                    p4[0] = log(p4[0] * p4[3] + p4[1] * p4[2]);  // internal node
                else
                    p4[0] = log(ppi[0] * p4[0] + ppi[1] * p4[1]); // root branch
            });
        }

    }
    // calcul du cout ..
    Cout = 0;
    for (int i = 0; i < n; ++i)
        Cout += pp[i * 4];

    int fs = n - ni;
    int nt = ni + fs;//  nb de feuilles/ tips
    Cout += nt * log(pb.psi * pb.p); // sampling of tips
    Cout += (nt - fs) * log(pb.la); // transmissions

    if (u != R(0.))
    {
        RR hidden_prob = pb.piE * pb.fU(ut, 0) + pb.piI * pb.fU(ut, 1);
        // calculate u if needed
        if (u < R(0.))
        {
            RR u_val = R(fs) / (RR(1.) - hidden_prob) * hidden_prob;
            Cout += u_val * log(hidden_prob);
//            cout << "Adding " << u_val << "hidden trees" << endl;
        }
        else
        {
            Cout += u * log(hidden_prob);
        }
    }
    auto end = high_resolution_clock::now();

    duratinit += duration_cast<duration<double>>(afterpb - start).count();
    duratbuild += duration_cast<duration<double>>(afterbuild - afterpb).count();
    duratopt += duration_cast<duration<double>>(end - start).count();

    return -Cout;
}


Solution *ErrRandDir(J_vdo &jvdoC, int nnum, int *num, int ndir, double *derr, int debug, R cpu_time, int nb_iter) {
    auto start = high_resolution_clock::now();

    jvdoC.nnewton = 0;
    normal_distribution<double> loi(0, 1.);
    vector<double> d(4, 0.);
    double r[4];
    fill(derr, derr + 8, 0.);
    int n2 = nnum * 2;
    int m2 = 1 << nnum; // vertex of hypercube ...
    int knan = 0, kdir = 0;

    for (int i = 0; i < ndir + n2 + m2; ++i) {

        if (i < n2) // face of hypercude
        {
            int k = i / 2;
            int s = (i % 2) * 2 - 1;// +- 1
            for (int j = 0; j < nnum; ++j)
                r[j] = (k == j) * s;
        } else if (i < n2 + m2)// vertex of hypercube
        {
            int k = i - n2;
            for (int j = 0; j < nnum; ++j) {
                int b = k & (1 << j);
                r[j] = b ? 1 : -1;
            }
        } else {
            for (int j = 0; j < nnum; ++j)
                r[j] = loi(rand_gen);
        }

        for (int j = 0; j < nnum; ++j)
            d[num[j]] = r[j];
        double ld = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2] + d[3] * d[3]);
        for (int j = 0; j < 4; ++j)
            d[j] /= ld;

        jvdoC.set(d);
        //  double rhom = jvdoC.Newton(-0.05);

        double rhop = jvdoC.Newton(0.01);
        if (isfinite(rhop))
            // to skip NaN
        {
            for (int j = 0; j < 4; ++j) {
                int j2 = 2 * j;
                double ej = d[j] * rhop;
                derr[j2] = min(derr[j2], ej);
                derr[j2 + 1] = max(derr[j2 + 1], ej);
            }
        } else knan++;
        kdir++;
    }
    auto end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end - start);

    return new Solution(-jvdoC.Jopt,
                        jvdoC.xopt[0], jvdoC.xopt[0] + derr[2 * 0], jvdoC.xopt[0] + derr[2 * 0 + 1],
                        jvdoC.xopt[1], jvdoC.xopt[1] + derr[2 * 1], jvdoC.xopt[1] + derr[2 * 1 + 1],
                        jvdoC.xopt[2], jvdoC.xopt[2] + derr[2 * 2], jvdoC.xopt[2] + derr[2 * 2 + 1],
                        jvdoC.xopt[3], jvdoC.xopt[3] + derr[2 * 3], jvdoC.xopt[3] + derr[2 * 3 + 1],
                        cpu_time + duration.count() / 1000.f, nb_iter + jvdoC.count);
}

vector<DataOde> &getForestDataODE(Forest &forest, vector<DataOde> &vdo, int _debug) {
    SetDataOde(forest, vdo);
    return vdo;
}


vector<tuple<double, int, int>> &getVS(vector<DataOde> &vdo, vector<tuple<double, int, int>> &vs) {

    for (int i = 0; i < vdo.size(); ++i) {
        DataOde &vdoi = vdo[i];
        vs.push_back(make_tuple(vdoi.T, i * 4 + vdoi.internal(), 0));
        vs.push_back(make_tuple(vdoi.T - vdoi.l1, i * 4 + 2, 0));
        if (vdoi.internal())
            vs.push_back(make_tuple(vdoi.T - vdoi.l2, i * 4 + 3, 0));
    }
    sort(vs.begin(), vs.end());

    int nc = 0;
    R sst = 0, sso = 0, t1 = get<0>(vs[0]), t0;
    int nk = 0;
    for (int i = 0; i < vs.size(); ++i) {//
        t0 = t1;
        t1 = get<0>(vs[i]);
        if (t1 - t0 > 1e-4) nk++;

        int from = get<1>(vs[i]);
        int ivdo = from / 4;
        int cas = from % 4;
        if (cas == 0) nc--;
        else if (cas == 1) nc -= 2;
        else nc++;
        if (cas < 2) vdo[ivdo].iT = nk; // END ..
        else if (cas == 2) vdo[ivdo].iT1 = nk;// BEGIN
        else if (cas == 3) vdo[ivdo].iT2 = nk;
        sst += (t1 - t0) * nc;
        sso += (t1 - t0);
        get<2>(vs[i]) = nk;
    }
    assert(nc == 0);
    return vs;
}


R calcLikelihood(Forest &forest, R mu, R lambda, R psi, R p, R pie, R u, R ut, int nt, int _debug) {
    if (nt < 1) {
        size_pool = thread::hardware_concurrency();
    } else {
        size_pool = nt;
    }

    vector<DataOde> vdo;
    vdo = getForestDataODE(forest, vdo, _debug);

    vector<tuple<double, int, int >> vs;
    vs = getVS(vdo, vs);

    R predefdata[] = {mu, lambda, psi, p};
    double eps = 1e-6;
    J_vdo jvdo(vdo, vs, predefdata, pie, u, ut, eps, _debug >= debugVal);

    nlopt::opt opt(nlopt::LD_MMA, 1);
    vector<double> lb(1, 1e-5), ub(1, HUGE_VAL), xj(1);
    ub[0] = mu + 1e-5;
    lb[0] = mu - 1e-5;
    xj[0] = mu;

    double minf;
    opt.set_lower_bounds(lb);
    opt.set_upper_bounds(ub);
    opt.set_initial_step(0.1);
    opt.set_min_objective(J_vdo::wrap, &jvdo);
    opt.set_xtol_rel(1e-5);
    nlopt::result result = opt.optimize(xj, minf);
    return -minf;
}



R calculateLikelihood(const string &treename, R mu, R lambda, R psi, R p, R pie, R u, R ut, int nt, int _debug) {
    debug = _debug;
    Forest forest(treename);
    return calcLikelihood(forest, mu, lambda, psi, p, pie, u, ut, nt, _debug);
}


Solution
*inferParameters(const string &treename, R *x0, const R *dub, R pie, R mu, R lambda, R psi, R p, R u, R ut,
                int nbdirerr, int nt, int debug_, int nStarts) {
    debug = debug_;

    Solution *s = nullptr;
    int num[4] = {0, 1, 2, 3}; // mu, la, psi , p

    double eps = 1e-6;
    int algo = 1, old = 0;
    if (nt < 1) {
        size_pool = thread::hardware_concurrency();
    } else {
        size_pool = nt;
    }
    nlopt::algorithm thealgo[] = {nlopt::LN_NEWUOA_BOUND, nlopt::LD_MMA, nlopt::LD_LBFGS, nlopt::GN_DIRECT_L};
    random_device reng;
    rand_gen.seed(reng());
    algo = algo % 4;
    algo = min(3, max(0, algo));


    Forest forest(treename);
    vector<DataOde> vdo;
    vdo = getForestDataODE(forest, vdo, debug_);

    vector<tuple<double, int, int >> vs;
    vs = getVS(vdo, vs);

    R cpu_time = 0.00;
    int nb_iter = 0;

    if (debug >= debugVal)
        cout << "forest file: " << treename << endl;
//             << ", mu  = " << mu << ", lambda = " << lambda << ", psi = " << psi << ", p = " << p << endl;
    R predefdata[] = {mu, lambda, psi, p};
    vector<double> xsol(4); // to store the solution ...
    int nbdata = 0;

    for (int k = 0; k < 4; ++k)
        if (predefdata[k] < 0) {
            num[nbdata++] = k;
        } else {
            if (debug >= infoVal) cout << " the value of " << BDEIParams[k] << " is fixed to " << predefdata[k] << endl;
        }
    assert(nbdata > 0);
    array<R, 4> sol;
    R Jsol;

    auto start = high_resolution_clock::now();

    J_vdo jvdo(vdo, vs, predefdata, pie, u, ut, eps, debug >= debugVal);
    nlopt::opt opt(thealgo[1], nbdata);
    vector<double> lb(nbdata, 1e-5), ub(nbdata, HUGE_VAL);
    for (int i = 0; i < nbdata; ++i)
        ub[i] = min(dub[num[i]], ub[i]);

    opt.set_lower_bounds(lb);
    opt.set_upper_bounds(ub);


    double bestRes =  std::numeric_limits<double>::infinity();
    vector<double> x(4); // copy min ...
    double minf;
    // if there are several starting points repeat the procedure for all of them and pick the best result
    for (int startI = 0; startI < nStarts; startI++) {
        opt.set_initial_step(0.1);
        opt.set_min_objective(J_vdo::wrap, &jvdo);
        opt.set_xtol_rel(1e-5);
        vector<double> xj(nbdata);
        for (int i = 0; i < nbdata; ++i) {
            xj[i] = x0[startI * 4 + num[i]];
        }
        if (debug >= infoVal)
        {
            cout << "Optimisation attempt " << startI + 1 << ":" << endl;
            cout << "  starting values for optimised parameters: " << xj << endl;
            cout << "  upper bounds : " << ub << endl;
        }
        bool ok = false;
        try {
            nlopt::result result = opt.optimize(xj, minf);
            ok = true;
            if (debug >= infoVal)
            {
                cout << "  optimised parameter values: " << xj << endl;
                cout << "  optimised log-likelihood: " << setprecision(10) << -minf << endl;
            }
        }
        catch (exception &e) {
            if (debug >= errorVal)  cout << "  nlopt failed: " << e.what() << endl;
        }
        nb_iter += jvdo.count;

        if ((minf < bestRes) && ok) {
            bestRes = minf;
            copy(predefdata, predefdata + 4, x.begin());
            for (int i = 0; i < nbdata; ++i)
                x[num[i]] = xj[i];
            xsol = x;
            copy(x.begin(), x.end(), x0);
        } else if (ok) {
            break;
        }
    }
    auto end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end - start);
    cpu_time = duration.count() / 1000.f;

    if (debug >= debugVal)
    {
        cout << " CPUtime  = " << cpu_time << " s" << endl;
    }

    for (int j = 0; j < 4; ++j)
        sol[j] = x[j];
    Jsol = minf;

    if (nbdirerr > 0) {
        if (debug >= infoVal)  cout << "CI estimation" << endl;
        double Jerr = khi_2_95[1];// err in on direction => 1 data  not nbdata
        double Js = Jsol + Jerr;
        vector<double> x(4), dir(4);
        copy(sol.begin(), sol.end(), x.begin());
        J_vdo jvdoC(vdo, vs, predefdata, pie, u, ut, eps, debug > 1, Js, x, Jsol);
        vector<double> errb(8, 0.);
        s = ErrRandDir(jvdoC, nbdata, num, nbdirerr, &errb[0], debug, cpu_time, nb_iter);
    } else {
        s = new Solution(-Jsol, sol[0], sol[1], sol[2], sol[3], cpu_time, nb_iter);
    }
    return s;
}


Solution
*inferParameters(const string &treename, R *x0, const R *dub, R pie, R mu, R lambda, R psi, R p, R u, R ut,
                int nbdirerr, int nt, int nstarts) {
    return inferParameters(treename, x0, dub, pie, mu, lambda, psi, p, u, ut, nbdirerr, nt, infoVal, nstarts);
}

Solution
*inferParameters(const string &treename, R *x0, const R *dub, R pie, R mu, R lambda, R psi, R p, R u, R ut, int nbdirerr, int nt) {
    return inferParameters(treename, x0, dub, pie, mu, lambda, psi, p, u, ut, nbdirerr, nt, infoVal, 1);
}
