#pragma once

#define Q0 120.0             // parameter for the SGP4 density function
#define S 1.01222928         // parameter for the SGP4 density function
#define Q0MS2T 1.88027916e-9 // (Q0 - S) ** 4
#define XJ2 1.082616e-3      // the second gravitational zonal harmonic of the Earth
#define XJ3 -0.253881e-5     // the third gravitational zonal harmonic of the Earth
#define XJ4 -1.65597e-6      // the fourth gravitational zonal harmonic of the Earth
#define XKE 0.743669161e-1   // sqrt(GM) where G is gravitational constant and M is the mass of the Earth
#define XKMPER 6378.135      // kilometers/Earth radii
#define XMNPDA 1440.0        // time units/day
#define AE 1.0               // distance units/Earth radii

#define DE2RA 0.174532925e-1 // degree to radian

#define PIO2 1.57079633
#define PI 3.14159265
#define X3PIO2 4.71238898
//#define TWOPI 6.2831853

#define E6A 1.0e-6
#define TOTHRD 0.6666666667

// #define CK2 5.413080E-4      // 0.5 * XJ2 * AE ** 2
// #define CK4 0.62098875E-6    // -0.375 * XJ4 * AE ** 4
const double CK2 = 0.5 * XJ2 * AE * AE;
const double CK4 = -0.375 * XJ4 * AE * AE * AE * AE;
const double A3OVK2 = -XJ3 * AE * AE * AE / CK2;

#define MU 398600.4416e9  // gravitational parameter of the Earth (GM)




typedef struct sgp4_init_params_t{
    double calc_time,
        M_dot, omega_dot, node_dot,
        eta,
        del_M0,
        sin_M0,
        sin_i0, cos_i0,
        a0_dp, n0_dp,
        x3thm1, x1mth2, x7thm1,
        C1, C2, C3, C4, C5, D2, D3, D4,
        node_coef, omega_coef, M_coef, L_coef, ay_coef,
        T2COF, T3COF, T4COF, T5COF;
    int simp_flag;
} sgp4_init_params_t;


typedef struct tle_set_t{
    double M,     // mean anomaly, deg
           node,  // ascending node latitude, deg
           omega, // argument of periapsis, deg
           e,     // eccentricity
           i,     // inclination, deg
           n,     // mean motion, revolutions per day
           B_star,
           epoch;
} tle_set_t;


typedef struct kepler_set_t{
    double M,     // mean anomaly, deg
           node,  // ascending node latitude, deg
           omega, // argument of periapsis, deg
           e,     // eccentricity
           i,     // inclination, deg
           n;     // mean motion, revolutions per day
} kepler_set_t;
