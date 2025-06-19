#pragma once

#include "stdio.h"
#include "math.h"
#include "geometry.c"

#define TWOPI 6.283185307
#define SEC_PER_DAY 86400.0
#define EPOCH_JAN1_12H_2000 2451545.0 // Jan  1.5 2000 = Jan  1 2000 12h UTC
#define OMEGA_E 1.00273790934 // earth rotation per sideral day


// read polar coefficients xp and yp on date MJD from "earth orientation parameters" file
int get_earth_positions(int MJD, double *xp, double *yp) {
    FILE *f = fopen("eop.txt", "r");
    char buf[104];
    int line_MJD;
    while (fgets(buf, 104, f)) {
        sscanf(buf + 11, "%d", &line_MJD);
        if (line_MJD == MJD) {
            sscanf(buf + 17, "%lf%lf", xp, yp);
            fclose(f);
            return 0;
        }
    }
    fclose(f);
    return -1;
}


int get_delta_UT1_UTC(int MJD, double *UT1_UTC) {
    FILE *f = fopen("eop.txt", "r");
    char buf[104];
    int line_MJD;
    while (fgets(buf, 104, f)) {
        sscanf(buf + 11, "%d", &line_MJD);
        if (line_MJD == MJD) {
            sscanf(buf + 17, "%lf", UT1_UTC);
            fclose(f);
            return 0;
        }
    }
    fclose(f);
    return -1;
}


int get_LOD(int MJD, double *LOD) {
    FILE *f = fopen("eop.txt", "r");
    char buf[104];
    int line_MJD;
    while (fgets(buf, 104, f)) {
        sscanf(buf + 11, "%d", &line_MJD);
        if (line_MJD == MJD) {
            sscanf(buf + 48, "%lf", LOD);
            fclose(f);
            return 0;
        }
    }
    fclose(f);
    return -1;
}


int get_delta_TAI_UTC(int MJD, double *TAI_UTC) {
    FILE *f = fopen("eop.txt", "r");
    char buf[104];
    int line_MJD;
    while (fgets(buf, 104, f)) {
        sscanf(buf + 11, "%d", &line_MJD);
        if (line_MJD == MJD) {
            sscanf(buf + 99, "%lf", TAI_UTC);
            fclose(f);
            return 0;
        }
    }
    fclose(f);
    return -1;
}


Matrix calculatePolarMotionMatrix(double julianDate) {
    int MJD = (int)(julianDate - 2400000.5);
    double xp, yp;
    get_earth_positions(MJD, &xp, &yp);
    xp *= 4.84813681e-6; // from arcseconds to radians
    yp *= 4.84813681e-6;
    double sin_xp = sin(xp);
    double sin_yp = sin(yp);
    double cos_xp = cos(xp);
    double cos_yp = cos(yp);
    Matrix polarMotionMatrix = {cos_xp,          0,        -sin_xp,
                                sin_xp * sin_yp, cos_yp,   cos_xp * sin_yp,
                                sin_xp * cos_yp, -sin_yp,  cos_xp * cos_yp};
    return polarMotionMatrix;
}


double JD(int year, double day) {
    year--;

    // Centuries are not leap years unless they divide by 400
    int A = (year / 100);
    int B = 2 - A + (A / 4);
 
    double jan01 = trunc(365.25 * year) +
                   trunc(30.6001 * 14)  +
                   1720994.5 + B;  // 1720994.5 = Oct 30, year -1

    return jan01 + day;
}


double MJD(int year, double day) {
    return JD(year, day) - 2400000.5;
}


double FromJan1_12h_2000(double julian_date) {
    return julian_date - EPOCH_JAN1_12H_2000;
}


double TLE_to_JD(double tle_time) {
    int year = (int)(tle_time * 1e-3);
    double day = tle_time - year * 1e3;
    return JD(year, day);
}


double JD_to_GMST(double julian_date) {
    double UT = fmod(julian_date + 0.5, 1.0);
    double TU = (FromJan1_12h_2000(julian_date) - UT) / 36525.0;
 
    double GMST = 24110.54841 + TU *
                  (8640184.812866 + TU * (0.093104 - TU * 6.2e-06));
 
    GMST = fmod(GMST + SEC_PER_DAY * OMEGA_E * UT, SEC_PER_DAY);
 
    if (GMST < 0.0) {
       GMST += SEC_PER_DAY;  // "wrap" negative modulo value
    }
 
    return  (TWOPI * (GMST / SEC_PER_DAY));
}


double JD_to_LMST(double julian_date, double lon) {
    return fmod(JD_to_GMST(julian_date) + lon, TWOPI);
}


int TEME_to_ECEF(Vector position_TEME, Vector velocity_TEME, double julian_date, Vector *position_ECEF, Vector *velocity_ECEF) {
    int MJD_UTC = (int)(julian_date - 2400000.5);

    double UT1_UTC;
    get_delta_UT1_UTC(MJD_UTC, &UT1_UTC);
    double gmst = JD_to_GMST(julian_date + UT1_UTC / SEC_PER_DAY);

    double TAI_UTC;
    get_delta_TAI_UTC(MJD_UTC, &TAI_UTC);
    double JD_TT = julian_date + (TAI_UTC + 32.184) / SEC_PER_DAY;
    double T_TT = FromJan1_12h_2000(JD_TT) / 36525.0;

    double meanLongitudeAscendingNodeMoon = (125.04452222 + T_TT * (-(5*360 + 134.1362608) + T_TT * (0.0020708 + T_TT * 2.2e-6))) * TWOPI / 360.0;
    meanLongitudeAscendingNodeMoon = fmod(meanLongitudeAscendingNodeMoon, TWOPI);
    if(julian_date > 2450449.5) {
        gmst += (0.00264 * sin(meanLongitudeAscendingNodeMoon) + 0.000063 * sin(2 * meanLongitudeAscendingNodeMoon)) * 0.5 * TWOPI / 648000;
    }

    Matrix PEF_TOD = {cos(gmst), -sin(gmst), 0,
                      sin(gmst), cos(gmst),  0,
                      0,         0,          1};

    Vector position_PEF = mul_mat_vec(PEF_TOD, position_TEME);
    Matrix polarMotionMatrix = calculatePolarMotionMatrix(julian_date);
    *position_ECEF = mul_mat_vec(polarMotionMatrix, position_PEF);

    double LOD;
    get_LOD(MJD_UTC, &LOD);
    double omegaEarth = 7.29211514670698e-05 * (1.0  - LOD / 86400.0);

    Vector velocity_PEF = mul_mat_vec(PEF_TOD, velocity_TEME);
    velocity_PEF.x += omegaEarth * position_PEF.y;
    velocity_PEF.y -= omegaEarth * position_PEF.x;
    *velocity_ECEF = mul_mat_vec(polarMotionMatrix, velocity_PEF);
}