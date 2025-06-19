#include "stdio.h"
#include "math.h"

#include "sgp4.c"


// #define ONTHRD 0.3333333333
// #define ONSEVN 0.1428571429
// #define ONENIN 0.1111111111
// #define ONTHRFAC 0.1666666667
// #define ONFOUFAC 0.4166666667
// #define ONFIVFAC 8.3333333333e-3
// #define ONSIXFAC 1.3888888889e-3
// #define ONSEVFAC 1.984126984e-4
// #define ONEIGFAC 2.48015873e-5
// #define ONNINFAC 2.755731922e-6


// double sinm(double x) {
//     double x2 = x * x;
//     return x * (1 + x2 * (-ONTHRFAC + x2 * (ONFIVFAC + x2 * (-ONSEVFAC + x2 * ONNINFAC))));
// }


// double cosm(double x) {
//     double x2 = x * x;
//     return 1 + x2 * (-ONFOUFAC + x2 * (ONSIXFAC + x2 * (-ONEIGFAC)));
// }


// double atanm(double x) {
//     double x2 = x * x;
//     return x * (1 + x2 * (-ONTHRD + x2 * (0.2 + x2 * (-ONSEVN + x2 * ONENIN))));
// }


// float FastSqrt(float x) {
//     int i = *(int *)&x;
//     i = 0x1fbd1df5 + (i >> 1);
//     float y = *(float *)&i;
//     //y *= 0.5f * (1.0f + x / (y * y));
//     return y;
// }


// return 1 if str begins with "exit"
int beg_exit(char *str) {
    return (str[0] == 'e') && (str[1] == 'x') && (str[2] == 'i') && (str[3] == 't');
}


// read input parameters for calculation and print result; return -1 if was exited, else 0
int trajectory_calc() {
    char buf[20];
    double t_start, t_finish, step;
    int frame;

    printf("insert range\n");
    scanf("%s", buf);
    if (beg_exit(buf)) {
        return -1;
    }
    sscanf(buf, "%lf", &t_start);
    scanf("%s", buf);
    if (beg_exit(buf)) {
        return -1;
    }
    sscanf(buf, "%lf", &t_finish);

    printf("insert step\n");
    scanf("%s", buf);
    if (beg_exit(buf)) {
        return -1;
    }
    sscanf(buf, "%lf", &step);
    if ((t_finish - t_start) * step < 0) { // уточнить этот момент
        step = -step;
    }

    printf("select frame\n");
    scanf("%s", buf);
    
    while (!beg_exit(buf)) {
        sscanf(buf, "%d", &frame);
        switch (frame) {
        case 0: // TEME
            SGP4_driver(t_start, t_finish, step, 0);
            return 0;
            break;
        case 1: // ECEF
            SGP4_driver(t_start, t_finish, step, 1);
            return 0;
            break;
        default:
            printf("unknown symbol\n");
            break;
        }
        printf("select frame\n");
        scanf("%s", buf);
    }
    return -1;
}


int direction_to_target() {

}


int direction_from_station() {

}


int input_TLE() {
    int i;
    char line1[71], line2[71];
    printf("input first line\n");
    for (i = 0; i < 70; i++) {
        scanf("%c", line1 + i);
    }
    printf("input second line\n");
    for (i = 0; i < 70; i++) {
        scanf("%c", line2 + i);
    }

    FILE *f = fopen("tle.txt", "w");
    for (i = 1; i < 70; i++) {
        fprintf(f, "%c", line1[i]);
    }
    fprintf(f, "\n");
    for (i = 1; i < 70; i++) {
        fprintf(f, "%c", line2[i]);
    }
    fclose(f);
    return 0;
}


// calculate TLE from current position and velocity and store it to memory
int update_TLE() {
    char buf[20];
    double x, y, z, vx, vy, vz;
    int frame;

    printf("insert coordinates\n");
    scanf("%s", buf);
    if (beg_exit(buf)) {
        return -1;
    }
    sscanf(buf, "%lf", &x);
    
    scanf("%s", buf);
    if (beg_exit(buf)) {
        return -1;
    }
    sscanf(buf, "%lf", &y);

    scanf("%s", buf);
    if (beg_exit(buf)) {
        return -1;
    }
    sscanf(buf, "%lf", &z);

    scanf("%s", buf);
    if (beg_exit(buf)) {
        return -1;
    }
    sscanf(buf, "%lf", &vx);

    scanf("%s", buf);
    if (beg_exit(buf)) {
        return -1;
    }
    sscanf(buf, "%lf", &vy);

    scanf("%s", buf);
    if (beg_exit(buf)) {
        return -1;
    }
    sscanf(buf, "%lf", &vz);


    printf("select frame\n");
    scanf("%s", buf);
    
    while (!beg_exit(buf)) {
        sscanf(buf, "%d", &frame);
        switch (frame) {
        case 0: // TEME
            
            return 0;
            break;
        case 1: // ECEF
            
            return 0;
            break;
        default:
            printf("unknown symbol\n");
            break;
        }
        printf("select frame\n");
        scanf("%s", buf);
    }
    return -1;
}


int main() {
    char buf[20];
    int mode;
    printf("select mode\n");
    scanf("%s", buf);
    
    while (!beg_exit(buf)) {
        sscanf(buf, "%d", &mode);
        switch (mode) {
        case 0:
            trajectory_calc();
            break;
        case 1:
            direction_to_target();
            break;
        case 2:
            direction_from_station();
            break;
        case 3:
            input_TLE();
            break;
        case 4:
            // update_TLE();
            Vector pos = {2328.970670, -5995.220833, 1719.970699}, vel = {2.912072 -0.983415 -7.090817};
            calculate_tle(pos, vel, 24275.98708465, 0.66816e-4);
            break;
        default:
            printf("unknown symbol\n");
            break;
        }

        printf("select mode\n");
        scanf("%s", buf);
    }

    return 0;
}

// 24275.98708465 24276.98708465

