#include "sgp4.c"

int main() {
    // Vector pos = {2328.970670, -5995.220833, 1719.970699}, vel = {2.912072, -0.983415, -7.090817};
    // kepler_set_t osc_elems;
    // osculating_elements(pos, vel, &osc_elems);
    // printf("%lf %lf %lf %lf %lf %lf\n", osc_elems.M, osc_elems.node, osc_elems.omega, osc_elems.e, osc_elems.i, osc_elems.n);
    

    // double **M = malloc(6 * sizeof(double *)), **N = malloc(6 * sizeof(double *));
    // int i, j;
    // for (i = 0; i < 6; i++) {
    //     M[i] = malloc(6 * sizeof(double));
    //     N[i] = malloc(6 * sizeof(double));
    //     for (j = 0; j < 6; j++) {
    //         scanf("%lf", M[i] + j);
    //     }
    // }
    // inv_matrix(6, M, N);
    // for (i = 0; i < 6; i++) {
    //     for (j = 0; j < 6; j++) {
    //         printf("%lf ", N[i][j]);
    //     }
    //     printf("\n");
    // }
    // for (i = 0; i < 5; i++) {
    //     free(M[i]);
    //     free(N[i]);
    // }
    // free(M);
    // free(N);


    // double tle_arr[6] = {110.5714, 115.9689, 52.6988, .0086731, 72.8435, 16.05824518}, osc_arr[6];
    // tle_to_oscul(tle_arr, 0.66816e-4, osc_arr);
    // for (int i = 0; i < 6; i++) {
    //     printf("%lf ", osc_arr[i]);
    // }
    // printf("\n");

    tle_set_t tle;
    Vector pos = {2328.970670, -5995.220833, 1719.970699}, vel = {2.912072, -0.983415, -7.090817};
    calculate_tle(pos, vel, 24275.98708465, 0.66816e-4, &tle);

    tle.node *= DE2RA;
    tle.omega *= DE2RA;
    tle.M *= DE2RA;
    tle.i *= DE2RA;
    tle.n *= TWOPI / XMNPDA;

    sgp4_init_params_t params;
    SGP4_init(&tle, &params);
    SGP4_continue(&tle, 0, params, &pos, &vel);

    pos = mul_vec_num(pos, XKMPER / AE);
    vel = mul_vec_num(vel, XKMPER / AE * XMNPDA / 86400.0);

    printf("%lf %lf %lf %lf %lf %lf\n", pos.x, pos.y, pos.z, vel.x, vel.y, vel.z);

    return 0;
}