#pragma once

#include "stdio.h"
#include "math.h"
#include "stdlib.h"

#include "common.h"

#include "geometry.c"
#include "frames.c"
#include "tle.c"


double mod_2pi(double x) {
    int i = (int)(x / TWOPI);
    x -= i * TWOPI;
    if (x < 0) x += TWOPI;
    return x;
}


double arctg2(double sinx, double cosx) {
    if (sinx == 0) return 0;
    if (cosx == 0) {
        if (sinx > 0) return PIO2;
        return X3PIO2;
    }
    if (cosx > 0) {
        if (sinx > 0) {
            return atan(sinx / cosx);
        }
        return TWOPI + atan(sinx / cosx);
    }
    return PI + atan(sinx / cosx);
}


int read_params(sgp4_init_params_t *params) {
    FILE *f = fopen("init.txt", "r");
    if (!f) {
        return -1;
    }
    fscanf(f, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%d",
        &params->calc_time, &params->M_dot, &params->omega_dot, &params->node_dot, &params->eta, &params->del_M0, &params->sin_M0,
        &params->sin_i0, &params->cos_i0, &params->a0_dp, &params->n0_dp, &params->x3thm1, &params->x1mth2, &params->x7thm1,
        &params->C1, &params->C2, &params->C3, &params->C4, &params->C5, &params->D2, &params->D3, &params->D4,
        &params->node_coef, &params->omega_coef, &params->M_coef, &params->L_coef, &params->ay_coef, &params->T2COF,
        &params->T3COF, &params->T4COF, &params->T5COF, &params->simp_flag);
    fclose(f);
}


void write_params(sgp4_init_params_t *params) {
    FILE *f = fopen("init.txt", "w");
    fprintf(f, "%.20lf\n%.20lf\n%.20lf\n%.20lf\n%.20lf\n%.20lf\n%.20lf\n%.20lf\n%.20lf\n%.20lf\n%.20lf\n%.20lf\n%.20lf\n%.20lf\n%.20lf\n%.20lf\n%.20lf\n%.20lf\n%.20lf\n%.20lf\n%.20lf\n%.20lf\n%.20lf\n%.20lf\n%.20lf\n%.20lf\n%.20lf\n%.20lf\n%.20lf\n%.20lf\n%.20lf\n%d",
        params->calc_time, params->M_dot, params->omega_dot, params->node_dot, params->eta, params->del_M0, params->sin_M0,
        params->sin_i0, params->cos_i0, params->a0_dp, params->n0_dp, params->x3thm1, params->x1mth2, params->x7thm1,
        params->C1, params->C2, params->C3, params->C4, params->C5, params->D2, params->D3, params->D4,
        params->node_coef, params->omega_coef, params->M_coef, params->L_coef, params->ay_coef, params->T2COF,
        params->T3COF, params->T4COF, params->T5COF, params->simp_flag);
    fclose(f);
}


void SGP4_init(tle_set_t *tle, sgp4_init_params_t *params) {

    params->del_M0 = 0;
    params->L_coef = params->ay_coef = 1;
    params->T2COF = params->T3COF = params->T4COF = params->T5COF = 0;

    double a1 = pow((XKE / tle->n), TOTHRD);
    params->cos_i0 = cos(tle->i);
    params->sin_i0 = sin(tle->i);
    double theta_sq = params->cos_i0 * params->cos_i0;
    params->x3thm1 = 3. * theta_sq - 1.;
    double e0_sq = tle->e * tle->e;
    double beta0_sq = 1. - e0_sq;
    double beta0 = sqrt(beta0_sq);
    double delta1 = 1.5 * CK2 * params->x3thm1 / (a1 * a1 * beta0 * beta0_sq);
    double a0 = a1 * (1. - delta1 * (.5 * TOTHRD + delta1 * (1. + 134. / 81. * delta1)));
    double delta0 = 1.5 * CK2 * params->x3thm1 / (a0 * a0 * beta0 * beta0_sq);
    params->n0_dp = tle->n / (1. + delta0);
    params->a0_dp = a0 / (1. - delta0);

    // INITIALIZATION

    params->simp_flag = ((params->a0_dp * (1. - tle->e) / AE) < (220. / XKMPER + AE));

    double S4 = S;
    double QOMS24 = Q0MS2T;
    double perigee = (params->a0_dp * (1. - tle->e) - AE) * XKMPER;

    if (perigee < 156.) {
        S4 = perigee - 78.;
        if (perigee <= 98)
        {
            S4 = 20.;
        }
        double temp = ((Q0 - S4) * AE / XKMPER);
        QOMS24 = temp * temp * temp * temp;
        S4 = S4 / XKMPER + AE;
    }

    double p_inv_sq = 1. / (params->a0_dp * params->a0_dp * beta0_sq * beta0_sq);
    double xi = 1. / (params->a0_dp - S4);
    params->eta = params->a0_dp * tle->e * xi;
    double eta_sq = params->eta * params->eta;
    double e_eta = tle->e * params->eta;
    double psi_sq = fabs(1. - eta_sq);
    double coef = QOMS24 * xi * xi * xi * xi;
    double coef1 = coef / sqrt(psi_sq * psi_sq * psi_sq * psi_sq * psi_sq * psi_sq * psi_sq);

    params->C2 = coef1 * params->n0_dp * (params->a0_dp * (1. + 1.5 * eta_sq + e_eta * (4. + eta_sq))
                    + .75 * CK2 * xi / psi_sq * params->x3thm1 * (8. + 3. * eta_sq * (8. + eta_sq)));
    params->C1 = tle->B_star * params->C2;
    params->C3 = coef * xi * A3OVK2 * params->n0_dp * AE * params->sin_i0 / tle->e;
    params->x1mth2 = 1. - theta_sq;
    params->C4 = 2. * params->n0_dp * coef1 * params->a0_dp * beta0_sq * (params->eta * (2. + .5 * eta_sq)
                    + tle->e * (.5 + 2. * eta_sq) - 2. * CK2 * xi / (params->a0_dp * psi_sq) * (-3. * params->x3thm1
                    * (1. - 2. * e_eta + eta_sq * (1.5 - .5 * e_eta)) + .75 * params->x1mth2 * (2. * eta_sq - e_eta
                    * (1. + eta_sq)) * cos(2. * tle->omega)));
    params->C5 = 2. * coef1 * params->a0_dp * beta0_sq * (1. + 2.75 * (eta_sq + e_eta) + e_eta * eta_sq);
    double theta_4 = theta_sq * theta_sq;

    double temp1 = 3. * CK2 * p_inv_sq * params->n0_dp;
    double temp2 = temp1 * CK2 * p_inv_sq;
    double temp3 = 1.25 * CK4 * p_inv_sq * p_inv_sq * params->n0_dp;
    params->M_dot = params->n0_dp + .5 * temp1 * beta0 * params->x3thm1 + .0625 * temp2 * beta0 * (13. - 78. * theta_sq + 137. * theta_4);
    double x1m5th = 1. - 5. * theta_sq;
    params->omega_dot = -.5 * temp1 * x1m5th + .0625 * temp2 * (7. - 114. * theta_sq + 395. * theta_4) + temp3 * (3. - 36. * theta_sq + 49. * theta_4);
    double XHDOT1 = -temp1 * params->cos_i0;
    params->node_dot = XHDOT1 + (.5 * temp2 * (4. - 19. * theta_sq) + 2. * temp3 * (3. - 7. * theta_sq)) * params->cos_i0;
    params->omega_coef = tle->B_star * params->C3 * cos(tle->omega);
    params->M_coef = -TOTHRD * coef * tle->B_star * AE / e_eta;
    params->node_coef = 3.5 * beta0_sq * XHDOT1 * params->C1;
    params->T2COF = 1.5 * params->C1;
    params->L_coef = .125 * A3OVK2 * params->sin_i0 * (3. + 5. * params->cos_i0) / (1. + params->cos_i0);

    params->ay_coef = .25 * A3OVK2 * params->sin_i0;
    double x1pecM0 = 1. + params->eta * cos(tle->M);
    params->del_M0 = x1pecM0 * x1pecM0 * x1pecM0;
    params->sin_M0 = sin(tle->M);
    params->x7thm1 = 7. * theta_sq - 1.;

    if (!(params->simp_flag)) {
        double C1_sq = params->C1 * params->C1;
        params->D2 = 4. * params->a0_dp * xi * C1_sq;
        double temp = params->D2 * xi * params->C1 / 3.;
        params->D3 = (17. * params->a0_dp + S4) * temp;
        params->D4 = .5 * temp * params->a0_dp * xi * (221. * params->a0_dp + 31. * S4) * params->C1;
        params->T3COF = params->D2 + 2. * C1_sq;
        params->T4COF = .25 * (3. * params->D3 + params->C1 * (12. * params->D2 + 10. * C1_sq));
        params->T5COF = .2 * (3. * params->D4 + 12. * params->C1 * params->D3 + 6. * params->D2 * params->D2 + 15. * C1_sq * (2. * params->D2 + C1_sq));
    }
}


void SGP4_continue(tle_set_t *tle, double t_since, sgp4_init_params_t params, Vector *position, Vector *velocity) {        

    double M_DF = tle->M + params.M_dot * t_since;
    double omega_DF = tle->omega + params.omega_dot * t_since;
    double node_DF = tle->node + params.node_dot * t_since;
    double omega = omega_DF;
    double M_p = M_DF;
    double t_sq = t_since * t_since;
    double node = node_DF + params.node_coef * t_sq;
    double temp_a = 1. - params.C1 * t_since;
    double temp_e = tle->B_star * params.C4 * t_since;
    double temp_L = params.T2COF * t_sq;

    if (!params.simp_flag) {
        double del_omega = params.omega_coef * t_since;
        double x1pecMDF = 1. + params.eta * cos(M_DF);
        double del_M = params.M_coef * (x1pecMDF * x1pecMDF * x1pecMDF - params.del_M0);
        double temp = del_omega + del_M;
        M_p = M_DF + temp;
        omega = omega_DF - temp;
        double t_cube = t_sq * t_since;
        double t_four = t_since * t_cube;
        temp_a = temp_a - params.D2 * t_sq - params.D3 * t_cube - params.D4 * t_four;
        temp_e = temp_e + tle->B_star * params.C5 * (sin(M_p) - params.sin_M0);
        temp_L = temp_L + params.T3COF * t_cube + t_four * (params.T4COF + t_since * params.T5COF);
    }

    double a = params.a0_dp * temp_a * temp_a;
    double E = tle->e - temp_e;
    double L = M_p + omega + node + params.n0_dp * temp_L;
    double beta = sqrt(1. - E * E);
    double sqrt_a = sqrt(a);
    double n = XKE / (a * sqrt_a);

    // LONG PERIOD PERIODICS

    double ax_N = E * cos(omega);
    double temp = 1. / (a * beta * beta);
    double L_L = temp * params.L_coef * ax_N;
    double ay_NL = temp * params.ay_coef;
    double L_T = L + L_L;
    double ay_N = E * sin(omega) + ay_NL;

    // SOLVE KEPLERS EQUATION

    double U = mod_2pi(L_T - node);
    double Epw = U;
    double sin_Epw, cos_Epw, temp2, temp3, temp4, temp5, temp6;
    do {
        temp2 = Epw;
        sin_Epw = sin(temp2);
        cos_Epw = cos(temp2);
        temp3 = ax_N * sin_Epw;
        temp4 = ay_N * cos_Epw;
        temp5 = ax_N * cos_Epw;
        temp6 = ay_N * sin_Epw;
        Epw = (U - temp4 + temp3 - temp2) / (1. - temp5 - temp6) + temp2;
    } while (fabs(Epw - temp2) > E6A);

    double e_cosE = temp5 + temp6;
    double e_sinE = temp3 - temp4;
    double e_L_sq = ax_N * ax_N + ay_N * ay_N;
    temp = 1. - e_L_sq;
    double p_L = a * temp;
    double r = a * (1. - e_cosE);
    double temp1 = 1. / r;
    double r_dot = XKE * sqrt_a * e_sinE * temp1;
    double rf_dot = XKE * sqrt(p_L) * temp1;
    temp2 = a * temp1;
    double beta_L = sqrt(temp);
    temp3 = 1. / (1. + beta_L);
    double cos_u = temp2 * (cos_Epw - ax_N + ay_N * e_sinE * temp3);
    double sin_u = temp2 * (sin_Epw - ay_N - ax_N * e_sinE * temp3);
    double u = arctg2(sin_u, cos_u);
    double sin_2u = 2. * sin_u * cos_u;
    double cos_2u = 2. * cos_u * cos_u - 1.;
    temp = 1. / p_L;
    temp1 = CK2 * temp;
    temp2 = temp1 * temp;

    // SHORT PERIODICS

    double r_k = r * (1. - 1.5 * temp2 * beta_L * params.x3thm1) + .5 * temp1 * params.x1mth2 * cos_2u;
    double u_k = u - .25 * temp2 * params.x7thm1 * sin_2u;
    double node_k = node + 1.5 * temp2 * params.cos_i0 * sin_2u;
    double i_k = tle->i + 1.5 * temp2 * params.cos_i0 * params.sin_i0 * cos_2u;
    double r_dot_k = r_dot - n * temp1 * params.x1mth2 * sin_2u;
    double rf_dot_k = rf_dot + n * temp1 * (params.x1mth2 * cos_2u + 1.5 * params.x3thm1);

    // ORIENTATION VECTORS

    double sin_u_k = sin(u_k);
    double cos_u_k = cos(u_k);
    double sin_i_k = sin(i_k);
    double cos_i_k = cos(i_k);
    double sin_node_k = sin(node_k);
    double cos_node_k = cos(node_k);

    double Mx = -sin_node_k * cos_i_k;
    double My = cos_node_k * cos_i_k;

    Vector U_vec = {Mx * sin_u_k + cos_node_k * cos_u_k,
                My * sin_u_k + sin_node_k * cos_u_k,
                sin_i_k * sin_u_k};
    Vector V_vec = {Mx * cos_u_k - cos_node_k * sin_u_k,
                My * cos_u_k - sin_node_k * sin_u_k,
                sin_i_k * cos_u_k};

    // POSITION AND VELOCITY

    *position = mul_vec_num(U_vec, r_k);
    *velocity = add_vec(mul_vec_num(U_vec, r_dot_k), mul_vec_num(V_vec, rf_dot_k));
}


void SGP4_driver(double start_time, double finish_time, double step, int frame) {
    //double M0, node0, omega0, e0, i0, n0, B_star, epoch;
    double t;
    tle_set_t tle;
    Vector position, velocity;

    switch (get_elems(&tle)) {
    case 1:
        printf("invalid TLE\n");
        break;
    case 2:
        printf("wrong checksum\n");
        break;
    default:
        break;
    }

    tle.node *= DE2RA;
    tle.omega *= DE2RA;
    tle.M *= DE2RA;
    tle.i *= DE2RA;
    tle.n *= TWOPI / XMNPDA;

    sgp4_init_params_t params;
    if ((read_params(&params) == -1) || (params.calc_time < tle.epoch)) {
        SGP4_init(&tle, &params);
        params.calc_time = tle.epoch;
        write_params(&params);
    }

    t = start_time;
    while (t <= finish_time) {
        SGP4_continue(&tle, (TLE_to_JD(t) - TLE_to_JD(tle.epoch)) * XMNPDA, params, &position, &velocity);

        position = mul_vec_num(position, XKMPER / AE);
        velocity = mul_vec_num(velocity, XKMPER / AE * XMNPDA / 86400.0);
        
        if (frame == 1) {
            TEME_to_ECEF(position, velocity, TLE_to_JD(t), &position, &velocity);
        }

        printf("%lf %lf %lf %lf %lf %lf %lf\n", t, position.x, position.y, position.z, velocity.x, velocity.y, velocity.z);
        t += step;
    }
}


void osculating_elements(Vector position, Vector velocity, kepler_set_t *oscul_elems) {
    //printf("%lf %lf %lf\n", velocity.x, velocity.y, velocity.z);
    position = mul_vec_num(position, 1000);  // SI conversion
    velocity = mul_vec_num(velocity, 1000);

    Vector h = cross_product(position, velocity);
    //printf("%lf %lf %lf\n", h.x, h.y, h.z);
    Vector K = {0, 0, 1};
    Vector n_vec = cross_product(K, h);
    double n = vec_len(n_vec);
    //printf("%lf %lf %lf\n", n_vec.x, n_vec.y, n_vec.z);

    double v = vec_len(velocity);
    double r = vec_len(position);

    //printf("%lf %lf %lf\n", velocity.x, velocity.y, velocity.z);
    //printf("%lf %lf\n", r, v);

    Vector e_vec = add_vec(mul_vec_num(position, v * v / MU - 1. / r),
                           mul_vec_num(velocity, -dot_product(position, velocity) / MU));
    oscul_elems->e = vec_len(e_vec);

    double E = 0.5 * v * v - MU / r;
    double a = -0.5 * MU / E;
    //printf("%lf\n", a);
    oscul_elems->n = sqrt(MU / (a * a * a)) * 86400. / TWOPI;
    oscul_elems->i = acos(h.z / vec_len(h)) / DE2RA;

    double epsilon = 1e-10;

    if (fabs(oscul_elems->i) < epsilon) {
        oscul_elems->node = 0;
        //printf("a\n");
        if (fabs(oscul_elems->e) < epsilon) {
            oscul_elems->omega = 0;
        }
        else {
            oscul_elems->omega = acos(e_vec.x / oscul_elems->e) / DE2RA;
        }
    }
    else {
        oscul_elems->node = acos(n_vec.x / n) / DE2RA;
        //printf("%lf %lf\n", n_vec.x, n);
        //printf("b\n");
        //printf("%lf\n", oscul_elems->node);
        if (n_vec.y < 0) {
            oscul_elems->node = 180. - oscul_elems->node;
            //printf("c\n");
        }
        oscul_elems->omega = acos(dot_product(n_vec, e_vec) / (n * oscul_elems->e)) / DE2RA;
    }

    double f;
    if (fabs(oscul_elems->e) < epsilon) {
        if (fabs(oscul_elems->i) < epsilon) {
            f = acos(position.x / r);
            if (velocity.x > 0) {
                f = TWOPI - f;
            }
        }
        else {
            f = acos(dot_product(n_vec, position) / (n * r));
            if (dot_product(n_vec, velocity) > 0) {
                f = TWOPI - f;
            }
        }
    }
    else {
        if (e_vec.z < 0) {
            oscul_elems->omega = TWOPI - oscul_elems->omega;
        }
        f = acos(dot_product(e_vec, position) / (oscul_elems->e * r));
        if (dot_product(position, velocity) < 0) {
            f = TWOPI - f;
        }
    }

    double ecc_anom = 2 * atan(tan(0.5 * f) * sqrt((1 - oscul_elems->e) / (1 + oscul_elems->e)));
    oscul_elems->M = (ecc_anom - oscul_elems->e * sin(ecc_anom)) / DE2RA;
}


double det(int n, double **M) {
    if (n == 1) {
        return M[0][0];
    }
    if (n == 2) {
        return M[0][0] * M[1][1] - M[1][0] * M[0][1];
    }

    int i, j;
    double **M_1_j = malloc((n - 1) * sizeof(double *));
    for (i = 0; i < n - 1; i++) {
        M_1_j[i] = malloc((n - 1) * sizeof(double));
        for (j = 0; j < n - 1; j++) {
            M_1_j[i][j] = M[i + 1][j + 1];
        }
    }

    double d = 0;
    for (j = 0; j < n - 1; j++) {
        d += (j % 2 ? 1 : -1) * M[0][j] * det(n - 1, M_1_j);

        for (i = 0; i < n - 1; i++) {
            M_1_j[i][j] = M[i + 1][j];
        }
    }
    d += (j % 2 ? 1 : -1) * M[0][j] * det(n - 1, M_1_j);

    for (i = 0; i < n - 1; i++) {
        free(M_1_j[i]);
    }
    free(M_1_j);

    return d * (n % 2 ? 1 : -1);
}


int get_minor(int n, double **M, int a, int b, double **res) {
    int i, j;
    for (i = 0; i < a; i++) {
        for (j = 0; j < b; j++) {
            res[i][j] = M[i][j];
        }
        for (j = b + 1; j < n; j++) {
            res[i][j - 1] = M[i][j];
        }
    }
    for (i = a + 1; i < n; i++) {
        for (j = 0; j < b; j++) {
            res[i - 1][j] = M[i][j];
        }
        for (j = b + 1; j < n; j++) {
            res[i - 1][j - 1] = M[i][j];
        }
    }
    return 0;
}


int inv_matrix(int n, double **M, double **res) {
    double d = det(n, M);
    int i, j;
    double **M_i_j = malloc((n - 1) * sizeof(double *));
    for (i = 0; i < n - 1; i++) {
        M_i_j[i] = malloc((n - 1) * sizeof(double));
    }
    //printf("e\n");

    //int k, l;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            get_minor(n, M, i, j, M_i_j);
            //printf("h\n");
            res[j][i] = ((i + j) % 2 ? -1 : 1) * det(n - 1, M_i_j) / d;

            // printf("%d %d %lf\n", i, j, det(n - 1, M_i_j));
            // for (k = 0; k < n - 1; k++) {
            //     for (l = 0; l < n - 1; l++) {
            //         printf("%lf ", M_i_j[k][l]);
            //     }
            //     printf("\n");
            // }
        }
    }

    //printf("f\n");

    for (i = 0; i < n - 1; i++) {
        free(M_i_j[i]);
    }
    free(M_i_j);

    return 0;
}


void tle_to_oscul(double *tle_arr, double B_star, double *oscul_arr) {
    tle_set_t tle_elems = {tle_arr[0] * DE2RA, tle_arr[1] * DE2RA, tle_arr[2] * DE2RA, tle_arr[3],
                           tle_arr[4] * DE2RA, tle_arr[5] * TWOPI / XMNPDA, B_star, 0};
    //tle_set_t tle_elems = {tle_arr[0], tle_arr[1], tle_arr[2], tle_arr[3], tle_arr[4], tle_arr[5], B_star, 0};
    sgp4_init_params_t params;
    Vector pos, vel;

    SGP4_init(&tle_elems, &params);
    SGP4_continue(&tle_elems, 0, params, &pos, &vel);
    pos = mul_vec_num(pos, XKMPER / AE);
    vel = mul_vec_num(vel, XKMPER / AE * XMNPDA / 86400.0);

    // printf("%lf %lf %lf\n", pos.x, pos.y, pos.z);
    // printf("%lf %lf %lf\n", vel.x, vel.y, vel.z);

    kepler_set_t oscul_elems;
    osculating_elements(pos, vel, &oscul_elems);

    oscul_arr[0] = oscul_elems.M;
    oscul_arr[1] = oscul_elems.node;
    oscul_arr[2] = oscul_elems.omega;
    oscul_arr[3] = oscul_elems.e;
    oscul_arr[4] = oscul_elems.i;
    oscul_arr[5] = oscul_elems.n;

    return;
}


void calculate_tle(Vector pos, Vector vel, double new_epoch, double B_star, tle_set_t *res) {
    
    double **M_mat = malloc(6 * sizeof(double *)), **M_inv = malloc(6 * sizeof(double *)),
           sum, tle_arr[6], osc_arr[6], osc_arr_0[6], osc_arr_true[6], delta_tle[6], delta_osc[6];
    int i, j;

    for (i = 0; i < 6; i++) {
        M_mat[i] = malloc(6 * sizeof(double));
        M_inv[i] = malloc(6 * sizeof(double));
    }

    kepler_set_t oscul;
    osculating_elements(pos, vel, &oscul);
    osc_arr_true[0] = oscul.M;
    osc_arr_true[1] = oscul.node;
    osc_arr_true[2] = oscul.omega;
    osc_arr_true[3] = oscul.e;
    osc_arr_true[4] = oscul.i;
    osc_arr_true[5] = oscul.n;

    for (i = 0; i < 6; i++) {
        tle_arr[i] = osc_arr_true[i];
    }

    double epsilon = 1e-16;
    Vector prev_pos = pos, prev_vel = vel;

    res->M = tle_arr[0] * DE2RA;
    res->node = tle_arr[1] * DE2RA;
    res->omega = tle_arr[2] * DE2RA;
    res->e = tle_arr[3];
    res->i = tle_arr[4] * DE2RA;
    res->n = tle_arr[5] * TWOPI / XMNPDA;

    int k = 0;
    double err = 1.;
    while ((fabs(err) > epsilon) && (k++ < 1000)) {
        tle_to_oscul(tle_arr, B_star, osc_arr_0);

        for (i = 0; i < 6; i++) {
            delta_tle[i] = 0.001 * tle_arr[i];
            tle_arr[i] += delta_tle[i];
            tle_to_oscul(tle_arr, B_star, osc_arr);
            for (j = 0; j < 6; j++) {
                delta_osc[j] = osc_arr[j] - osc_arr_0[j];
                M_mat[j][i] = delta_osc[j] / delta_tle[i];
            }
            tle_arr[i] -= delta_tle[i];
        }

        inv_matrix(6, M_mat, M_inv);

        for (i = 0; i < 6; i++) {
            sum = 0;
            for (j = 0; j < 6; j++) {
                sum += M_inv[i][j] * (osc_arr_0[j] - osc_arr_true[j]);
            }
            tle_arr[i] -= sum;
        }

        err = 0;
        for (i = 0; i < 6; i++) {
            err += osc_arr_0[i] - osc_arr_true[i];
        }
    }

    res->M = tle_arr[0];
    res->node = tle_arr[1];
    res->omega = tle_arr[2];
    res->e = tle_arr[3];
    res->i = tle_arr[4];
    res->n = tle_arr[5];
    res->B_star = B_star;
    res->epoch = new_epoch;

    for (i = 0; i < 6; i++) {
        free(M_mat[i]);
        free(M_inv[i]);
    }
    free(M_mat);
    free(M_inv);

}
