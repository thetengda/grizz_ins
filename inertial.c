//
// Created by tdhuang on 2022/4/12.
//
#include "inertial.h"
#include "inertial_math.h"

extern size_t INS_decode_binary(const char *ins, void **data) {

    FILE *fp = fopen(ins, "rb");

    if (!fp) {
        printf("file not exist   %s", ins);
        return 0;
    }
    fseek(fp, 0, SEEK_END);
    size_t len = ftell(fp);
    rewind(fp);
    *data = malloc(len);
    if (!(*data)) {
        printf("INS file too big, malloc failed");
        return 0;
    }
    if (!fread(*data, len, 1, fp)) {
        printf("read wrong, can not get data");
    }
    fclose(fp);
    return len;
}

extern int INS_update_param(INS_data_t *data, INS_result_t *state, INS_result_t *last) {
    int i;

    state->sin_lat = sin(state->lat);
    state->sin_lon = sin(state->lon);
    state->cos_lat = cos(state->lat);
    state->cos_lon = cos(state->lon);

    state->sin_roll = sin(state->roll);
    state->cos_roll = cos(state->roll);
    state->sin_pitch = sin(state->pitch);
    state->cos_pitch = cos(state->pitch);
    state->sin_heading = sin(state->heading);
    state->cos_heading = cos(state->heading);

    state->sin_half_roll = sin(0.5 * state->roll);
    state->cos_half_roll = cos(0.5 * state->roll);
    state->sin_half_pitch = sin(0.5 * state->pitch);
    state->cos_half_pitch = cos(0.5 * state->pitch);
    state->sin_half_heading = sin(0.5 * state->heading);
    state->cos_half_heading = cos(0.5 * state->heading);

    state->gg = (EA * EGA * state->cos_lat * state->cos_lat + EB * EGB * state->sin_lat * state->sin_lat) /
                sqrt(EA * EA * state->cos_lat * state->cos_lat + EB * EB * state->sin_lat * state->sin_lat)
                * (1.0 - (2.0 / EA) * (1.0 + EF + EWM - 2.0 * EF * state->sin_lat * state->sin_lat) * state->height +
                   3.0 * state->height * state->height / (EA * EA));

    state->w_ie_n[0] = EW * state->cos_lat;
    state->w_ie_n[1] = 0;
    state->w_ie_n[2] = -EW * state->sin_lat;
    state->rm = EA * (1.0 - SQR(EE)) / pow(1.0 - SQR(EE * state->sin_lat), 1.5);
    state->rn = EA / sqrt(1.0 - SQR(EE * state->sin_lat));
    state->w_en_n[0] = state->ve / (state->rn + state->height);
    state->w_en_n[1] = -state->vn / (state->rm + state->height);
    state->w_en_n[2] = -state->ve * state->sin_lat / state->cos_lat / (state->rn + state->height);

    if (last->gg == 0.0 && all_zero(last->w_ie_n, 3) && all_zero(last->w_en_n, 3)) {
        *last = *state;
    }

    state->gg = state->gg + (state->gg - last->gg) / 2;
    for (i = 0; i < 3; ++i) {
        state->w_ie_n[i] = state->w_ie_n[i] + (state->w_ie_n[i] - last->w_ie_n[i]) / 2;
        state->w_en_n[i] = state->w_en_n[i] + (state->w_en_n[i] - last->w_en_n[i]) / 2;
    }

    state->dt = data->sec - (data - 1)->sec;
    return 1;
}

static int update_att_b(INS_data_t *data, INS_result_t *state, double *q_b_b_k_k_1) {

    double phi_k[3];
    double tmp1[3], tmp2[3], tmpd1, tmpd2;

    double theta_k[3] = {data->gx, data->gy, data->gz},
            theta_k_1[3] = {(data - 1)->gx, (data - 1)->gy, (data - 1)->gz};

    cross3(theta_k_1, theta_k, tmp1);
    multiply_number(tmp1, 1.0 / 12.0, 3, 1, tmp2);
    plus(theta_k, tmp2, 3, 1, phi_k);
    tmpd1 = norm(phi_k, 3);
    tmpd2 = sin(0.5 * tmpd1) / (0.5 * tmpd1) * 0.5;
    q_b_b_k_k_1[0] = cos(0.5 * tmpd1);
    multiply_number(phi_k, tmpd2, 3, 1, q_b_b_k_k_1 + 1);

    return 1;
}

static int update_att_n(INS_result_t *state, double *q_n_n_k_k_1) {
    double tmpd1, tmpd2;

    plus(state->w_ie_n, state->w_en_n, 3, 1, state->zeta_k);
    multiply_number(state->zeta_k, state->dt, 3, 1, state->zeta_k);

    tmpd1 = norm(state->zeta_k, 3);
    tmpd2 = -sin(0.5 * tmpd1) / (0.5 * tmpd1) * 0.5;
    q_n_n_k_k_1[0] = cos(0.5 * tmpd1);
//    q_n_n_k_k_1[0] = 1.0 - SQR(tmpd1) / 8.0;
//    tmpd2 = -0.5;
    multiply_number(state->zeta_k, tmpd2, 3, 1, q_n_n_k_k_1 + 1);

    return 1;

}

static int update_att_bn(INS_result_t *state, double *q_b_n_k_1_k_1) {
    q_b_n_k_1_k_1[0] = state->cos_half_heading * state->cos_half_pitch * state->cos_half_roll +
                       state->sin_half_heading * state->sin_half_pitch * state->sin_half_roll;
    q_b_n_k_1_k_1[1] = state->cos_half_heading * state->cos_half_pitch * state->sin_half_roll -
                       state->sin_half_heading * state->sin_half_pitch * state->cos_half_roll;
    q_b_n_k_1_k_1[2] = state->cos_half_heading * state->sin_half_pitch * state->cos_half_roll +
                       state->sin_half_heading * state->cos_half_pitch * state->sin_half_roll;
    q_b_n_k_1_k_1[3] = state->sin_half_heading * state->cos_half_pitch * state->cos_half_roll -
                       state->cos_half_heading * state->sin_half_pitch * state->sin_half_roll;

    return 1;

}

extern int INS_update_att(INS_data_t *data, INS_result_t *state, INS_result_t *last) {
    double q_b_b_k_k_1[4], q_n_n_k_k_1[4], q_b_n_k_1_k_1[4], q_b_b_k_k[4];
    double tmp1[4], tmp2[4];
    //qbb(k,k-1)
    update_att_b(data, state, q_b_b_k_k_1);
    //qnn(k,k-1)
    update_att_n(state, q_n_n_k_k_1);
    //qbn(k-1,k-1)
    update_att_bn(state, q_b_n_k_1_k_1);
    //qbn(k,k)
    quaternion_multiply(q_n_n_k_k_1, q_b_n_k_1_k_1, tmp1);
    quaternion_multiply(tmp1, q_b_b_k_k_1, tmp2);
    normalize(tmp2, 4, q_b_b_k_k);

    state->roll = atan2(2.0 * (q_b_b_k_k[0] * q_b_b_k_k[1] + q_b_b_k_k[2] * q_b_b_k_k[3]),
                        1.0 - 2.0 * (q_b_b_k_k[1] * q_b_b_k_k[1] + q_b_b_k_k[2] * q_b_b_k_k[2]));
    state->pitch = asin(2.0 * (q_b_b_k_k[0] * q_b_b_k_k[2] - q_b_b_k_k[1] * q_b_b_k_k[3]));
    state->heading = atan2(2.0 * (q_b_b_k_k[0] * q_b_b_k_k[3] + q_b_b_k_k[1] * q_b_b_k_k[2]),
                           1.0 - 2.0 * (q_b_b_k_k[2] * q_b_b_k_k[2] + q_b_b_k_k[3] * q_b_b_k_k[3]));
    return 1;
}

static int update_vel_fb(INS_data_t *data, INS_result_t *state, double *dvel_f_b_k_k_1) {
    double dvel_k[3] = {data->ax, data->ay, data->az},
            dvel_k_1[3] = {(data - 1)->ax, (data - 1)->ay, (data - 1)->az},
            theta_k[3] = {data->gx, data->gy, data->gz},
            theta_k_1[3] = {(data - 1)->gx, (data - 1)->gy, (data - 1)->gz};
    double tmp1[3 * 3], tmp2[3], tmp3[3];

    cross3(theta_k, dvel_k, tmp1);
    multiply_number(tmp1, 0.5, 3, 1, tmp1);

    cross3(theta_k_1, dvel_k, tmp2);
    cross3(dvel_k_1, theta_k, tmp3);
    plus(tmp2, tmp3, 3, 1, tmp2);
    multiply_number(tmp2, 1.0 / 12.0, 3, 1, tmp2);

    plus(dvel_k, tmp1, 3, 1, tmp1);
    plus(tmp1, tmp2, 3, 1, dvel_f_b_k_k_1);

    return 1;
}

static int update_vel_bn(INS_result_t *state, double *c_b_n_k_1_k_1) {
    c_b_n_k_1_k_1[MID(3, 0, 0)] = state->cos_pitch * state->cos_heading;
    c_b_n_k_1_k_1[MID(3, 0, 1)] =
            -state->cos_roll * state->sin_heading + state->sin_roll * state->sin_pitch * state->cos_heading;
    c_b_n_k_1_k_1[MID(3, 0, 2)] =
            state->sin_roll * state->sin_heading + state->cos_roll * state->sin_pitch * state->cos_heading;
    c_b_n_k_1_k_1[MID(3, 1, 0)] = state->cos_pitch * state->sin_heading;
    c_b_n_k_1_k_1[MID(3, 1, 1)] =
            state->cos_roll * state->cos_heading + state->sin_roll * state->sin_pitch * state->sin_heading;
    c_b_n_k_1_k_1[MID(3, 1, 2)] =
            -state->sin_roll * state->cos_heading + state->cos_roll * state->sin_pitch * state->sin_heading;
    c_b_n_k_1_k_1[MID(3, 2, 0)] = -state->sin_pitch;
    c_b_n_k_1_k_1[MID(3, 2, 1)] = state->sin_roll * state->cos_pitch;
    c_b_n_k_1_k_1[MID(3, 2, 2)] = state->cos_roll * state->cos_pitch;
    return 1;
}

static int update_vel_fn(INS_result_t *state, double *dvel_f_b_k_k_1, double *c_b_n_k_1_k_1, double *dvel_f_n_k_k) {
    double tmp1[3 * 3];
    double zeta_k_cross[3 * 3];
    cross3_matrix(state->zeta_k, zeta_k_cross);
    multiply_number(zeta_k_cross, 0.5, 3, 3, zeta_k_cross);
    identity_subtract(zeta_k_cross, 1.0, 3, zeta_k_cross);

    multiply(zeta_k_cross, c_b_n_k_1_k_1, 3, 3, 3, tmp1);
    multiply(tmp1, dvel_f_b_k_k_1, 3, 3, 1, dvel_f_n_k_k);

    return 1;
}

static int update_vel_gn(INS_result_t *state, INS_result_t *last, double *dvel_g_n_k_k) {
    double g_ln[3] = {0, 0, state->gg};
    double tmp1[3], tmp2[3];
    double vel_k_1[3] = {last->vn, last->ve, last->vd};

    multiply_number(state->w_ie_n, 2.0, 3, 1, tmp1);
    plus(tmp1, state->w_en_n, 3, 1, tmp1);
    cross3(tmp1, vel_k_1, tmp2);
    subtract(g_ln, tmp2, 3, 1, tmp1);
    multiply_number(tmp1, state->dt, 3, 1, dvel_g_n_k_k);
}

extern int INS_update_vel(INS_data_t *data, INS_result_t *state, INS_result_t *last) {

    double dvel_f_b_k_k_1[3], c_b_n_k_1_k_1[3 * 3], dvel_f_n_k_k[3], dvel_g_n_k_k[3];

    update_vel_fb(data, state, dvel_f_b_k_k_1);
    update_vel_bn(state, c_b_n_k_1_k_1);
    update_vel_fn(state, dvel_f_b_k_k_1, c_b_n_k_1_k_1, dvel_f_n_k_k);
    update_vel_gn(state, last, dvel_g_n_k_k);

    state->vn = last->vn + dvel_f_n_k_k[0] + dvel_g_n_k_k[0];
    state->ve = last->ve + dvel_f_n_k_k[1] + dvel_g_n_k_k[1];
    state->vd = last->vd + dvel_f_n_k_k[2] + dvel_g_n_k_k[2];

    return 1;
}

extern int INS_update_pos(INS_data_t *data, INS_result_t *state, INS_result_t *last) {
    double h_bar, lat_bar, dh, dlat, dlon;
    dh = -0.5 * (state->vd + last->vd) * state->dt;
    h_bar = state->height + 0.5 * dh;
    state->height += dh;
    dlat = 0.5 * (state->vn + last->vn) / (state->rm + h_bar) * state->dt;
    lat_bar = state->lat + 0.5 * dlat;
    state->lat += dlat;
    state->rn = EA / sqrt(1.0 - SQR(EE * sin(lat_bar)));
    dlon = 0.5 * (state->ve + last->ve) / ((state->rn + h_bar) * cos(lat_bar)) * state->dt;
    state->lon += dlon;

    state->sec = data->sec;
    return 1;
}

extern int INS_output(FILE *fp, INS_result_t *state) {
    return fprintf(fp, "%20.10lf %20.10lf %20.10lf %20.10lf %20.10lf %20.10lf %20.10lf %20.10lf %20.10lf %20.10lf\n",
                   state->sec, state->lat / D2R, state->lon / D2R, state->height, state->vn, state->ve, state->vd,
                   state->roll / D2R,
                   state->pitch / D2R, state->heading / D2R);
}

extern int INS_state(FILE *fp, INS_state_t *state) {
    return fprintf(fp, "%20.10lf %20.10lf %20.10lf %20.10lf %20.10lf %20.10lf %20.10lf %20.10lf %20.10lf %20.10lf\n",
                   state->sec, state->lat, state->lon , state->height, state->vn, state->ve, state->vd,
                   state->roll , state->pitch , state->heading );
}