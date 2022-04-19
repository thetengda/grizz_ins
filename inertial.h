//
// Created by tdhuang on 2022/4/12.
//

#ifndef GRIZZ_INS_INERTIAL_H
#define GRIZZ_INS_INERTIAL_H

#include <stdio.h>
#include <stdint.h>
#include <malloc.h>
#include <math.h>

#define INS_MAX_IN_LINE 1024
//长半径
#define EA 6378137.0
//短半径
#define EB 6356752.3141
//偏心率
#define EE 0.0818191910428
//重力常数
#define EGM 3.986005E14
//扁率
#define EF (1 / 298.257222101)
//赤道重力加速度
#define EGA 9.7803267715
//两极重力加速度
#define EGB 9.8321863685
//地球自转角速度
#define EW 0.00007292115

#define EWM (EW*EW*EA*EA*EB/EGM)

#define D2R (M_PI/180.0)

typedef struct INS_data_t {
    double sec, gx, gy, gz, ax, ay, az;
} INS_data_t;

typedef struct INS_state_t {
    double sec, lat, lon, height, vn, ve, vd, roll, pitch, heading;
} INS_state_t;
typedef struct INS_result_t {
    //state
    double sec, lat, lon, height, vn, ve, vd, roll, pitch, heading;
    //auxiliary
    double sin_lat, cos_lat, sin_lon, cos_lon, gg, dt;
    double sin_roll, cos_roll, sin_pitch, cos_pitch, sin_heading, cos_heading;
    double sin_half_roll, cos_half_roll, sin_half_pitch, cos_half_pitch, sin_half_heading, cos_half_heading;
    double   rm, rn, w_ie_n[3], w_en_n[3], zeta_k[3];
//    double vel_k_1[3];
} INS_result_t;

/// decode binary INS data
/// \param fp
/// \param data
/// \return
extern size_t INS_decode_binary(const char *fp, void **data);

/// update gg, and others
/// \param state
/// \return
extern int INS_update_param(INS_data_t *data, INS_result_t *state, INS_result_t *last);

/// attitude update
/// \param data
/// \param state
/// \return
extern int INS_update_att(INS_data_t *data, INS_result_t *state, INS_result_t *last);

/// velocity update
/// \param data
/// \param state
/// \return
extern int INS_update_vel(INS_data_t *data, INS_result_t *state, INS_result_t *last);

/// position update
/// \param data
/// \param state
/// \return
extern int INS_update_pos(INS_data_t *data, INS_result_t *state, INS_result_t *last);

/// output result
/// \param fp
/// \param state
/// \return
extern int INS_output(FILE *fp, INS_result_t *state);

extern int INS_state(FILE *fp, INS_state_t *state);

#endif //GRIZZ_INS_INERTIAL_H
