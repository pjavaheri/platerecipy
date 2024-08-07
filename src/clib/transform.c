/**
 * @file transform.c
 * @author Pejvak Javaheri; pejvak.javaheri@mail.utoronto.ca
 * @brief Definitions for image transformations
 */

#include "transform.h"

float great_circle_angle_32bit(
    float lon1,
    float lat1,
    float lon2,
    float lat2
) {
    const float dlon = fabsf(lon2 - lon1);
    const float y = sqrtf(
        powf(cosf(lat2)*sinf(dlon), 2) 
        + powf(cosf(lat1)*sinf(lat2)-sinf(lat1)*cosf(lat2)*cosf(dlon), 2)
    );
    const float x = sinf(lat1)*sinf(lat2) + cosf(lat1)*cosf(lat2)*cosf(dlon);

    return atan2f(y, x);
}

double great_circle_angle_64bit(
    double lon1,
    double lat1,
    double lon2,
    double lat2
) {
    const double dlon = fabs(lon2 - lon1);
    const double y = sqrt(
        powf(cos(lat2)*sin(dlon), 2) 
        + powf(cos(lat1)*sin(lat2)-sin(lat1)*cos(lat2)*cos(dlon), 2)
    );
    const double x = sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(dlon);

    return atan2(y, x);
}

void fused_distance_threshold_transform_32bit(
    bool * arr,
    int32_t i_max,
    int32_t j_max,
    float threshold,
    bool * arr_out
) {
    const float dlon = 2.*PI / ((float) i_max); 
    const float dlat = PI / ((float) j_max); 

    const int di = (int) threshold / dlon;
    const int dj = (int) threshold / dlat;

    for (int i = 0; i < i_max; i++) {
        for (int j = 0; j < j_max; j++) {
            if (arr[i*j_max + j]) {
                // if the neighborhood is nice! 
                if (((di/2 + i < i_max) && (0 <= i - di/2))
                    && ((dj/2 + j < j_max) && (0 <= j - dj/2))) {
                    for (int ii = i - di/2; ii <= i + di/2; ii++) {
                        for (int jj = j - dj/2; jj <= j + dj/2; jj++) {
                            if (
                                (!arr_out[ii*j_max + jj])
                                && great_circle_angle_32bit(
                                    i*dlon,
                                    PI_2 - j*dlat,  // converting colatitude to
                                                    // latitude
                                    ii*dlon,
                                    PI_2 - jj*dlat  // converting colatitude to
                                                    // latitude
                                ) < threshold
                            ) {
                                arr_out[ii*j_max + jj] = true;
                            }
                        }
                    }
                }
            } else {
                arr_out[i*j_max + j] = false;
            }
        }
    }
}

void fused_distance_threshold_transform_64bit(
    bool * arr,
    int64_t i_max,
    int64_t j_max,
    double threshold,
    bool * arr_out
) {
    const double dlon = 2.*PI / ((double) i_max); 
    const double dlat = PI / ((double) j_max); 

    const int di = (int) threshold / dlon;
    const int dj = (int) threshold / dlat;

    for (int i = 0; i < i_max; i++) {
        for (int j = 0; j < j_max; j++) {
            if (arr[i*j_max + j]) {
                // if the neighborhood is nice! 
                if (((di/2 + i < i_max) && (0 <= i - di/2))
                    && ((dj/2 + j < j_max) && (0 <= j - dj/2))) {
                    for (int ii = i - di/2; ii <= i + di/2; ii++) {
                        for (int jj = j - dj/2; jj <= j + dj/2; jj++) {
                            if (
                                (!arr_out[ii*j_max + jj])
                                && great_circle_angle_64bit(
                                    i*dlon,
                                    PI_2 - j*dlat,  // converting colatitude to
                                                    // latitude
                                    ii*dlon,
                                    PI_2 - jj*dlat  // converting colatitude to
                                                    // latitude
                                ) < threshold
                            ) {
                                arr_out[ii*j_max + jj] = true;
                            }
                        }
                    }
                }
            } else {
                arr_out[i*j_max + j] = false;
            }
        }
    }
}