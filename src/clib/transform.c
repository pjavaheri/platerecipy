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