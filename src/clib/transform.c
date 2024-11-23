/**
 * @file transform.c
 * @author Pejvak Javaheri; pejvak.javaheri@mail.utoronto.ca
 * @brief Definitions for spherical distance transforms
 */

#include "transform.h"

void single_plate_interior_distance_transform_64bit(
    double *    xs,
    double *    ys,
    double *    zs,
    int64_t     num_points,
    double      R,
    double *    arr_out
) {
    // first, only working with Cartesian distance
    const double ONE_OVER_TWO_R_SQUARED = 1./(2.*R*R);
    const double MAX_CARTESIAN_DISTANCE = 4.*R*R;
    
    for (int i = 0; i < num_points; i++) {
        if (arr_out[i] == -1.) {
            arr_out[i] = MAX_CARTESIAN_DISTANCE;
            for (int j = 0; j < num_points; j++) {
                if (arr_out[j] == -2.) {
                    // the point [j] is not on the plate

                    // crude checking whether to even bother with
                    // spherical distance
                    const double dx = xs[j] - xs[i];
                    const double dy = ys[j] - ys[i];
                    const double dz = zs[j] - zs[i];
                    const double d2 = dx*dx + dy*dy + dz*dz;
                    if (d2 < arr_out[i]) {
                        arr_out[i] = d2;
                    }
                }
            }
            // converting the minimum Cartesian distance to geodesic
            arr_out[i] = acos(1. - ONE_OVER_TWO_R_SQUARED*arr_out[i]);
        }    
    }
}

void full_plate_interior_distance_transform_64bit(
    double *    xs,
    double *    ys,
    double *    zs,
    int64_t *   plate_IDs,
    int64_t     num_points,
    double      R,
    double *    arr_out
) {
    // first, only working with Cartesian distance
    const double ONE_OVER_TWO_R_SQUARED = 1./(2.*R*R);
    const double MAX_CARTESIAN_DISTANCE = 4.*R*R;
    
    for (int i = 0; i < num_points; i++) {
        arr_out[i] = MAX_CARTESIAN_DISTANCE;
        for (int j = 0; j < num_points; j++) {
            if (plate_IDs[j] != -plate_IDs[i]) {
                // the point [j] is not on the plate

                // crude checking whether to even bother with
                // spherical distance
                const double dx = xs[j] - xs[i];
                const double dy = ys[j] - ys[i];
                const double dz = zs[j] - zs[i];
                const double d2 = dx*dx + dy*dy + dz*dz;
                if (d2 < arr_out[i]) {
                    arr_out[i] = d2;
                }
            }
        }
        // converting the minimum Cartesian distance to geodesic
        arr_out[i] = acos(1. - ONE_OVER_TWO_R_SQUARED*arr_out[i]);
    }    
}

/**
 * (internal)
 * Argument struct for the threaded function:
 * `single_plate_interior_distance_transform_64bit_threaded_func`
 * The parameters are exactly as listed for:
 * `full_plate_interior_distance_transform_64bit_threaded_func`
 * 
 * @warning: `plate_IDs` is only used for 
 *           `full_plate_interior_distance_transform_64bit_threaded_func`
 */
struct plate_interior_distance_transform_64bit_threaded_args {
    double *    xs;
    double *    ys;
    double *    zs;
    int64_t *   plate_IDs;
    int64_t     num_points;
    double      R;
    double *    arr_out;
    int64_t     num_threads;
    int64_t     i_thread;
};

/**
 * (internal)
 * Threaded function for `single_plate_interior_distance_transform_64bit_threaded`
 * 
 * @param args the argument struct containing function parameters
 */
void * single_plate_interior_distance_transform_64bit_threaded_func(
    void * args
) {
    /**
     * Assumes that plate_IDs are initialized with -1. and 0. :
     *      -> -1. corresponds to the plate of interest
     *      -> -2. corresponds to other plates
     */

    struct plate_interior_distance_transform_64bit_threaded_args * targs =
        (struct plate_interior_distance_transform_64bit_threaded_args *) args;
    

    int i_start = (targs->num_points / targs->num_threads)*(targs->i_thread);
    int i_end   = (targs->num_points / targs->num_threads)*(targs->i_thread +1);

    // to prevent missing out due to division round downs
    if (targs->i_thread == targs->num_threads - 1) {
        i_end = targs->num_points;
    }

    // first, only working with Cartesian distance
    const double ONE_OVER_TWO_R_SQUARED = 1./(2. * targs->R * targs->R);
    const double MAX_CARTESIAN_DISTANCE = 4. * targs->R * targs->R;
    
    for (int i = 0; i < targs->num_points; i++) {
        if (targs->arr_out[i] == -1.) {
            targs->arr_out[i] = MAX_CARTESIAN_DISTANCE;
            for (int j = 0; j < targs->num_points; j++) {
                if (targs->arr_out[j] == -2.) {
                    // the point [j] is not on the plate

                    // crude checking whether to even bother with
                    // spherical distance
                    const double dx = targs->xs[j] - targs->xs[i];
                    const double dy = targs->ys[j] - targs->ys[i];
                    const double dz = targs->zs[j] - targs->zs[i];
                    const double d2 = dx*dx + dy*dy + dz*dz;
                    if (d2 < targs->arr_out[i]) {
                        targs->arr_out[i] = d2;
                    }
                }
            }
            // converting the minimum Cartesian distance to geodesic
            targs->arr_out[i] = acos(
                1. - ONE_OVER_TWO_R_SQUARED * targs->arr_out[i]
            );
        }    
    }
    return NULL;
}


void single_plate_interior_distance_transform_64bit_threaded(
    double *    xs,
    double *    ys,
    double *    zs,
    int64_t     num_points,
    double      R,
    double *    arr_out,
    int64_t     num_threads
) {
    pthread_t tids[num_threads];
    struct plate_interior_distance_transform_64bit_threaded_args targss[num_threads];

    for (int i_thread=0; i_thread < num_threads; i_thread++) {
        targss[i_thread].xs             = xs;
        targss[i_thread].ys             = ys;
        targss[i_thread].zs             = zs;
        targss[i_thread].num_points     = num_points;
        targss[i_thread].R              = R;
        targss[i_thread].arr_out        = arr_out;
        targss[i_thread].num_threads    = num_threads;
        targss[i_thread].i_thread       = i_thread;

        // creating threads
        pthread_create(
            &tids[i_thread], 
            NULL, 
            single_plate_interior_distance_transform_64bit_threaded_func, 
            &targss[i_thread]
        );
    }
    
    // waiting for all threads to be done
    for (int i_thread=0; i_thread < num_threads; i_thread++) {
        pthread_join(tids[i_thread], NULL);
    }
}

/**
 * (internal)
 * Threaded function for `full_plate_interior_distance_transform_64bit_threaded`
 * 
 * @param args the argument struct containing function parameters
 */
void * full_plate_interior_distance_transform_64bit_threaded_func(
    void * args
) {
    /**
     * Assumes that plate_IDs are initialized with -1. and 0. :
     *      -> -1. corresponds to the plate of interest
     *      -> -2. corresponds to other plates
     */

    struct plate_interior_distance_transform_64bit_threaded_args * targs =
        (struct plate_interior_distance_transform_64bit_threaded_args *) args;
    

    int i_start = (targs->num_points / targs->num_threads)*(targs->i_thread);
    int i_end   = (targs->num_points / targs->num_threads)*(targs->i_thread +1);

    // to prevent missing out due to division round downs
    if (targs->i_thread == targs->num_threads - 1) {
        i_end = targs->num_points;
    }

    // first, only working with Cartesian distance
    const double ONE_OVER_TWO_R_SQUARED = 1./(2. * targs->R * targs->R);
    const double MAX_CARTESIAN_DISTANCE = 4. * targs->R * targs->R;
    
    for (int i = 0; i < targs->num_points; i++) {
        targs->arr_out[i] = MAX_CARTESIAN_DISTANCE;
        for (int j = 0; j < targs->num_points; j++) {
            if (targs->plate_IDs[j] != targs->plate_IDs[i]) {
                // the point [j] is not on the plate

                // crude checking whether to even bother with
                // spherical distance
                const double dx = targs->xs[j] - targs->xs[i];
                const double dy = targs->ys[j] - targs->ys[i];
                const double dz = targs->zs[j] - targs->zs[i];
                const double d2 = dx*dx + dy*dy + dz*dz;
                if (d2 < targs->arr_out[i]) {
                    targs->arr_out[i] = d2;
                }
            }
        }
        // converting the minimum Cartesian distance to geodesic
        targs->arr_out[i] = acos(
            1. - ONE_OVER_TWO_R_SQUARED * targs->arr_out[i]
        );
    }
    return NULL;
}

void full_plate_interior_distance_transform_64bit_threaded(
    double *    xs,
    double *    ys,
    double *    zs,
    int64_t *   plate_IDs,
    int64_t     num_points,
    double      R,
    double *    arr_out,
    int64_t     num_threads
) {
    pthread_t tids[num_threads];
    struct plate_interior_distance_transform_64bit_threaded_args targss[num_threads];

    for (int i_thread=0; i_thread < num_threads; i_thread++) {
        targss[i_thread].xs             = xs;
        targss[i_thread].ys             = ys;
        targss[i_thread].zs             = zs;
        targss[i_thread].plate_IDs      = plate_IDs;
        targss[i_thread].num_points     = num_points;
        targss[i_thread].R              = R;
        targss[i_thread].arr_out        = arr_out;
        targss[i_thread].num_threads    = num_threads;
        targss[i_thread].i_thread       = i_thread;

        // creating threads
        pthread_create(
            &tids[i_thread], 
            NULL, 
            full_plate_interior_distance_transform_64bit_threaded_func, 
            &targss[i_thread]
        );
    }
    
    // waiting for all threads to be done
    for (int i_thread=0; i_thread < num_threads; i_thread++) {
        pthread_join(tids[i_thread], NULL);
    }
}

int fused_distance_threshold_transform_32bit(
    float *     xs,
    float *     ys,
    float *     zs,
    bool *      arr,
    int32_t     num_points,
    float       R,
    float       threshold,
    bool *      arr_out
) {
    // first, only working with Cartesian distance
    const float THRESHOLD_D_SQUARED = 2.*R*R*(1. - cosf(threshold));
    
    for (int i = 0; i < num_points; i++) {
        if (arr[i]) {
            arr_out[i] = true;
            for (int j = 0; j < num_points; j++) {
                if (!arr_out[j]) {
                    // the point [j] is not on the plate

                    // crude checking whether to even bother with
                    // spherical distance
                    const double dx = xs[j] - xs[i];
                    const double dy = ys[j] - ys[i];
                    const double dz = zs[j] - zs[i];
                    const double d2 = dx*dx + dy*dy + dz*dz;
                    
                    arr_out[j] = d2 < THRESHOLD_D_SQUARED;
                }
            }
        }    
    }
    return 0;
}

int fused_distance_threshold_transform_64bit(
    double *    xs,
    double *    ys,
    double *    zs,
    bool *      arr,
    int64_t     num_points,
    double      R,
    double       threshold,
    bool *      arr_out
) {
    // first, only working with Cartesian distance
    const double THRESHOLD_D_SQUARED = 2.*R*R*(1. - cos(threshold));
    const double MAX_CARTESIAN_DISTANCE = 4.*R*R;
    
    for (int i = 0; i < num_points; i++) {
        if (arr[i]) {
            arr_out[i] = true;
            for (int j = 0; j < num_points; j++) {
                if (!arr_out[j]) {
                    // the point [j] is not on the plate

                    // crude checking whether to even bother with
                    // spherical distance
                    const double dx = xs[j] - xs[i];
                    const double dy = ys[j] - ys[i];
                    const double dz = zs[j] - zs[i];
                    const double d2 = dx*dx + dy*dy + dz*dz;

                    arr_out[j] = d2 < THRESHOLD_D_SQUARED;
                }
            }
        }    
    }
    return 0;
}

int gridded_fused_distance_threshold_transform_32bit(
    float *     xs,
    float *     ys,
    float *     zs,
    bool *      arr,
    int32_t     i_max,
    int32_t     j_max,
    float       R,
    float       threshold,
    bool *      arr_out
) {
    // first, only working with Cartesian distance
    const float THRESHOLD_D_SQUARED = 2.*R*R*(1. - cosf(threshold));
    
    const float dlat = PI / ((float) i_max); 
    const int di = (int) (threshold / dlat) + 2;        // +2 for a conservative

    for (int i = 0; i < i_max; i++) {
        for (int j = 0; j < j_max; j++) {
            const int ref_index = i*j_max + j;
            if (arr[ref_index]) {
                arr_out[ref_index] = true;

                int i_start = 0;
                int i_end   = i_max-1;

                if (i - di > 0)         i_start = i - di;
                if (i + di < i_max-1)   i_end   = i + di;

                for (int ii = i_start; ii <= i_end; ii++) {
                    for (int jj = 0; jj < j_max; jj++) {
                        const int index = ii*j_max + jj;
                        if (!arr_out[index]) {
                            // the point [j] is not on the plate

                            // crude checking whether to even bother with
                            // spherical distance
                            const float dx = xs[index] - xs[ref_index];
                            const float dy = ys[index] - ys[ref_index];
                            const float dz = zs[index] - zs[ref_index];
                            const float d2 = dx*dx + dy*dy + dz*dz;
                            
                            arr_out[index] = d2 < THRESHOLD_D_SQUARED;
                        }
                    }
                }
            }    
        }
    }
    return 0;
}

int gridded_fused_distance_threshold_transform_64bit(
    double *    xs,
    double *    ys,
    double *    zs,
    bool *      arr,
    int64_t     i_max,
    int64_t     j_max,
    double      R,
    double      threshold,
    bool *      arr_out
) {
    // first, only working with Cartesian distance
    const double THRESHOLD_D_SQUARED = 2.*R*R*(1. - cos(threshold));
    
    const double dlat = PI / ((double) i_max); 
    const int di = (int) (threshold / dlat) + 2;        // +2 for a conservative

    for (int i = 0; i < i_max; i++) {
        for (int j = 0; j < j_max; j++) {
            const int ref_index = i*j_max + j;
            if (arr[ref_index]) {
                arr_out[ref_index] = true;

                int i_start = 0;
                int i_end   = i_max-1;

                if ((i - di) > 0)         i_start = i - di;
                if ((i + di) < i_max-1)   i_end   = i + di;

                for (int ii = i_start; ii <= i_end; ii++) {
                    for (int jj = 0; jj < j_max; jj++) {
                        const int index = ii*j_max + jj;
                        if (!arr_out[index]) {
                            // the point [j] is not on the plate

                            // crude checking whether to even bother with
                            // spherical distance
                            const double dx = xs[index] - xs[ref_index];
                            const double dy = ys[index] - ys[ref_index];
                            const double dz = zs[index] - zs[ref_index];
                            const double d2 = dx*dx + dy*dy + dz*dz;
                            
                            arr_out[index] = d2 < THRESHOLD_D_SQUARED;
                        }
                    }
                }
            }    
        }
    }
    return 0;
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//                    LEGACY CODE
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

double great_circle_angle_64bit(
    double lon1,
    double lat1,
    double lon2,
    double lat2
) {
    const double dlon = fabs(lon2 - lon1);
    
    /*
    const double y = sqrt(
        pow(cos(lat2)*sin(dlon), 2) 
        + pow(cos(lat1)*sin(lat2) - sin(lat1)*cos(lat2)*cos(dlon), 2)
    );
    const double x = sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(dlon);
    */

    // faster than individual calls
    double sin_lat1, cos_lat1, sin_lat2, cos_lat2, sin_dlon, cos_dlon;

    sincos(dlon, &sin_dlon, &cos_dlon);
    sincos(lat1, &sin_lat1, &cos_lat1);
    sincos(lat2, &sin_lat2, &cos_lat2);

    const double y = sqrt(
        pow(cos_lat2*sin_dlon, 2) 
        + pow(cos_lat1*sin_lat2 - sin_lat1*cos_lat2*cos_dlon, 2)
    );
    const double x = sin_lat1*sin_lat2 + cos_lat1*cos_lat2*cos_dlon;

    //return acos(sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(dlon));
    return atan2(y, x);
}


void sph_fused_distance_threshold_transform_32bit(
    bool * arr,
    int32_t i_max,
    int32_t j_max,
    float threshold,
    bool * arr_out
) {
    const float dlon = 2.*PI / ((float) i_max); 
    const float dlat = PI / ((float) j_max); 

    const int di = (int) (threshold / dlon) + 2;        // +2 for a conservative
    const int dj = (int) (threshold / dlat) + 2;        // estimate

    for (int i = 0; i < i_max; i++) {
        for (int j = 0; j < j_max; j++) {
            if (arr[i*j_max + j]) {
                arr_out[i*j_max + j] = true;
 
                // determining the patch to work on
                
                // inclusive lower bound index for longitude
                const int ii_min = 0;
                // inclusive upper bound index for longitude
                const int ii_max = i_max-1;
                // inclusive lower bound index for latitude
                int jj_min = 0;
                // inclusive upper bound index for longitude
                int jj_max = j_max-1;

                // whether to omit far away regions
                if (0 <= j - dj) {
                    jj_min = j - dj;
                } 
                if (j + dj < j_max) {
                    jj_max= j + dj;
                }

                for (int ii = ii_min; ii <= ii_max; ii++) {
                    for (int jj = jj_min; jj <= jj_max; jj++) {
                        if (
                            (!arr_out[ii*j_max + jj])
                            && (
                                great_circle_angle_64bit(
                                    ((float) i)*dlon,
                                    PI_2 - ((float) j)*dlat,    // converting 
                                                                // colatitude 
                                                                // to latitude
                                    ((float) ii)*dlon,
                                    PI_2 - ((float) jj)*dlat    // converting 
                                                                // colatitude 
                                                                // to latitude
                                ) < threshold
                            )
                        ) {
                            arr_out[ii*j_max + jj] = true;
                        }
                    }
                }
                
            }
        }
    }
}

void sph_fused_distance_threshold_transform_64bit(
    bool * arr,
    int64_t i_max,
    int64_t j_max,
    double threshold,
    bool * arr_out
) {
    const double dlon = 2.*PI / ((double) i_max); 
    const double dlat = PI / ((double) j_max); 

    const int di = (int) (threshold / dlon) + 2;        // +2 for a conservative
    const int dj = (int) (threshold / dlat) + 2;        // estimate

    for (int i = 0; i < i_max; i++) {
        for (int j = 0; j < j_max; j++) {
            if (arr[i*j_max + j]) {
                arr_out[i*j_max + j] = true;
 
                // determining the patch to work on
                
                // inclusive lower bound index for longitude
                const int ii_min = 0;
                // inclusive upper bound index for longitude
                const int ii_max = i_max-1;
                // inclusive lower bound index for latitude
                int jj_min = 0;
                // inclusive upper bound index for longitude
                int jj_max = j_max-1;

                // whether to omit far away regions
                if (0 <= j - dj) {
                    jj_min = j - dj;
                } 
                if (j + dj < j_max) {
                    jj_max= j + dj;
                }

                for (int ii = ii_min; ii <= ii_max; ii++) {
                    for (int jj = jj_min; jj <= jj_max; jj++) {
                        if (
                            (!arr_out[ii*j_max + jj])
                            && (
                                great_circle_angle_64bit(
                                    ((double) i)*dlon,
                                    PI_2 - ((double) j)*dlat,   // converting 
                                                                // colatitude 
                                                                // to latitude
                                    ((double) ii)*dlon,
                                    PI_2 - ((double) jj)*dlat   // converting 
                                                                // colatitude 
                                                                // to latitude
                                ) < threshold
                            )
                        ) {
                            arr_out[ii*j_max + jj] = true;
                        }
                    }
                }
                
            }
        }
    }
}

/**
 * (internal)
 * Argument struct for the threaded function:
 * `sph_fused_distance_threshold_transform_32bit_threaded_func()`
 * The parameters are exactly as listed for:
 * `sph_fused_distance_threshold_transform_32bit_threaded()`
 */
struct sph_fused_distance_threshold_transform_32bit_threaded_args {
    bool *      arr;
    int32_t     i_max;
    int32_t     j_max;
    float       threshold;
    bool *      arr_out;
    int32_t     i_start;        // starting i index for the thread
    int32_t     i_end;          // ending i index for the thread
    int32_t     j_start;        // starting j index for the thread
    int32_t     j_end;          // starting j index for the thread
    int32_t     i_thread;
};

/**
 * (internal)
 * Function pointer for the threaded function:
 * `sph_fused_distance_threshold_transform_32bit_threaded_func()`
 * 
 * @param args a pointer to the argument struct
 */
void * sph_fused_distance_threshold_transform_32bit_threaded_func(
    void * args
) {
    struct sph_fused_distance_threshold_transform_32bit_threaded_args * targs 
    = (struct sph_fused_distance_threshold_transform_32bit_threaded_args *) args;

    const float dlon = 2.*PI / ((float) targs->i_max); 
    const float dlat = PI / ((float) targs->j_max); 

    const int di = (int) (targs->threshold / dlon) + 2;        // +2 for a conservative
    const int dj = (int) (targs->threshold / dlat) + 2;        // estimate

    for (int i = targs->i_start; i < targs->i_end; i++) {
        for (int j = targs->j_start; j < targs->j_end; j++) {
            if (targs->arr[i*(targs->j_max) + j]) {
                targs->arr_out[i*(targs->j_max) + j] = true;
 
                // determining the patch to work on
                
                // inclusive lower bound index for longitude
                const int ii_min = 0;
                // inclusive upper bound index for longitude
                const int ii_max = targs->i_max-1;
                // inclusive lower bound index for latitude
                int jj_min = 0;
                // inclusive upper bound index for longitude
                int jj_max = targs->j_max-1;

                // whether to omit far away regions
                if (0 <= j - dj) {
                    jj_min = j - dj;
                } 
                if (j + dj < targs->j_max) {
                    jj_max= j + dj;
                }

                for (int ii = ii_min; ii <= ii_max; ii++) {
                    for (int jj = jj_min; jj <= jj_max; jj++) {
                        if (
                            (!targs->arr_out[ii*(targs->j_max) + jj])
                            && (
                                great_circle_angle_64bit(
                                    ((float) i)*dlon,
                                    PI_2 - ((float) j)*dlat,    // converting 
                                                                // colatitude 
                                                                // to latitude
                                    ((float) ii)*dlon,
                                    PI_2 - ((float) jj)*dlat    // converting 
                                                                // colatitude 
                                                                // to latitude
                                ) < targs->threshold
                            )
                        ) {
                            targs->arr_out[ii*(targs->j_max) + jj] = true;
                        }
                    }
                }
                
            }
        }
    }

    return NULL;
}

void sph_fused_distance_threshold_transform_32bit_threaded(
    bool * arr,
    int32_t i_max,
    int32_t j_max,
    float threshold,
    bool * arr_out,
    int32_t num_threads
) {
    pthread_t tids[num_threads];
    struct sph_fused_distance_threshold_transform_32bit_threaded_args targss[num_threads];

    for (int i_thread=0; i_thread < num_threads; i_thread++) {
        targss[i_thread].arr = arr;
        targss[i_thread].i_max = i_max;
        targss[i_thread].j_max = j_max;
        targss[i_thread].threshold = threshold;
        targss[i_thread].arr_out = arr_out;
        targss[i_thread].i_start = 0;
        targss[i_thread].i_end = i_max;

        // division of labour is more efficient to be done row by ro since
        // the array is row major
        targss[i_thread].j_start = j_max / num_threads * i_thread;
        targss[i_thread].j_end = j_max / num_threads * (i_thread+1);
        targss[i_thread].i_thread = i_thread;

        // ensuring the whole array is covered
        if (i_thread == num_threads - 1) {
            targss[num_threads-1].j_end = j_max;
        }

        // creating threads
        pthread_create(
            &tids[i_thread], 
            NULL, 
            sph_fused_distance_threshold_transform_32bit_threaded_func, 
            &targss[i_thread]
        );
    }
    
    // waiting for all threads to be done
    for (int i_thread=0; i_thread < num_threads; i_thread++) {
        pthread_join(tids[i_thread], NULL);
    }
}

/**
 * (internal)
 * Argument struct for the threaded function:
 * `sph_fused_distance_threshold_transform_64bit_threaded_func()`
 * The parameters are exactly as listed for:
 * `sph_fused_distance_threshold_transform_64bit_threaded()`
 */
struct sph_fused_distance_threshold_transform_64bit_threaded_args {
    bool *      arr;
    int64_t     i_max;
    int64_t     j_max;
    double       threshold;
    bool *      arr_out;
    int64_t     i_start;        // starting i index for the thread
    int64_t     i_end;          // ending i index for the thread
    int64_t     j_start;        // starting j index for the thread
    int64_t     j_end;          // starting j index for the thread
    int64_t     i_thread;
};

/**
 * (internal)
 * Function pointer for the threaded function:
 * `sph_fused_distance_threshold_transform_64bit_threaded_func()`
 * 
 * @param args a pointer to the argument struct
 */
void * sph_fused_distance_threshold_transform_64bit_threaded_func(
    void * args
) {
    struct sph_fused_distance_threshold_transform_64bit_threaded_args * targs 
    = (struct sph_fused_distance_threshold_transform_64bit_threaded_args *) args;

    const double dlon = 2.*PI / ((double) targs->i_max); 
    const double dlat = PI / ((double) targs->j_max); 

    const int di = (int) (targs->threshold / dlon) + 2;        // +2 for a conservative
    const int dj = (int) (targs->threshold / dlat) + 2;        // estimate

    for (int i = targs->i_start; i < targs->i_end; i++) {
        for (int j = targs->j_start; j < targs->j_end; j++) {
            if (targs->arr[i*(targs->j_max) + j]) {
                targs->arr_out[i*(targs->j_max) + j] = true;
 
                // determining the patch to work on
                
                // inclusive lower bound index for longitude
                const int ii_min = 0;
                // inclusive upper bound index for longitude
                const int ii_max = targs->i_max-1;
                // inclusive lower bound index for latitude
                int jj_min = 0;
                // inclusive upper bound index for longitude
                int jj_max = targs->j_max-1;

                // whether to omit far away regions
                if (0 <= j - dj) {
                    jj_min = j - dj;
                } 
                if (j + dj < targs->j_max) {
                    jj_max= j + dj;
                }

                for (int ii = ii_min; ii <= ii_max; ii++) {
                    for (int jj = jj_min; jj <= jj_max; jj++) {
                        if (
                            (!targs->arr_out[ii*(targs->j_max) + jj])
                            && (
                                great_circle_angle_64bit(
                                    ((double) i)*dlon,
                                    PI_2 - ((double) j)*dlat,   // converting 
                                                                // colatitude 
                                                                // to latitude
                                    ((double) ii)*dlon,
                                    PI_2 - ((double) jj)*dlat   // converting 
                                                                // colatitude 
                                                                // to latitude
                                ) < targs->threshold
                            )
                        ) {
                            targs->arr_out[ii*(targs->j_max) + jj] = true;
                        }
                    }
                }
                
            }
        }
    }
    return NULL;
}

void sph_fused_distance_threshold_transform_64bit_threaded(
    bool * arr,
    int64_t i_max,
    int64_t j_max,
    double threshold,
    bool * arr_out,
    int64_t num_threads
) {
    pthread_t tids[num_threads];
    struct sph_fused_distance_threshold_transform_64bit_threaded_args targss[num_threads];

    for (int i_thread=0; i_thread < num_threads; i_thread++) {
        targss[i_thread].arr = arr;
        targss[i_thread].i_max = i_max;
        targss[i_thread].j_max = j_max;
        targss[i_thread].threshold = threshold;
        targss[i_thread].arr_out = arr_out;
        targss[i_thread].i_start = 0;
        targss[i_thread].i_end = i_max;

        // division of labour is more efficient to be done row by ro since
        // the array is row major
        targss[i_thread].j_start = j_max / num_threads * i_thread;
        targss[i_thread].j_end = j_max / num_threads * (i_thread+1);
        targss[i_thread].i_thread = i_thread;

        // ensuring the whole array is covered
        if (i_thread == num_threads - 1) {
            targss[num_threads-1].j_end = j_max;
        }

        // creating threads
        pthread_create(
            &tids[i_thread], 
            NULL, 
            sph_fused_distance_threshold_transform_64bit_threaded_func, 
            &targss[i_thread]
        );
    }
    
    // waiting for all threads to be done
    for (int i_thread=0; i_thread < num_threads; i_thread++) {
        pthread_join(tids[i_thread], NULL);
    }
}
