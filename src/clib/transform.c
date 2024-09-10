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
        + powf(cosf(lat1)*sinf(lat2) - sinf(lat1)*cosf(lat2)*cosf(dlon), 2)
    );
    const float x = sinf(lat1)*sinf(lat2) + cosf(lat1)*cosf(lat2)*cosf(dlon);

    return acosf(sinf(lat1)*sinf(lat2) + cosf(lat1)*cosf(lat2)*cosf(dlon));
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
        pow(cos(lat2)*sin(dlon), 2) 
        + pow(cos(lat1)*sin(lat2) - sin(lat1)*cos(lat2)*cos(dlon), 2)
    );
    const double x = sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(dlon);

    return acos(sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(dlon));
    return atan2(y, x);
}

void sph_distance_transform_32bit(
    bool * arr,
    int32_t i_max,
    int32_t j_max,
    float * arr_out
) {
    const float dlon = 2.*PI / ((float) i_max); 
    const float dlat = PI / ((float) j_max); 

    for (int i = 0; i < i_max; i++) {
        for (int j = 0; j < j_max; j++) {
            if (arr[i*j_max + j]) {
                arr_out[i*j_max + j] = 0.;
 
                // determining the patch to work on
                
                // inclusive lower bound index for longitude
                const int ii_min = 0;
                // inclusive upper bound index for longitude
                const int ii_max = i_max-1;
                // inclusive lower bound index for latitude
                const int jj_min = 0;
                // inclusive upper bound index for longitude
                const int jj_max = j_max-1;

                for (int ii = ii_min; ii <= ii_max; ii++) {
                    for (int jj = jj_min; jj <= jj_max; jj++) {
                        if (arr_out[ii*j_max + jj] != 0.f) {
                            const float distance = great_circle_angle_32bit(
                                ((float) i)*dlon,
                                PI_2 - ((float) j)*dlat,    // converting 
                                                            // colatitude 
                                                            // to latitude
                                ((float) ii)*dlon,
                                PI_2 - ((float) jj)*dlat    // converting 
                                                            // colatitude 
                                                            // to latitude
                            );
                            if (
                                (arr_out[ii*j_max + jj] == -1.f) 
                                    || (distance < arr_out[ii*j_max + jj])
                             ) {
                                arr_out[ii*j_max + jj] = distance;
                            }
                        }
                    }
                }
                
            }
        }
    }
}


void sph_distance_transform_64bit(
    bool * arr,
    int64_t i_max,
    int64_t j_max,
    double * arr_out
) {
    const double dlon = 2.*PI / ((double) i_max); 
    const double dlat = PI / ((double) j_max); 

    for (int i = 0; i < i_max; i++) {
        for (int j = 0; j < j_max; j++) {
            if (arr[i*j_max + j]) {
                arr_out[i*j_max + j] = 0.;
 
                // determining the patch to work on
                
                // inclusive lower bound index for longitude
                const int ii_min = 0;
                // inclusive upper bound index for longitude
                const int ii_max = i_max-1;
                // inclusive lower bound index for latitude
                const int jj_min = 0;
                // inclusive upper bound index for longitude
                const int jj_max = j_max-1;

                for (int ii = ii_min; ii <= ii_max; ii++) {
                    for (int jj = jj_min; jj <= jj_max; jj++) {
                        if (arr_out[ii*j_max + jj] != 0.) {
                            const double distance = great_circle_angle_64bit(
                                ((double) i)*dlon,
                                PI_2 - ((double) j)*dlat,   // converting 
                                                            // colatitude 
                                                            // to latitude
                                ((double) ii)*dlon,
                                PI_2 - ((double) jj)*dlat   // converting 
                                                            // colatitude 
                                                            // to latitude
                            );
                            if (
                                (arr_out[ii*j_max + jj] == -1.) 
                                    || (distance < arr_out[ii*j_max + jj])
                             ) {
                                arr_out[ii*j_max + jj] = distance;
                            }
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
 * `sph_distance_transform_32bit_threaded_func()`
 * The parameters are exactly as listed for:
 * `sph_distance_transform_32bit_threaded()`
 */
struct sph_distance_transform_32bit_threaded_args {
    bool *      arr;
    int32_t     i_max;
    int32_t     j_max;
    float *      arr_out;
    int32_t     i_start;        // starting i index for the thread
    int32_t     i_end;          // ending i index for the thread
    int32_t     j_start;        // starting j index for the thread
    int32_t     j_end;          // starting j index for the thread
    int32_t     i_thread;
};

/**
 * (internal)
 * Function pointer for the threaded function:
 * `sph_distance_transform_32bit_threaded_func()`
 * 
 * @param args a pointer to the argument struct
 */
void * sph_distance_transform_32bit_threaded_func(
    void * args
) {
    struct sph_distance_transform_32bit_threaded_args * targs 
    = (struct sph_distance_transform_32bit_threaded_args *) args;

    const float dlon = 2.*PI / ((float) targs->i_max); 
    const float dlat = PI / ((float) targs->j_max);

    for (int i = targs->i_start; i < targs->i_end; i++) {
        for (int j = targs->j_start; j < targs->j_end; j++) {
            if (targs->arr[i*(targs->j_max) + j]) {
                targs->arr_out[i*(targs->j_max) + j] = 0.;
 
                // determining the patch to work on
                
                // inclusive lower bound index for longitude
                const int ii_min = 0;
                // inclusive upper bound index for longitude
                const int ii_max = targs->i_max-1;
                // inclusive lower bound index for latitude
                int jj_min = 0;
                // inclusive upper bound index for longitude
                int jj_max = targs->j_max-1;

                for (int ii = ii_min; ii <= ii_max; ii++) {
                    for (int jj = jj_min; jj <= jj_max; jj++) {
                        if (targs->arr_out[ii*(targs->j_max) + jj] != 0.f) {
                            const float distance = great_circle_angle_32bit(
                                ((float) i)*dlon,
                                PI_2 - ((float) j)*dlat,    // converting 
                                                            // colatitude 
                                                            // to latitude
                                ((float) ii)*dlon,
                                PI_2 - ((float) jj)*dlat    // converting 
                                                            // colatitude 
                                                            // to latitude
                            );

                            if (
                                (targs->arr_out[ii*(targs->j_max) + jj] == -1.f)
                                    || (distance < targs->arr_out[ii*(targs->j_max) + jj]) 
                            ) {
                                targs->arr_out[ii*(targs->j_max) + jj] = distance;
                            }
                        }
                    }
                }
                
            }
        }
    }

    return NULL;
}

void sph_distance_transform_32bit_threaded(
    bool * arr,
    int32_t i_max,
    int32_t j_max,
    float * arr_out,
    int32_t num_threads
) {
    pthread_t tids[num_threads];
    struct sph_distance_transform_32bit_threaded_args targss[num_threads];

    for (int i_thread=0; i_thread < num_threads; i_thread++) {
        targss[i_thread].arr = arr;
        targss[i_thread].i_max = i_max;
        targss[i_thread].j_max = j_max;
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
            sph_distance_transform_32bit_threaded_func, 
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
 * `sph_distance_transform_64bit_threaded_func()`
 * The parameters are exactly as listed for:
 * `sph_distance_transform_64bit_threaded()`
 */
struct sph_distance_transform_64bit_threaded_args {
    bool *      arr;
    int64_t     i_max;
    int64_t     j_max;
    double *    arr_out;
    int64_t     i_start;        // starting i index for the thread
    int64_t     i_end;          // ending i index for the thread
    int64_t     j_start;        // starting j index for the thread
    int64_t     j_end;          // starting j index for the thread
    int64_t     i_thread;
};

/**
 * (internal)
 * Function pointer for the threaded function:
 * `sph_distance_transform_64bit_threaded_func()`
 * 
 * @param args a pointer to the argument struct
 */
void * sph_distance_transform_64bit_threaded_func(
    void * args
) {
    struct sph_distance_transform_64bit_threaded_args * targs 
    = (struct sph_distance_transform_64bit_threaded_args *) args;

    const double dlon = 2.*PI / ((double) targs->i_max); 
    const double dlat = PI / ((double) targs->j_max);

    for (int i = targs->i_start; i < targs->i_end; i++) {
        for (int j = targs->j_start; j < targs->j_end; j++) {
            if (targs->arr[i*(targs->j_max) + j]) {
                targs->arr_out[i*(targs->j_max) + j] = 0.;
 
                // determining the patch to work on
                
                // inclusive lower bound index for longitude
                const int ii_min = 0;
                // inclusive upper bound index for longitude
                const int ii_max = targs->i_max-1;
                // inclusive lower bound index for latitude
                int jj_min = 0;
                // inclusive upper bound index for longitude
                int jj_max = targs->j_max-1;

                for (int ii = ii_min; ii <= ii_max; ii++) {
                    for (int jj = jj_min; jj <= jj_max; jj++) {
                        if (targs->arr_out[ii*(targs->j_max) + jj] != 0.) {
                            const double distance = great_circle_angle_64bit(
                                ((double) i)*dlon,
                                PI_2 - ((double) j)*dlat,   // converting 
                                                            // colatitude 
                                                            // to latitude
                                ((double) ii)*dlon,
                                PI_2 - ((double) jj)*dlat   // converting 
                                                            // colatitude 
                                                            // to latitude
                            );

                            if (
                                (targs->arr_out[ii*(targs->j_max) + jj] == -1.)
                                 || (distance < targs->arr_out[ii*(targs->j_max) + jj]) 
                            ) {
                                targs->arr_out[ii*(targs->j_max) + jj] = distance;
                            }
                        }
                    }
                }
                
            }
        }
    }

    return NULL;
}

void sph_distance_transform_64bit_threaded(
    bool * arr,
    int64_t i_max,
    int64_t j_max,
    double * arr_out,
    int64_t num_threads
) {
    pthread_t tids[num_threads];
    struct sph_distance_transform_64bit_threaded_args targss[num_threads];

    for (int i_thread=0; i_thread < num_threads; i_thread++) {
        targss[i_thread].arr = arr;
        targss[i_thread].i_max = i_max;
        targss[i_thread].j_max = j_max;
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
            sph_distance_transform_64bit_threaded_func, 
            &targss[i_thread]
        );
    }
    
    // waiting for all threads to be done
    for (int i_thread=0; i_thread < num_threads; i_thread++) {
        pthread_join(tids[i_thread], NULL);
    }
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
                                great_circle_angle_32bit(
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
                                great_circle_angle_32bit(
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



