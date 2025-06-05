/**
 * @file transform.c
 * @author Pejvak Javaheri; pejvak.javaheri@mail.utoronto.ca
 * @brief Definitions for spherical distance transforms
 */

#include "transform.h"

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//      single_plate_interior_distance_transform
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/**
 * The following set of functions are to be used when the spherical distance 
 * transform is to be applied on the interior of only a single plates, indicated 
 * by the initial -1 values for plate interior and -2 values for point not on 
 * the plate.
 */


int single_plate_interior_distance_transform_64bit(
    double *    xs,
    double *    ys,
    double *    zs,
    int32_t     num_points,
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
    return 0;
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
    int32_t *   plate_IDs;
    int32_t     num_points;
    double      R;
    double *    arr_out;
    int32_t     num_threads;
    int32_t     i_thread;
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
    
    for (int i = i_start; i < i_end; i++) {
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


int single_plate_interior_distance_transform_64bit_threaded(
    double *    xs,
    double *    ys,
    double *    zs,
    int32_t     num_points,
    double      R,
    double *    arr_out,
    int32_t     num_threads
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

    return 0;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//       full_plate_interior_distance_transform
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/**
 * The following set of functions are to be used when the spherical distance 
 * transform is to be applied on the interior of all plates, indicated in the
 * `plate_IDs` array.
 */

int full_plate_interior_distance_transform_64bit(
    double *    xs,
    double *    ys,
    double *    zs,
    int32_t *   plate_IDs,
    int32_t     num_points,
    double      R,
    double *    arr_out
) {
    FILE * f = fopen("stdout.txt", "w");
    fprintf(f, "in C: num_points = %d\n", num_points);
    fclose(f);

    // first, only working with Cartesian distance
    const double ONE_OVER_TWO_R_SQUARED = 1./(2.*R*R);
    const double MAX_CARTESIAN_DISTANCE = 4.*R*R;
    
    for (int i = 0; i < num_points; i++) {
        arr_out[i] = MAX_CARTESIAN_DISTANCE;
        for (int j = 0; j < num_points; j++) {
            if (plate_IDs[j] != plate_IDs[i]) {
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
    return 0;
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
    
    for (int i = i_start; i < i_end; i++) {
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

int full_plate_interior_distance_transform_64bit_threaded(
    double *    xs,
    double *    ys,
    double *    zs,
    int32_t *   plate_IDs,
    int32_t     num_points,
    double      R,
    double *    arr_out,
    int32_t     num_threads
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

    return 0;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     fused_distance_threshold_transform
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/**
 * The following set of functions are to be used when the fused distance 
 * transform is applied on any array of booleans. In all, the original 
 * array `arr` must be different from the output `arr_out`.
 */

int fused_distance_threshold_transform_64bit(
    double *    xs,
    double *    ys,
    double *    zs,
    bool *      arr,
    int32_t     num_points,
    double      R,
    double      threshold,
    bool *      arr_out
) {
    // first, only working with Cartesian distance
    const double THRESHOLD_D_SQUARED = 2.*R*R*(1. - cos(threshold));
    
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


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     gridded_fused_distance_threshold_transform
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/**
 * The following set of functions are to be used when the fused distance 
 * transform is applied on a 2D gridded array of booleans. In all, the original 
 * array `arr` must be different from the output `arr_out`.
 */

int gridded_fused_distance_threshold_transform_64bit(
    double *    xs,
    double *    ys,
    double *    zs,
    bool *      arr,
    int32_t     i_max,
    int32_t     j_max,
    double      R,
    double      threshold,
    bool *      arr_out
) {
    // as a spherical grid, the first and last columns are the same

    // first, only working with Cartesian distance
    const double THRESHOLD_D_SQUARED = 2.*R*R*(1. - cos(threshold));
    
    const double dlat = PI / ((double) (i_max - 1)); 
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

        // the first and last columns are the same
        //arr_out[i*j_max + (j_max-1)] = arr_out[i*j_max];
    }

    return 0;
}

/**
 * (internal)
 * Argument struct for the threaded function:
 * `gridded_fused_distance_threshold_transform_64bit_threaded_func`
 * The parameters are exactly as listed for:
 * `gridded_fused_distance_threshold_transform_64bit_threaded_func`
 */
struct gridded_fused_distance_threshold_transform_64bit_threaded_args {
    double *    xs;
    double *    ys;
    double *    zs;
    bool *      arr;
    int32_t     i_max;
    int32_t     j_max;
    double      R;
    double      threshold;
    bool *      arr_out;
    int32_t     num_threads;
    int32_t     i_thread;
};

/**
 * (internal)
 * Threaded function for `single_plate_interior_distance_transform_64bit_threaded`
 * 
 * @param args the argument struct containing function parameters
 */
void * gridded_fused_distance_threshold_transform_64bit_threaded_func(
    void * args
) {
    /**
     * Assumes that plate_IDs are initialized with -1. and 0. :
     *      -> -1. corresponds to the plate of interest
     *      -> -2. corresponds to other plates
     */

    struct gridded_fused_distance_threshold_transform_64bit_threaded_args * targs =
        (struct gridded_fused_distance_threshold_transform_64bit_threaded_args *) args;
    
    int i_start = (targs->i_max / targs->num_threads)*(targs->i_thread);
    int i_end   = (targs->i_max / targs->num_threads)*(targs->i_thread +1);

    // to prevent missing out due to division round downs
    if (targs->i_thread == targs->num_threads - 1) {
        i_end = targs->i_max;
    }

    // as a spherical grid, the first and last columns are the same

    // first, only working with Cartesian distance
    const double THRESHOLD_D_SQUARED = 2.*targs->R*targs->R*(1. - cos(targs->threshold));
    
    const double dlat = PI / ((double) targs->i_max); 
    const int di = (int) (targs->threshold / dlat) + 2;        // +2 for a conservative

    for (int i = i_start; i < i_end; i++) {
        for (int j = 0; j < targs->j_max; j++) {
            const int ref_index = i*targs->j_max + j;
            if (targs->arr[ref_index]) {
                targs->arr_out[ref_index] = true;

                int ii_start = 0;
                int ii_end   = targs->i_max-1;

                if ((i - di) > 0)               ii_start = i - di;
                if ((i + di) < targs->i_max-1)  ii_end   = i + di;

                for (int ii = ii_start; ii <= ii_end; ii++) {
                    for (int jj = 0; jj < targs->j_max; jj++) {
                        const int index = ii*targs->j_max + jj;
                        if (!targs->arr_out[index]) {
                            // the point [j] is not on the plate

                            // crude checking whether to even bother with
                            // spherical distance
                            const double dx = targs->xs[index] - targs->xs[ref_index];
                            const double dy = targs->ys[index] - targs->ys[ref_index];
                            const double dz = targs->zs[index] - targs->zs[ref_index];
                            const double d2 = dx*dx + dy*dy + dz*dz;
                            
                            targs->arr_out[index] = d2 < THRESHOLD_D_SQUARED;
                        }
                    }
                }
            }    
        }

        // the first and last columns are the same
        //targs->arr_out[i*targs->j_max + (targs->j_max-1)] = 
        //        targs->arr_out[i*targs->j_max];
    }

    return NULL;
}


int gridded_fused_distance_threshold_transform_64bit_threaded(
    double *    xs,
    double *    ys,
    double *    zs,
    bool *      arr,
    int32_t     i_max,
    int32_t     j_max,
    double      R,
    double      threshold,
    bool *      arr_out,
    int32_t     num_threads
) {
    pthread_t tids[num_threads];
    struct gridded_fused_distance_threshold_transform_64bit_threaded_args targss[num_threads];

    for (int i_thread=0; i_thread < num_threads; i_thread++) {
        targss[i_thread].xs             = xs;
        targss[i_thread].ys             = ys;
        targss[i_thread].zs             = zs;
        targss[i_thread].arr            = arr;
        targss[i_thread].i_max          = i_max;
        targss[i_thread].j_max          = j_max;
        targss[i_thread].R              = R;
        targss[i_thread].threshold      = threshold;
        targss[i_thread].arr_out        = arr_out;
        targss[i_thread].num_threads    = num_threads;
        targss[i_thread].i_thread       = i_thread;

        // creating threads
        pthread_create(
            &tids[i_thread], 
            NULL, 
            gridded_fused_distance_threshold_transform_64bit_threaded_func, 
            &targss[i_thread]
        );
    }
    
    // waiting for all threads to be done
    for (int i_thread=0; i_thread < num_threads; i_thread++) {
        pthread_join(tids[i_thread], NULL);
    }
    return 0;
}