#include "segmentation.h"

/**
 * (internal)
 * 
 * Find the index that stores the largest value.
 * 
 * @param a array of values
 * @param i_max size of the array
 * 
 * @returns index
 */
int argmax(
    double * a,
    int i_max
) {
    int index = 0;
    double maxval = a[0];

    for (int i = 1; i < i_max; i++) {
        if (a[i] > maxval) {
            maxval = a[i];
            index = i;
        }
    }

    return index;
}

int populate_edges(
    int32_t *   edges,
    int32_t     n_i,
    int32_t     n_j
) {
    // num_edges must be equal to (n_i-1)*n_j + n_i*(n_j-1)
    // edges is a flattened 2D array with shape [num_edges, 2]
    
    int l = 0;

    for (int i = 0; i < n_i; i++) {
        for (int j = 0; j < n_j-1; j++) {
            edges[l] = i*n_j + j;
            l++;
            edges[l] = i*n_j + (j+1);
            l++;
        }
    }

    for (int i = 0; i < n_i-1; i++) {
        for (int j = 0; j < n_j; j++) {
            edges[l] = i*n_j + j;
            l++;
            edges[l] = (i+1)*n_j + j;
            l++;
        }
    }

    return 0;
}

int populate_edges_sph(
    int32_t *   edges,
    int32_t     n_i,
    int32_t     n_j
) {
    // num_edges must be equal to (2*n_i-3)*n_j
    // edges is a flattened 2D array with shape [num_edges, 2]
    
    int l = 0;

    // north pole
    for (int j = 0; j < n_j; j++) {
        edges[l] = n_j-1;
        l++;
        edges[l] = n_j + j;
        l++;
    }
    // interior horizontal
    for (int i = 1; i < n_i-1; i++) {
        for (int j = 0; j < n_j-1; j++) {
            edges[l] = i*n_j + j;
            l++;
            edges[l] = i*n_j + (j+1);
            l++;
        }
    }
    // interior vertical
    for (int i = 1; i < n_i-2; i++) {
        for (int j = 0; j < n_j; j++) {
            edges[l] = i*n_j + j;
            l++;
            edges[l] = (i+1)*n_j + j;
            l++;
        }
    }
    // interior horizontal wraparound
    for (int i = 1; i < n_i-1; i++) {
        edges[l] = i*n_j;
        l++;
        edges[l] = i*n_j + (n_j-1);
        l++;
    }
    // south pole
    for (int j = 0; j < n_j; j++) {
        edges[l] = (n_i-2)*n_j + j;
        l++;
        edges[l] = (n_i-1)*n_j;
        l++;
    }

    return 0;
}


int get_Laplacian_from_edges(
    double *    image,
    int32_t *   edges,
    int32_t     n_i,
    int32_t     n_j,
    double      beta,
    int32_t *   rows,
    int32_t *   columns,
    double *    values
) {
    // num_edges must be equal to (n_i-1)*n_j + n_i*(n_j-1)
    // image is a flattened 2D array with shape [n_i, n_j]
    // edges is a flattened 2D array with shape [num_edges, 2]
    // rows, columns, and array should be num_edges*2 long 
    int num_edges = (n_i-1)*n_j + n_i*(n_j-1);
    const int edges_size = 2 * num_edges;

    for (int l = 0; l < edges_size-1; l = l + 2){
        const int32_t vertex1 = edges[l];
        const int32_t vertex2 = edges[l+1];
        const double contrast = image[vertex1] - image[vertex2];
        const double value = -exp(-beta*(contrast*contrast));

        rows[l]         = vertex1;
        rows[l+1]       = vertex2;
        columns[l]      = vertex2;
        columns[l+1]    = vertex1;
        values[l]       = value;
        values[l+1]     = value;
    }

    return 0;
}


int get_Laplacian_from_edges_psph(
    double *    image,
    int32_t *   edges,
    int32_t     n_i,
    int32_t     n_j,
    double      theta_min,
    double      theta_max,
    double      phi_min,
    double      phi_max,
    double      beta,
    int32_t *   rows,
    int32_t *   columns,
    double *    values
) {
    // num_edges must be equal to (n_i-1)*n_j + n_i*(n_j-1)
    // image is a flattened 2D array with shape [n_i, n_j]
    // edges is a flattened 2D array with shape [num_edges, 2]
    // rows, columns, and array should be num_edges*2 long 
    int num_edges = (n_i-1)*n_j + n_i*(n_j-1);
    const int edges_size = 2 * num_edges;

    const double d_theta    = (theta_max - theta_min) / ((double) (n_i - 1));
    const double d_phi      = (phi_max - phi_min) / ((double) n_j);
    
    for (int l = 0; l < edges_size-1; l = l + 2){
        const int32_t vertex1 = edges[l];
        const int32_t vertex2 = edges[l+1];
        const double contrast = image[vertex1] - image[vertex2];

        const int i1 = vertex1/n_j;
        const int j1 = vertex1%n_j;
        const int i2 = vertex2/n_j;
        const int j2 = vertex2%n_j;

        double metric_correction = 1.;
        if (i1 == i2) {
            // same latitude and not a polar point
            const double theta = theta_min + ((double) i1) * d_theta;
            const double distance =  sin(theta)*d_phi;
            metric_correction = 1./distance;
        } else if (j1 == j2) {
            // same longitude
            const double distance = d_theta;
            metric_correction = 1./distance;
        }
        if (((isinf(metric_correction)) || (isnan(metric_correction))) 
            || (metric_correction > BIG)) {
            metric_correction = BIG;
        }

        const double value = -exp(-beta*contrast*contrast)*metric_correction;

        rows[l]         = vertex1;
        rows[l+1]       = vertex2;
        columns[l]      = vertex2;
        columns[l+1]    = vertex1;
        values[l]       = value;
        values[l+1]     = value;
    }
    

    return 0;
}


int get_Laplacian_from_edges_sph(
    double *    image,
    int32_t *   edges,
    int32_t     n_i,
    int32_t     n_j,
    double      beta,
    //double      max_beta_spatial_correction,
    int32_t *   rows,
    int32_t *   columns,
    double *    values
) {
    // num_edges must be equal to (2*n_i-3)*n_j
    // image is a flattened 2D array with shape [n_i, n_j]
    // edges is a flattened 2D array with shape [num_edges, 2]
    // rows, columns, and array should be num_edges*2 long 
    int num_edges = (2*n_i-3)*n_j;
    const int edges_size = 2 * num_edges;

    const double d_theta    = PI / ((double) (n_i - 1));
    const double d_phi      = 2. * PI / ((double) n_j);
    
    for (int l = 0; l < edges_size-1; l = l + 2){
        const int32_t vertex1 = edges[l];
        const int32_t vertex2 = edges[l+1];
        const double contrast = image[vertex1] - image[vertex2];

        const int i1 = vertex1/n_j;
        const int j1 = vertex1%n_j;
        const int i2 = vertex2/n_j;
        const int j2 = vertex2%n_j;

        double metric_correction = 1.;
        if (i1 == i2) {
            // same latitude and not a polar point
            const double theta = ((double) i1) * d_theta;
            const double distance =  sin(theta)*d_phi;
            metric_correction = 1./distance;
        } else if (j1 == j2) {
            // same longitude
            const double distance = d_theta;
            metric_correction = 1./distance;
        }
        if (((isinf(metric_correction)) || (isnan(metric_correction))) 
            || (metric_correction > BIG)) {
            metric_correction = BIG;
        }

        const double value = -exp(-beta*contrast*contrast)*metric_correction;

        rows[l]         = vertex1;
        rows[l+1]       = vertex2;
        columns[l]      = vertex2;
        columns[l+1]    = vertex1;
        values[l]       = value;
        values[l+1]     = value;
    }
    

    return 0;
}


int get_original_to_ordered_mapping(
    int32_t *   labels,
    int32_t     n_i,
    int32_t     n_j,
    int32_t *   org2ord,
    int32_t *   ord2org
) {
    int k = 0;
    for (int i = 0; i < n_i; i++) {
        for (int j = 0; j < n_j; j++) {
            if (labels[i*n_j + j] > 0) {
                ord2org[k] = n_j*i + j;
                org2ord[n_j*i + j] = k;
                k++;
            }
        }
    }

    for (int i = 0; i < n_i; i++) {
        for (int j = 0; j < n_j; j++) {
            if (labels[i*n_j + j] == 0) {
                ord2org[k] = n_j*i + j   ; 
                org2ord[n_j*i + j] = k;
                k++;
            }
        }
    }

    return 0;
}

int get_original_to_ordered_mapping_sph(
    int32_t *   labels,
    int32_t     n_i,
    int32_t     n_j,
    int32_t *   org2ord,
    int32_t *   ord2org
) {
    int k = 0;
    // first the labelled ones
    const int shift = n_j-1;
    
    // north pole
    if (labels[n_j-1] > 0) {
        ord2org[k] = n_j-1 - shift;
        org2ord[n_j-1 - shift] = k;
        k++;
    }
    // interior
    for (int i = 1; i < n_i-1; i++) {
        for (int j = 0; j < n_j; j++) {
            if (labels[i*n_j + j] > 0) {
                ord2org[k] = n_j*i + j - shift;
                org2ord[n_j*i + j - shift] = k;
                k++;
            }
        }
    }
    // south pole
    if (labels[(n_i-1)*n_j] > 0) {
        ord2org[k] = (n_i-1)*n_j - shift;
        org2ord[(n_i-1)*n_j - shift] = k;
        k++;
    }

    // then the unlabelled
    // north pole
    if (labels[n_j-1] == 0) {
        ord2org[k] = n_j-1 - shift;
        org2ord[n_j-1 - shift] = k;
        k++;
    }
    // interior
    for (int i = 1; i < n_i-1; i++) {
        for (int j = 0; j < n_j; j++) {
            if (labels[i*n_j + j] == 0) {
                ord2org[k] = n_j*i + j - shift;
                org2ord[n_j*i + j - shift] = k;
                k++;
            }
        }
    }
    // south pole
    if (labels[(n_i-1)*n_j] == 0) {
        ord2org[k] = (n_i-1)*n_j - shift;
        org2ord[(n_i-1)*n_j - shift] = k;
        k++;
    }

    return 0;
}


int order_Laplacian(
    int32_t *   org2ord,
    int32_t *   rows,
    int32_t *   columns,
    int32_t     n_i,
    int32_t     n_j,
    int32_t *   rows_ord,
    int32_t *   columns_ord
) {
    int num_edges = (n_i-1)*n_j + n_i*(n_j-1);
    const int edges_size = 2 * num_edges;

    for (int i = 0; i < edges_size; i++) {
        rows_ord[i] = org2ord[rows[i]];
        columns_ord[i] = org2ord[columns[i]];
    }

    return 0;
}

int order_Laplacian_sph(
    int32_t *   org2ord,
    int32_t *   rows,
    int32_t *   columns,
    int32_t     n_i,
    int32_t     n_j,
    int32_t *   rows_ord,
    int32_t *   columns_ord
) {
    int num_edges = (2*n_i-3)*n_j;
    const int edges_size = 2 * num_edges;

    const int shift = n_j - 1;

    for (int i = 0; i < edges_size; i++) {
        rows_ord[i] = org2ord[rows[i] - shift];
        columns_ord[i] = org2ord[columns[i] - shift];
    }

    return 0;
}


int get_ordered_Laplacian_vectors(
    double *    data,
    int32_t *   labels,
    int32_t     n_i,
    int32_t     n_j,
    double      beta,
    int32_t *   rows_ord,
    int32_t *   columns_ord,
    double *    values,
    int32_t *   ord2org
) {
    // total size of the problem
    const int N = n_i * n_j;

    // number of connecting edges
    int num_edges = (n_i - 1)*n_j + n_i*(n_j - 1);
    
    // the vector of edges containing pairs of [from, to]
    int32_t * edges   = (int32_t *) malloc((2*num_edges) * sizeof(int32_t));

    // vectors to contain non-diagonal Laplacian entries
    int32_t * rows    = (int32_t *) malloc((2*num_edges) * sizeof(int32_t));
    int32_t * columns = (int32_t *) malloc((2*num_edges) * sizeof(int32_t));
    
    // reference bookkeeping vectors to go back and forth between original 
    // (row-major) and ordered (labeled followed by unlabeled).
    int32_t * org2ord = malloc(N * sizeof(int32_t));
    
    populate_edges(edges, n_i, n_j);
    get_Laplacian_from_edges(data, edges, n_i, n_j, beta, rows, columns, values);
    free(edges);
    get_original_to_ordered_mapping(labels, n_i, n_j, org2ord, ord2org);
    order_Laplacian(org2ord, rows, columns, n_i, n_j, rows_ord, columns_ord); 
    free(org2ord);
    free(rows);
    free(columns);

    return 0;
}


int get_ordered_Laplacian_vectors_psph(
    double *    data,
    int32_t *   labels,
    int32_t     n_i,
    int32_t     n_j,
    double      theta_min,
    double      theta_max,
    double      phi_min,
    double      phi_max,
    double      beta,
    int32_t *   rows_ord,
    int32_t *   columns_ord,
    double *    values,
    int32_t *   ord2org
) {
    // total size of the problem
    const int N = n_i * n_j;

    // number of connecting edges
    int num_edges = (n_i - 1)*n_j + n_i*(n_j - 1);
    
    // the vector of edges containing pairs of [from, to]
    int32_t * edges   = (int32_t *) malloc((2*num_edges) * sizeof(int32_t));

    // vectors to contain non-diagonal Laplacian entries
    int32_t * rows    = (int32_t *) malloc((2*num_edges) * sizeof(int32_t));
    int32_t * columns = (int32_t *) malloc((2*num_edges) * sizeof(int32_t));
    
    // reference bookkeeping vectors to go back and forth between original 
    // (row-major) and ordered (labeled followed by unlabeled).
    int32_t * org2ord = malloc(N * sizeof(int32_t));
    
    populate_edges(edges, n_i, n_j);
    get_Laplacian_from_edges_psph(
        data, edges, n_i, n_j, 
        theta_min, theta_max, phi_min, phi_max, 
        beta, rows, columns, values
    );
    free(edges);
    get_original_to_ordered_mapping(labels, n_i, n_j, org2ord, ord2org);
    order_Laplacian(org2ord, rows, columns, n_i, n_j, rows_ord, columns_ord); 
    free(org2ord);
    free(rows);
    free(columns);

    return 0;
}



int get_ordered_Laplacian_vectors_sph(
    double *    data,
    int32_t *   labels,
    int32_t     n_i,
    int32_t     n_j,
    double      beta,
    int32_t *   rows_ord,
    int32_t *   columns_ord,
    double *    values,
    int32_t *   ord2org
) {

    // total size of the problem
    const int N = (n_i-2) * n_j + 2;

    // number of connecting edges
    int num_edges = (2*n_i - 3)*n_j;
    
    // the vector of edges containing pairs of [from, to]
    int32_t * edges   = (int32_t *) malloc((2*num_edges) * sizeof(int32_t));

    // vectors to contain non-diagonal Laplacian entries
    int32_t * rows    = (int32_t *) malloc((2*num_edges) * sizeof(int32_t));
    int32_t * columns = (int32_t *) malloc((2*num_edges) * sizeof(int32_t));
    
    // reference bookkeeping vectors to go back and forth between original 
    // (row-major) and ordered (labeled followed by unlabeled).
    int32_t * org2ord = malloc(N * sizeof(int32_t));
    
    
    populate_edges_sph(edges, n_i, n_j);
    get_Laplacian_from_edges_sph(data, edges, n_i, n_j, beta, rows, columns, values);
    free(edges);
    get_original_to_ordered_mapping_sph(labels, n_i, n_j, org2ord, ord2org);
    order_Laplacian_sph(org2ord, rows, columns, n_i, n_j, rows_ord, columns_ord);
    free(org2ord);
    free(rows);
    free(columns);

    return 0;
}



int get_ordered_boundary_matrix(
    int32_t *   labels,
    int32_t *   ord2org,
    int32_t     n_i,
    int32_t     n_j,
    int32_t     num_labelled,
    int32_t     largest_label,
    double *    M
) {
    for (int k_ord = 0; k_ord < num_labelled; k_ord++){
        for (int s = 0; s < largest_label; s++) {
            const int k_org = ord2org[k_ord];
            const int i_org = k_org/n_j;
            const int j_org = k_org%n_j;
            if (labels[i_org*n_j + j_org] == s+1){
                M[k_ord*largest_label + s] = 1.;
            }
        }
    }
    return 0;
}

int get_ordered_boundary_matrix_sph(
    int32_t *   labels,
    int32_t *   ord2org,
    int32_t     n_i,
    int32_t     n_j,
    int32_t     num_labelled,
    int32_t     largest_label,
    double *    M
) {
    // considering the shift in indicies since the first node considered is the 
    // north pole at the last column of the first row.
    const int shift = n_j - 1;

    for (int k_ord = 0; k_ord < num_labelled; k_ord++){
        for (int s = 0; s < largest_label; s++) {
            const int k_org = ord2org[k_ord] + shift;
            const int i_org = k_org/n_j;
            const int j_org = k_org%n_j;
            if (labels[i_org*n_j + j_org] == s+1){
                M[k_ord*largest_label + s] = 1.;
            }
        }
    }
    return 0;
}

int get_IDs_and_probs_from_X(
    double *    X,
    int32_t     X_i_max,
    int32_t     X_j_max,
    int32_t *   labels,
    int32_t *   ord2org,
    int32_t     n_i,
    int32_t     n_j,
    int32_t     num_labelled,
    int32_t *   IDs,
    double *    probs
) {
    const int N = n_i * n_j;

    for (int k_ord = 0; k_ord < N; k_ord++) {
        const int k_org = ord2org[k_ord];
        if (k_ord < num_labelled) {
            IDs[k_org] = labels[k_org];
            probs[k_org*X_j_max + labels[k_org]-1] = 1.;
        } else {
            IDs[k_org] = argmax(
                X + (k_ord-num_labelled)*X_j_max,
                X_j_max
            ) + 1;
            for (int j = 0; j < X_j_max; j++) {
                probs[k_org*X_j_max + j] = X[(k_ord-num_labelled)*X_j_max + j];
            }
        }
    }


    return 0;
}


int get_IDs_and_probs_from_X_sph(
    double *    X,
    int32_t     X_i_max,
    int32_t     X_j_max,
    int32_t *   labels,
    int32_t *   ord2org,
    int32_t     n_i,
    int32_t     n_j,
    int32_t     num_labelled,
    int32_t *   IDs,
    double *    probs
) {
    const int N = (n_i-2) * n_j + 2;
    const int shift = n_j - 1;

    for (int k_ord = 0; k_ord < N; k_ord++) {
        const int k_org = ord2org[k_ord] + shift;
        if (k_ord < num_labelled) {
            IDs[k_org] = labels[k_org];
            probs[k_org*X_j_max + labels[k_org]-1] = 1.;
        } else {
            IDs[k_org] = argmax(
                X + (k_ord-num_labelled)*X_j_max,
                X_j_max
            ) + 1;
            for (int j = 0; j < X_j_max; j++) {
                probs[k_org*X_j_max + j] = X[(k_ord-num_labelled)*X_j_max + j];
            }
        }
    }

    return 0;
}

int get_segment_boundaries(
    int32_t *   IDs,
    int32_t     n_i,
    int32_t     n_j,
    bool        is_wraparound,
    bool *      boundaries
){
    int i;
    int j;
    for (i = 1; i < n_i-1; i++) {
        for (j = 1; j < n_j-1; j++) {
            //const int ref_ID = IDs;
            const int sum = IDs[(i-1)*n_j + (j-1)] + IDs[(i-1)*n_j + (j  )] + IDs[(i-1)*n_j + (j+1)]
                          + IDs[(i  )*n_j + (j-1)] + IDs[(i  )*n_j + (j  )] + IDs[(i  )*n_j + (j+1)]
                          + IDs[(i+1)*n_j + (j-1)] + IDs[(i+1)*n_j + (j  )] + IDs[(i+1)*n_j + (j+1)];
            boundaries[i*n_j + j] = (sum != 9*IDs[i*n_j + j]);
        }
    }
    i = 1;
    for (j = 1; j < n_j-1; j++) {
        //const int ref_ID = IDs;
        const int sum = IDs[(i  )*n_j + (j-1)] + IDs[(i  )*n_j + (j  )] + IDs[(i  )*n_j + (j+1)]
                      + IDs[(i+1)*n_j + (j-1)] + IDs[(i+1)*n_j + (j  )] + IDs[(i+1)*n_j + (j+1)];
        boundaries[i*n_j + j] = (sum != 6*IDs[i*n_j + j]);
    }
    i = n_i-1;
    for (j = 1; j < n_j-1; j++) {
        //const int ref_ID = IDs;
        const int sum = IDs[(i-1)*n_j + (j-1)] + IDs[(i-1)*n_j + (j  )] + IDs[(i-1)*n_j + (j+1)]
                      + IDs[(i  )*n_j + (j-1)] + IDs[(i  )*n_j + (j  )] + IDs[(i  )*n_j + (j+1)];
        boundaries[i*n_j + j] = (sum != 6*IDs[i*n_j + j]);
    }
    if (is_wraparound) {
        int sum;
        i = 0;
        j = 0;
        sum = IDs[(i  )*n_j + (j  )] + IDs[(i  )*n_j + (j+1)]
            + IDs[(i+1)*n_j + (j  )] + IDs[(i+1)*n_j + (j+1)];
        boundaries[i*n_j + j] = (sum != 4*IDs[i*n_j + j]);

        i = 0;
        j = n_j-1;
        sum = IDs[(i  )*n_j + (j  )] + IDs[(i  )*n_j + (j-1)]
            + IDs[(i+1)*n_j + (j  )] + IDs[(i+1)*n_j + (j-1)];
        boundaries[i*n_j + j] = (sum != 4*IDs[i*n_j + j]);

        i = n_i-1;
        j = 0;
        sum = IDs[(i  )*n_j + (j  )] + IDs[(i  )*n_j + (j+1)]
            + IDs[(i-1)*n_j + (j  )] + IDs[(i-1)*n_j + (j+1)];
        boundaries[i*n_j + j] = (sum != 4*IDs[i*n_j + j]);

        i = n_i-1;
        j = n_j-1;
        sum = IDs[(i  )*n_j + (j  )] + IDs[(i  )*n_j + (j-1)]
            + IDs[(i-1)*n_j + (j  )] + IDs[(i-1)*n_j + (j-1)];
        boundaries[i*n_j + j] = (sum != 4*IDs[i*n_j + j]);

        j = 0;
        for (i = 1; i < n_i-1; i++) {
            //const int ref_ID = IDs;
            const int sum = IDs[(i-1)*n_j + (n_j-1)] + IDs[(i-1)*n_j + (j  )] + IDs[(i-1)*n_j + (j+1)]
                          + IDs[(i  )*n_j + (n_j-1)] + IDs[(i  )*n_j + (j  )] + IDs[(i  )*n_j + (j+1)]
                          + IDs[(i+1)*n_j + (n_j-1)] + IDs[(i+1)*n_j + (j  )] + IDs[(i+1)*n_j + (j+1)];
            boundaries[i*n_j + j] = (sum != 9*IDs[i*n_j + j]);
        }

        j = n_j-1;
        for (i = 1; i < n_i-1; i++) {
            //const int ref_ID = IDs;
            const int sum = IDs[(i-1)*n_j + (j-1)] + IDs[(i-1)*n_j + (j  )] + IDs[(i-1)*n_j + (0  )]
                          + IDs[(i  )*n_j + (j-1)] + IDs[(i  )*n_j + (j  )] + IDs[(i  )*n_j + (0  )]
                          + IDs[(i+1)*n_j + (j-1)] + IDs[(i+1)*n_j + (j  )] + IDs[(i+1)*n_j + (0  )];
            boundaries[i*n_j + j] = (sum != 9*IDs[i*n_j + j]);
        }
    } else {
        int sum;
        i = 0;
        j = 0;
        sum = IDs[(i  )*n_j + (j  )] + IDs[(i  )*n_j + (j+1)]
            + IDs[(i+1)*n_j + (j  )] + IDs[(i+1)*n_j + (j+1)]
            + IDs[(i  )*n_j + (n_j-1)] + IDs[(i+1)*n_j + (n_j-1)];
        boundaries[i*n_j + j] = (sum != 6*IDs[i*n_j + j]);

        i = 0;
        j = n_j-1;
        sum = IDs[(i  )*n_j + (j  )] + IDs[(i  )*n_j + (j-1)]
            + IDs[(i+1)*n_j + (j  )] + IDs[(i+1)*n_j + (j-1)]
            + IDs[(i  )*n_j + (0  )] + IDs[(i+1)*n_j + (0  )];
        boundaries[i*n_j + j] = (sum != 6*IDs[i*n_j + j]);

        i = n_i-1;
        j = 0;
        sum = IDs[(i  )*n_j + (j  )] + IDs[(i  )*n_j + (j+1)]
            + IDs[(i-1)*n_j + (j  )] + IDs[(i-1)*n_j + (j+1)]
            + IDs[(i  )*n_j + (n_j-1)] + IDs[(i-1)*n_j + (n_j-1)];
        boundaries[i*n_j + j] = (sum != 6*IDs[i*n_j + j]);

        i = n_i-1;
        j = n_j-1;
        sum = IDs[(i  )*n_j + (j  )] + IDs[(i  )*n_j + (j-1)]
            + IDs[(i-1)*n_j + (j  )] + IDs[(i-1)*n_j + (j-1)]
            + IDs[(i  )*n_j + (0  )] + IDs[(i-1)*n_j + (0  )];
        boundaries[i*n_j + j] = (sum != 6*IDs[i*n_j + j]);

        j = 0;
        for (i = 1; i < n_i-1; i++) {
            //const int ref_ID = IDs;
            const int sum = IDs[(i-1)*n_j + (j  )] + IDs[(i-1)*n_j + (j+1)]
                          + IDs[(i  )*n_j + (j  )] + IDs[(i  )*n_j + (j+1)]
                          + IDs[(i+1)*n_j + (j  )] + IDs[(i+1)*n_j + (j+1)]
                          + IDs[(i+1)*n_j + (j  )] + IDs[(i+1)*n_j + (j+1)];
            boundaries[i*n_j + j] = (sum != 9*IDs[i*n_j + j]);
        }

        j = n_j-1;
        for (i = 1; i < n_i-1; i++) {
            //const int ref_ID = IDs;
            const int sum = IDs[(i-1)*n_j + (j-1)] + IDs[(i-1)*n_j + (j  )] 
                          + IDs[(i  )*n_j + (j-1)] + IDs[(i  )*n_j + (j  )] 
                          + IDs[(i+1)*n_j + (j-1)] + IDs[(i+1)*n_j + (j  )];
            boundaries[i*n_j + j] = (sum != 9*IDs[i*n_j + j]);
        }
        
    }
    return 0;
}