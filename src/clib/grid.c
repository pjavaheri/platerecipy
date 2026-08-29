/**
 * @file grid.c
 * @author Pejvak Javaheri; pejvak.javaheri@mail.utoronto.ca
 * @brief Definitions for grid structures
 */

#include "grid.h"

#ifdef _WIN32
    #ifndef PYINIT_LIBPLATERECIPY_GRID
        #define PYINIT_LIBPLATERECIPY_GRID
        void PyInit_libplaterecipy_grid() {};
    #endif
#endif

// ~~~~~ Callback mechanism for the Python logger ~~~~~

void set_grid_h_logger(
    grid_h_log_func func
) {
    grid_h_logger = func;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




// ~~~~~~~~~~~~~~~~ Internal functions ~~~~~~~~~~~~~~~~

/**
 * (Internal)
 * Calculate and set num_edges by going through all the nodes.
 *
 * @warning This assumes num_neighs for all the nodes are up-to-date.
 */
void set_num_edges(
    Map * map
) {
    map->num_edges = 0;
    for (int i = 0; i < map->num_nodes; i++) {
        map->num_edges += map->nodes[i].num_neighs;
    }
    map->num_edges /= 2;
}

/**
 * (Internal)
 * Allocated a shared neighs array of size neighs_array_length,
 * and assign sub-arrays based on each node's num_neighs.
 *
 * @warning This assumes num_neighs for all the nodes are up-to-date.
 * @warning The allocation is done using malloc (i.e., not zeroed).
 */
void alloc_and_assign_neighs(
    Map * map,
    int neighs_array_length
) {
    // a single array will be shared by all nodes
    map->nodes[0].neighs = (int32_t *) malloc(neighs_array_length * sizeof(int32_t));
    
    // we set the neighs array for each node to point at a 
    // sub-array of the larger shared array
    for (int i = 1; i < map->num_nodes; i++) {
        map->nodes[i].neighs = map->nodes[i-1].neighs + map->nodes[i-1].num_neighs;
    } 
}

/**
 * (Internal)
 * Allocated a shared edge_lengths array of size num_edges,
 * and assign sub-arrays based on each node's num_neighs.
 *
 * @warning This assumes num_edges is up-to-date.
 * @warning The allocation is done using malloc (i.e., not zeroed).
 */
void alloc_and_init_edge_lengths(
    Map * map
) {
    // a single array is shared by all the nodes
    map->nodes[0].edge_lengths = (double *) malloc(2*map->num_edges * sizeof(double));
    for (int i = 1; i < map->num_nodes; i++) {
        map->nodes[i].edge_lengths = map->nodes[i-1].edge_lengths + map->nodes[i-1].num_neighs;
        for (int j = 0; j < map->nodes[i].num_neighs; j++) {
            map->nodes[i].edge_lengths[j] = 1.;
        }
    }
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




// ~~~~~~~~~~~~~~~~ Exposed functions ~~~~~~~~~~~~~~~~~

Node * get_node_at(
    Map *   map, 
    int32_t i
) {
    if ((i >= 0) && (i < map->num_nodes)) {
        return &(map->nodes[i]);
    }
    if (grid_h_logger) grid_h_logger("ERROR: Requesting for a node at an illegal index.");
    return NULL;
}


int32_t get_neigh_at(
    Node *  nodes, 
    int32_t i
) {
    if ((i >= 0) && (i < nodes->num_neighs)) {
        return nodes->neighs[i];
    }
    if (grid_h_logger) grid_h_logger("ERROR: Requesting for a neighbour at an illegal index.");
    return -1;
}


void set_npy_idx(
    Map *       map,
    int32_t *   npy_idxs
) {
    for (int i = 0; i < map->num_nodes; i++) {
        map->nodes[i].npy_idx = npy_idxs[i];
    }
}


Map * alloc_map_from_cells(
    int32_t     num_nodes,
    int32_t *   cells,
    int32_t     num_cells
) {
    if (grid_h_logger) grid_h_logger("Allocating a custom map from cell/node adjacency ...");

    Map * map = (Map *) calloc(1, sizeof(Map));
    map->num_nodes = num_nodes;
    map->i_max = -1;
    map->j_max = -1;
    map->grid_type = GTYPE_CUSTOM;
    map->nodes = (Node *) calloc(map->num_nodes, sizeof(Node)); 

    // getting an estimate on the shared array size
    int estimated_num_edges = 0;
    int c = 0;
    while (c < num_cells) {
        const int n = cells[c];
        c++;
        for (int i = 0; i < n-1; i++) {
            const int idx1 = cells[c+i];
            const int idx2 = cells[c+i+1];
            map->nodes[idx1].num_neighs++;
            map->nodes[idx2].num_neighs++;
            estimated_num_edges++;
        }
        const int idx1 = cells[c];
        const int idx2 = cells[c+n-1];
        map->nodes[idx1].num_neighs++;
        map->nodes[idx2].num_neighs++;
        estimated_num_edges++;
        c += n;
    }

    // a single array will be shared by all nodes
    // we assume a larger array length than we need to
    const int neighs_array_length = 2*estimated_num_edges;
    alloc_and_assign_neighs(map, neighs_array_length);
    // initializing the array by a delimiter, so that we can 
    // later keep track of how many actuall neighboours each point 
    // has
    for (int i = 0; i < neighs_array_length; i++) {
        map->nodes[0].neighs[i] = -1; // array delimiter
    }
    
    // going over the points again to correct the sub-array lengths
    c = 0;
    while (c < num_cells) {
        const int n = cells[c];
        c++;
        for (int i = 0; i < n-1; i++) {
            const int idx1 = cells[c+i];
            const int idx2 = cells[c+i+1];
            for (int j = 0; j < map->nodes[idx1].num_neighs; j++) {
                if (idx2 == map->nodes[idx1].neighs[j]) {
                    // there's a duplicate
                    map->nodes[idx1].num_neighs--;
                    break;
                }
                if (-1 == map->nodes[idx1].neighs[j]) {
                    // this is a new connection
                    map->nodes[idx1].neighs[j] = idx2;
                    break;
                }
            }
            for (int j = 0; j < map->nodes[idx2].num_neighs; j++) {
                if (idx1 == map->nodes[idx2].neighs[j]) {
                    // there's a duplicate
                    map->nodes[idx2].num_neighs--;
                    break;
                }
                if (-1 == map->nodes[idx2].neighs[j]) {
                    // this is a new connection
                    map->nodes[idx2].neighs[j] = idx1;
                    break;
                }
            }
        }
        const int idx1 = cells[c];
        const int idx2 = cells[c+n-1];
        for (int j = 0; j < map->nodes[idx1].num_neighs; j++) {
            if (idx2 == map->nodes[idx1].neighs[j]) {
                // there's a duplicate
                map->nodes[idx1].num_neighs--;
                break;
            }
            if (-1 == map->nodes[idx1].neighs[j]) {
                // this is a new connection
                map->nodes[idx1].neighs[j] = idx2;
                break;
            }
        }
        for (int j = 0; j < map->nodes[idx2].num_neighs; j++) {
            if (idx1 == map->nodes[idx2].neighs[j]) {
                // there's a duplicate
                map->nodes[idx2].num_neighs--;
                break;
            }
            if (-1 == map->nodes[idx2].neighs[j]) {
                // this is a new connection
                map->nodes[idx2].neighs[j] = idx1;
                break;
            }
        }

        c += n;
    }

    // with each map->node[i].num_edges corrected, we can now
    // calculate the number of edges
    set_num_edges(map);

    // with num_edges calculated, we initialize the edge_lengths
    alloc_and_init_edge_lengths(map);

    return map;
}


Map * alloc_rect_map(
    int32_t i_max,
    int32_t j_max,
    int32_t grid_type
) {
    if (grid_h_logger) grid_h_logger("Allocating a basic rectilinear map ...");

    int i, j, index;

    Map * map = (Map *) calloc(1, sizeof(Map));
    map->num_nodes = i_max*j_max;
    map->nodes = (Node *) calloc(1, map->num_nodes * sizeof(Node));
    map->i_max = i_max;
    map->j_max = j_max;
    map->grid_type = grid_type;

    for (i = 0; i < map->num_nodes; i++) {
        map->nodes[i].idx = i;
        map->nodes[i].npy_idx = -1;     // unknown
        map->nodes[i].ord_idx = -1;     // unknown
    }

    
    // first, we set num_neighs for sub-array assignment
    // top left corner
    i = 0;
    j = 0;
    index = 0;
    map->nodes[index].num_neighs = 2;
    
    // top right corner
    i = 0;
    j = j_max-1;
    index = j;
    map->nodes[index].num_neighs = 2;

    // bottom left corner
    i = i_max-1;
    j = 0;
    index = i*j_max;
    map->nodes[index].num_neighs = 2;

    // top left corner
    i = i_max-1;
    j = j_max-1;
    index = map->num_nodes-1;
    map->nodes[index].num_neighs = 2;

    // top side nodes
    i = 0;
    for (j = 1; j < j_max-1; j++) {
        index = i*j_max + j;
        map->nodes[index].num_neighs = 3;
    }

    // left side nodes
    j = 0;
    for (i = 1; i < i_max-1; i++) {
        index = i*j_max + j;
        map->nodes[index].num_neighs = 3;
    }
    
    // right side nodes
    j = j_max-1;
    for (i = 1; i < i_max-1; i++) {
        index = i*j_max + j;
        map->nodes[index].num_neighs = 3;
    }

    // bottom side nodes
    i = i_max-1;
    for (j = 1; j < j_max-1; j++) {
        index = i*j_max + j;
        map->nodes[index].num_neighs = 3;
    }
    

    // interior nodes have 4 connections
    for (i = 1; i < i_max-1; i++) {
        for (j = 1; j < j_max-1; j++) {
            index = i*j_max + j;
            map->nodes[index].num_neighs = 4;
        }
    }

    // assigning sub arrays
    // allocating a large array to be shared by all nodes
    const int neighs_array_length = (i_max-2)*(j_max-2)*4 + 2*(i_max-2)*3 + 2*(j_max-2)*3 + 4*2;
    alloc_and_assign_neighs(map, neighs_array_length);

    // then, we set nieghs each node
    // top left corner
    i = 0;
    j = 0;
    index = 0;
    map->nodes[index].neighs[0] = i*j_max + (j+1);    // right
    map->nodes[index].neighs[1] = (i+1)*j_max + j;    // bottom
    
    // top right corner
    i = 0;
    j = j_max-1;
    index = j;
    map->nodes[index].neighs[0] = i*j_max + (j-1);    // left
    map->nodes[index].neighs[1] = (i+1)*j_max + j;    // bottom

    // bottom left corner
    i = i_max-1;
    j = 0;
    index = i*j_max;
    map->nodes[index].neighs[0] = i*j_max + (j+1);    // right
    map->nodes[index].neighs[1] = (i-1)*j_max + j;    // top

    // top left corner
    i = i_max-1;
    j = j_max-1;
    index = map->num_nodes-1;
    map->nodes[index].neighs[0] = (i-1)*j_max + j;    // top
    map->nodes[index].neighs[1] = i*j_max + (j-1);    // left

    // top side nodes
    i = 0;
    for (j = 1; j < j_max-1; j++) {
        index = i*j_max + j;
        map->nodes[index].neighs[0] = i*j_max + (j+1);    // right
        map->nodes[index].neighs[1] = i*j_max + (j-1);    // left
        map->nodes[index].neighs[2] = (i+1)*j_max + j;    // bottom
    }

    // left side nodes
    j = 0;
    for (i = 1; i < i_max-1; i++) {
        index = i*j_max + j;
        map->nodes[index].neighs[0] = i*j_max + (j+1);    // right
        map->nodes[index].neighs[1] = (i-1)*j_max + j;    // top
        map->nodes[index].neighs[2] = (i+1)*j_max + j;    // bottom
    }
    
    // right side nodes
    j = j_max-1;
    for (i = 1; i < i_max-1; i++) {
        index = i*j_max + j;
        map->nodes[index].neighs[0] = (i-1)*j_max + j;    // top
        map->nodes[index].neighs[1] = i*j_max + (j-1);    // left
        map->nodes[index].neighs[2] = (i+1)*j_max + j;    // bottom
    }

    // bottom side nodes
    i = i_max-1;
    for (j = 1; j < j_max-1; j++) {
        index = i*j_max + j;
        map->nodes[index].neighs[0] = i*j_max + (j+1);    // right
        map->nodes[index].neighs[1] = (i-1)*j_max + j;    // top
        map->nodes[index].neighs[2] = i*j_max + (j-1);    // left
    }
    

    // interior nodes have 4 connections
    for (i = 1; i < i_max-1; i++) {
        for (j = 1; j < j_max-1; j++) {
            index = i*j_max + j;
            map->nodes[index].neighs[0] = i*j_max + (j+1);    // right
            map->nodes[index].neighs[1] = (i-1)*j_max + j;    // top
            map->nodes[index].neighs[2] = i*j_max + (j-1);    // left
            map->nodes[index].neighs[3] = (i+1)*j_max + j;    // bottom
        }
    }


    set_num_edges(map);
    alloc_and_init_edge_lengths(map);
    return map;
}

Map * alloc_sph_map(
    int32_t i_max,
    int32_t j_max
) {
    /**
    Note: the node ordering is different here compared to the numpy array.
    Here, the nodes are ordered by first interior, then south pole, then north pole
    but in numpy, it's first north pole, then interior, then south pole.
     */
    
    if (grid_h_logger) grid_h_logger("Allocating a basic spherical map ...");

    const int n_i = i_max - 2; // exluding the polar rows
    const int n_j = j_max;  // wrap-around column is already trimmed off

    Map * map = (Map *) calloc(1, sizeof(Map));
    map->num_nodes = n_i*n_j + 2;
    map->i_max = n_i;
    map->j_max = n_j;
    map->nodes = (Node *) calloc(map->num_nodes, sizeof(Node));

    map->grid_type = GTYPE_SPHERICAL;

    const int north_pole_idx = map->num_nodes - 1;
    const int south_pole_idx = map->num_nodes - 2;
    
    int i, j, index;

    for (i = 0; i < map->num_nodes; i++) {
        map->nodes[i].idx = i;
        map->nodes[i].npy_idx = -1;     // unknown
        map->nodes[i].ord_idx = -1;     // unknown
    }

    
    // first, setting num_neighs
    // top left corner
    for (int i = 0; i < map->num_nodes - 2; i++) {
        map->nodes[i].num_neighs = 4;
    }

    // north pole
    index = north_pole_idx;
    map->nodes[index].num_neighs = n_j;

    // south pole
    index = south_pole_idx;
    map->nodes[index].num_neighs = n_j;

    // assigning sub arrays
    // sharing the same array 
    const int neighs_array_length = (n_i*n_j)*4 + 2*n_j;
    alloc_and_assign_neighs(map, neighs_array_length);

    // then, setting num_neighs
    // top left corner
    i = 0;
    j = 0;
    index = 0;
    map->nodes[index].neighs[0] = i*n_j + (j+1);      // right
    map->nodes[index].neighs[1] = north_pole_idx;       // top
    map->nodes[index].neighs[2] = i*n_j + (n_j-1);  // left
    map->nodes[index].neighs[3] = (i+1)*n_j + j;      // bottom
    
    // top right corner
    i = 0;
    j = n_j-1;
    index = j;
    map->nodes[index].neighs[0] = i*n_j;              // right
    map->nodes[index].neighs[1] = north_pole_idx;       // top
    map->nodes[index].neighs[2] = i*n_j + (j-1);      // left
    map->nodes[index].neighs[3] = (i+1)*n_j + j;      // bottom

    // bottom left corner
    i = n_i-1;
    j = 0;
    index = i*n_j;
    map->nodes[index].neighs[0] = i*n_j + (j+1);      // right
    map->nodes[index].neighs[1] = (i-1)*n_j + j;      // top
    map->nodes[index].neighs[2] = i*n_j + (n_j-1);  // left
    map->nodes[index].neighs[3] = south_pole_idx;       // bottom


    // bottom right corner
    i = n_i-1;
    j = n_j-1;
    index = i*n_j + j;
    map->nodes[index].neighs[0] = i*n_j;              // right
    map->nodes[index].neighs[1] = (i-1)*n_j + j;      // top
    map->nodes[index].neighs[2] = i*n_j + (j-1);      // left
    map->nodes[index].neighs[3] = south_pole_idx;       // bottom

    // top side nodes
    i = 0;
    for (j = 1; j < n_j-1; j++) {
        index = i*n_j + j;
        map->nodes[index].neighs[0] = i*n_j + (j+1);    // right
        map->nodes[index].neighs[1] = north_pole_idx;     // top
        map->nodes[index].neighs[2] = i*n_j + (j-1);    // left
        map->nodes[index].neighs[3] = (i+1)*n_j + j;    // bottom
    }

    // left side nodes
    j = 0;
    for (i = 1; i < n_i-1; i++) {
        index = i*n_j + j;
        map->nodes[index].neighs[0] = i*n_j + (j+1);    // right
        map->nodes[index].neighs[1] = (i-1)*n_j + j;    // top
        map->nodes[index].neighs[2] = i*n_j + (n_j-1);    // left
        map->nodes[index].neighs[3] = (i+1)*n_j + j;    // bottom
    }
    
    // right side nodes
    j = n_j-1;
    for (i = 1; i < n_i-1; i++) {
        index = i*n_j + j;
        map->nodes[index].neighs[0] = i*n_j;            // right
        map->nodes[index].neighs[1] = (i-1)*n_j + j;    // top
        map->nodes[index].neighs[2] = i*n_j + (j-1);    // left
        map->nodes[index].neighs[3] = (i+1)*n_j + j;    // bottom
    }

    // bottom side nodes
    i = n_i-1;
    for (j = 1; j < n_j-1; j++) {
        index = i*n_j + j;
        map->nodes[index].neighs[0] = i*n_j + (j+1);    // right
        map->nodes[index].neighs[1] = (i-1)*n_j + j;    // top
        map->nodes[index].neighs[2] = i*n_j + (j-1);    // left
        map->nodes[index].neighs[3] = south_pole_idx;     // bottom
    }
    

    // interior nodes have 4 connections
    for (i = 1; i < n_i-1; i++) {
        for (j = 1; j < n_j-1; j++) {
            index = i*n_j + j;
            map->nodes[index].neighs[0] = i*n_j + (j+1);    // right
            map->nodes[index].neighs[1] = (i-1)*n_j + j;    // top
            map->nodes[index].neighs[2] = i*n_j + (j-1);    // left
            map->nodes[index].neighs[3] = (i+1)*n_j + j;    // bottom
        }
    }

    // north pole
    index = north_pole_idx;
    for (int j = 0; j < n_j; j++) {
        map->nodes[index].neighs[j] = j;
    }

    // south pole
    index = south_pole_idx;
    for (int j = 0; j < n_j; j++) {
        map->nodes[index].neighs[j] = (n_i-1)*n_j + j;
    }


    set_num_edges(map);
    alloc_and_init_edge_lengths(map);
    return map;
}

void set_cartesian_edge_length(
    Map *       map,
    double *    xs,
    double *    ys,
    double *    zs,
    double      R
) {
    for (int i = 0; i < map->num_nodes; i++) {
        const double x0 = xs[map->nodes[i].npy_idx];
        const double y0 = ys[map->nodes[i].npy_idx];
        const double z0 = zs[map->nodes[i].npy_idx];

        for (int j = 0; j < map->nodes[i].num_neighs; j++) {
            const double dxj = xs[map->nodes[map->nodes[i].neighs[j]].npy_idx]-x0;
            const double dyj = ys[map->nodes[map->nodes[i].neighs[j]].npy_idx]-y0;
            const double dzj = zs[map->nodes[map->nodes[i].neighs[j]].npy_idx]-z0;
            map->nodes[i].edge_lengths[j] = sqrt(dxj*dxj + dyj*dyj + dzj*dzj)/R;
        }
    }
}

void set_spherical_edge_length(
    Map *       map,
    double *    xs,
    double *    ys,
    double *    zs,
    double      R
) {
    const double R2 = R*R;
    for (int i = 0; i < map->num_nodes; i++) {
        const double x0 = xs[map->nodes[i].npy_idx];
        const double y0 = ys[map->nodes[i].npy_idx];
        const double z0 = zs[map->nodes[i].npy_idx];

        for (int j = 0; j < map->nodes[i].num_neighs; j++) {
            const double xj = xs[map->nodes[map->nodes[i].neighs[j]].npy_idx];
            const double yj = ys[map->nodes[map->nodes[i].neighs[j]].npy_idx];
            const double zj = zs[map->nodes[map->nodes[i].neighs[j]].npy_idx];
            map->nodes[i].edge_lengths[j] = acos(((xj*x0) + (yj*y0) + (zj*z0))/R2);
        }
    }
}

void free_map(Map * map) {
    if (grid_h_logger) grid_h_logger("Deallocating map memory ...");

    // both arrays are actually a single 
    // large shared array
    free(map->nodes[0].neighs);
    free(map->nodes[0].edge_lengths);
    
    free(map->nodes);
    free(map);
}


void gen_cells_from_map(
    Map *       map,
    int32_t *   cells
) {
    if (grid_h_logger) grid_h_logger("Generating cells from map ...");

    if ((map->grid_type == GTYPE_PLANAR) || (map->grid_type == GTYPE_PARTIAL_SPHERICAL)) {
        if (grid_h_logger) grid_h_logger("... it is a rect map");
        int c = 0;

        for (int i = 0; i < map->i_max-1; i++) {
            for (int j = 0; j < map->j_max-1; j++) {
                const int idx1 = i*map->j_max + j;
                const int idx2 = (i+1)*map->j_max + j;
                const int idx3 = (i+1)*map->j_max + (j+1);
                const int idx4 = i*map->j_max + (j+1);
                cells[c++] = 4;
                cells[c++] = map->nodes[idx1].npy_idx;
                cells[c++] = map->nodes[idx2].npy_idx;
                cells[c++] = map->nodes[idx3].npy_idx;
                cells[c++] = map->nodes[idx4].npy_idx;
            }
        }

    } else if (map->grid_type == GTYPE_SPHERICAL) {
        // ensuring the first point (north pole) on the numpy array 
        // has index 0
        const int idx_npy_offset = map->nodes[map->num_nodes-1].npy_idx;

        int c = 0;

        for (int i = 0; i < map->i_max-1; i++) {
            for (int j = 0; j < map->j_max-1; j++) {
                const int idx1 = i*map->j_max + j;
                const int idx2 = (i+1)*map->j_max + j;
                const int idx3 = (i+1)*map->j_max + (j+1);
                const int idx4 = i*map->j_max + (j+1);
                cells[c++] = 4;
                cells[c++] = map->nodes[idx1].npy_idx - idx_npy_offset;
                cells[c++] = map->nodes[idx2].npy_idx - idx_npy_offset;
                cells[c++] = map->nodes[idx3].npy_idx - idx_npy_offset;
                cells[c++] = map->nodes[idx4].npy_idx - idx_npy_offset;
            }
            // wrap-around cells
            const int idx1 = i*map->j_max + (map->j_max-1);
            const int idx2 = (i+1)*map->j_max + (map->j_max-1);
            const int idx3 = (i+1)*map->j_max + 0;
            const int idx4 = i*map->j_max + 0;
            cells[c++] = 4;
            cells[c++] = map->nodes[idx1].npy_idx - idx_npy_offset;
            cells[c++] = map->nodes[idx2].npy_idx - idx_npy_offset;
            cells[c++] = map->nodes[idx3].npy_idx - idx_npy_offset;
            cells[c++] = map->nodes[idx4].npy_idx - idx_npy_offset;
        }

        // south pole
        const int south_pole_idx = map->num_nodes - 2;
        for (int j = 0; j < map->j_max-1; j++) {
            const int idx1 = (map->i_max-1)*map->j_max + j;
            const int idx3 = (map->i_max-1)*map->j_max + j+1;
            cells[c++] = 3;
            cells[c++] = map->nodes[idx1].npy_idx - idx_npy_offset;
            cells[c++] = map->nodes[south_pole_idx].npy_idx - idx_npy_offset;
            cells[c++] = map->nodes[idx3].npy_idx - idx_npy_offset;
        }
        cells[c++] = 3;
        cells[c++] = map->nodes[(map->i_max-1)*map->j_max + (map->j_max-1)].npy_idx - idx_npy_offset;
        cells[c++] = map->nodes[south_pole_idx].npy_idx - idx_npy_offset;
        cells[c++] = map->nodes[(map->i_max-1)*map->j_max].npy_idx - idx_npy_offset;

        // north pole
        const int north_pole_idx = map->num_nodes - 1;
        for (int j = 0; j < map->j_max-1; j++) {
            const int idx2 = j;
            const int idx3 = j+1;
            cells[c++] = 3;
            cells[c++] = map->nodes[north_pole_idx].npy_idx - idx_npy_offset;
            cells[c++] = map->nodes[idx2].npy_idx - idx_npy_offset;
            cells[c++] = map->nodes[idx3].npy_idx - idx_npy_offset;
        }
        cells[c++] = 3;
        cells[c++] = map->nodes[north_pole_idx].npy_idx - idx_npy_offset;
        cells[c++] = map->nodes[map->j_max-1].npy_idx - idx_npy_offset;
        cells[c++] = map->nodes[0].npy_idx - idx_npy_offset;

    
    } else {
        // unknown grid 
        if (grid_h_logger) grid_h_logger("ERROR: unkown map.");
        cells[0] = -1;
    }
    if (grid_h_logger) grid_h_logger("... done.");
}