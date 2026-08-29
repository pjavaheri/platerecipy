/**
 * @file grid.h
 * @author Pejvak Javaheri; pejvak.javaheri@mail.utoronto.ca
 * @brief Header file for grid structures
 */

#ifndef CLIB_GRID_H
#define CLIB_GRID_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>

#ifdef _WIN32
    #define CLIB_EXPORT __declspec(dllexport)
#else
    #define CLIB_EXPORT
#endif


// ~~~~~ Callback mechanism for the Python logger ~~~~~

typedef void (*grid_h_log_func)(const char *);
static grid_h_log_func grid_h_logger = NULL;
CLIB_EXPORT void set_grid_h_logger(grid_h_log_func func);

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





// ~~~~~~~~~~~~~~~ Grid type identifiers ~~~~~~~~~~~~~~~

// Custom/unknown grid
const int32_t GTYPE_CUSTOM               = -1;

// Flat rectilinear grid
const int32_t GTYPE_PLANAR               = 0;

// Partial rectilinear spherical grid
const int32_t GTYPE_PARTIAL_SPHERICAL    = 1;

// Spherical grid (with wraparound and polar nodes)
const int32_t GTYPE_SPHERICAL            = 2;

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



// ~~~~~~~~~~~~~~~~~~~~ Exposed API ~~~~~~~~~~~~~~~~~~~~

// Node structure
typedef struct Node {
    int32_t     idx;            // main index along Map.nodes[idx]
    int32_t     npy_idx;        // the index on the numpy array
    int32_t     ord_idx;        // the index after RW reordering
    int32_t *   neighs;         // pointer to a sub-array of neighbouring indices
    double *    edge_lengths;   // pointer to a sub-array of connection lengths to neighbouring indices
    int32_t     num_neighs;     // number of neigbours
} Node;


// Map/Grid structure
typedef struct Map {
    Node *      nodes;          // list of connected nodes comprising the grid
    int32_t     num_nodes;      // number of nodes
    int32_t     num_edges;      // number of connections (i.e., gridlines)
    int32_t     grid_type;      // grid type (e.g., GTYPE_######)
    int32_t     i_max;          // first-index extent of the 2D layout, if NA, is -1
    int32_t     j_max;          // second-index extent of the 2D layout, if NA, is -1
} Map;


/**
 * Returns a node at index i in the map/grid.
 *
 * @note This is mostly intended to be used from Python API, as from C,
 * accessing the i-th node is trivial (i.e., map->node[i]).
 *
 * @param map a pointer to a map struct
 * @param i index of the node
 * 
 * @returns a pointer to the node struct or NULL if out-of-bound
 */
CLIB_EXPORT Node * get_node_at(
    Map *       map,    // a pointer to a map struct
    int32_t     i       // index of the node
);


/**
 * Returns a node at index i in the node list/array.
 *
 * @note This is mostly intended to be used from Python API, as from C,
 * accessing the i-th node is trivial (i.e., node[i]).
 *
 * @param nodes a pointer to a node array
 * @param i index of the node
 * 
 * @returns index of the neighbour or -1 if out-of-bound
 */
CLIB_EXPORT int32_t get_neigh_at(
    Node *      nodes,      // a pointer to a node array
    int32_t     i           // index of the node
);


/**
 * Sets numpy indices to the nodes that comprise the grid.
 *
 * @param map a pointer to map
 * @param npy_idxs a pointer to an array of indicies
 */
CLIB_EXPORT void set_npy_idx(
    Map *       map,
    int32_t *   npy_idxs
); 


/**
 * Allocate and initialize a map by reading all the cells 
 * within a mesh, by reading a cells array compressed as:
 * cells: 3, 0, 1, 2, 4, 12, 13, 14, 15, ...
 * where 3 shows the cell has three nodes (0, 1, and 2), 
 * and so on.
 *
 * @param num_nodes number of nodes in the mesh (also one greater than the largest node id)
 * @param cells a pointer to the cells array
 * @param num_cells lengths of the cells array (and not the number of cells)
 */
CLIB_EXPORT Map * alloc_map_from_cells(
    int32_t     num_nodes,
    int32_t *   cells,
    int32_t     num_cells
);


/**
    Basic Rectilinear Map

    connection pattern:
        
    O -- O -- O -- O -- O
    |    |    |    |    |
    O -- O -- O -- O -- O
    |    |    |    |    |
    O -- O -- O -- O -- O
    |    |    |    |    |
    O -- O -- O -- O -- O
 */
CLIB_EXPORT Map * alloc_rect_map(
    int32_t i_max,
    int32_t j_max,
    int32_t grid_type
);

/**
    Basic Spherical Map

    connection pattern:
        
          / -- / -- / -- / -- O
          |    |    |    |    |
    ... - O -- O -- O -- O -- O - ...
          |    |    |    |    |
    ... - O -- O -- O -- O -- O - ...
          |    |    |    |    |
          O -- / -- / -- / -- /
 */
CLIB_EXPORT Map * alloc_sph_map(
    int32_t i_max,
    int32_t j_max
);

/**
 * Free all heap-allocated arrays embedded directly or indirectly 
 * within the map struct.
 * 
 * @param map
 */
CLIB_EXPORT void free_map(
    Map * map
);


/**
 * Populates cells array with cells that comprise 
 * either known (e.g., rect or sph) maps. Otherwise,
 * the output is filled with -1.
 *  
 * 
 * @param map a pointer to a Map struct 
 * @param cells a pointer to the output array
 * 
 * @warning cells must be alreay allocated and will be overwritten
 */
CLIB_EXPORT void gen_cells_from_map(
    Map *       map,
    int32_t *   cells
);


#endif