/**
 * @file legacyvtk.c
 * @author Pejvak Javaheri; pejvak.javaheri@mail.utoronto.ca
 * @brief Definitions for legacy VTK ouput
 */

#include "legacyvtk.h"

/**
 * (Internal)
 * Array mean calculator.
 */
double mean(double * a, int N) {
    double avg = 0.;
    for (int i = 0; i < N; i++) {
        avg += a[i];
    }
    return avg/N;
}

int make_rectangular_vtk_grid(
    char * adr,
    double * xs,
    double * ys,
    double * zs,
    int32_t i_max,
    int32_t j_max
) {
    FILE * fptr = fopen(adr, "w");

    if (fptr == NULL) {
        printf("Error opening file!\n");
        exit(1);
    }

    const int num_pts       = i_max*j_max;
    const int num_faces     = (i_max-1)*(j_max-1);
    const int num_floats    = num_faces*5;

    fprintf(
        fptr,
        "# vtk DataFile Version 3.0\nTwo triangles\nASCII\nDATASET POLYDATA\n\n"
    );
    fprintf(
        fptr, 
        "POINTS %d float\n",
        num_pts
    );

    for (int i = 0; i < i_max; i++) {
        for (int j = 0; j < j_max; j++) {
            fprintf(
                fptr,
                "%e %e %e\n",
                xs[i*j_max + j], ys[i*j_max + j], zs[i*j_max + j]
            );
        }
    }

    fprintf(
        fptr, 
        "\nPOLYGONS %d %d\n",
        num_faces, num_floats
    );

    for (int i = 0; i < i_max-1; i++) {
        for (int j = 0; j < j_max-1; j++) {
            fprintf(
                fptr,
                "4 %d %d %d %d \n",
                (i    )*j_max + (j    ), 
                (i + 1)*j_max + (j    ), 
                (i + 1)*j_max + (j + 1), 
                (i    )*j_max + (j + 1)
            );
        }
    }

    fprintf(
        fptr, 
        "\nPOINT_DATA %d\n",
        num_pts
    );

    fclose(fptr);
    return 0;
}

int add_rectangular_vtk_field(
    char * adr,
    char * field_name,
    double * field,
    int32_t i_max,
    int32_t j_max
) {
    FILE * fptr = fopen(adr, "a");

    if (fptr == NULL) {
        printf("Error opening file!\n");
        exit(1);
    }
    fprintf(
        fptr,
        "\nSCALARS %s float 1\nLOOKUP_TABLE default\n", 
        field_name
    );
    for (int i = 0; i < i_max; i++) {
        for (int j = 0; j < j_max; j++) {
            fprintf(
                fptr,
                "%e\n", 
                field[i*j_max + j]
            );
        }
    }

    fclose(fptr);

    return 0;
}


int make_spherical_vtk_grid(
    char * adr,
    double * xs,
    double * ys,
    double * zs,
    int32_t i_max,
    int32_t j_max
) {
    FILE * fptr = fopen(adr, "w");

    if (fptr == NULL) {
        printf("Error opening file!\n");
        exit(1);
    }

    const double npolex = mean(xs, i_max);
    const double npoley = mean(ys, i_max);
    const double npolez = mean(zs, i_max);
    const double spolex = mean(xs + (i_max-1)*j_max, i_max);
    const double spoley = mean(ys + (i_max-1)*j_max, i_max);
    const double spolez = mean(zs + (i_max-1)*j_max, i_max);

    const int N = i_max - 2;
    const int M = j_max - 1;

    const int num_pts       = N*M + 2;
    const int num_faces     = (N + 1)*M;
    const int num_floats    = (N-1)*(M)*5 + 2*M*4;

    const int npole_ID = N*M;
    const int spole_ID = npole_ID + 1;

    fprintf(
        fptr,
        "# vtk DataFile Version 3.0\nTwo triangles\nASCII\nDATASET POLYDATA\n\n"
    );
    fprintf(
        fptr, 
        "POINTS %d float\n",
        num_pts
    );

    for (int i = 1; i < i_max-1; i++) {
        for (int j = 0; j < j_max-1; j++) {
            fprintf(
                fptr,
                "%e %e %e\n",
                xs[i*j_max + j], ys[i*j_max + j], zs[i*j_max + j]
            );
        }
    }
    fprintf(fptr, "%e %e %e\n", npolex, npoley, npolez);
    fprintf(fptr, "%e %e %e\n", spolex, spoley, spolez);

    fprintf(
        fptr, 
        "\nPOLYGONS %d %d\n",
        num_faces, num_floats
    );

    for (int ii = 0; ii < N-1; ii++) {
        for (int jj = 0; jj < M-1; jj++) {
            fprintf(
                fptr,
                "4 %d %d %d %d \n",
                (ii    )*M + (jj    ), 
                (ii + 1)*M + (jj    ), 
                (ii + 1)*M + (jj + 1), 
                (ii    )*M + (jj + 1)
            );
        }
        const int jj = M - 1;
        fprintf(
            fptr,
            "4 %d %d %d %d \n",
            (ii    )*M + (jj    ), 
            (ii + 1)*M + (jj    ), 
            (ii + 1)*M + (0     ), 
            (ii    )*M + (0     )
        );
    }

    // north pole
    for (int jj = 0; jj < M-1; jj++) {
        const int ii = 0;
        fprintf(
            fptr,
            "3 %d %d %d \n", 
            npole_ID, 
            (ii    )*M + (jj    ), 
            (ii    )*M + (jj + 1)
        );
    }
    fprintf(
        fptr,
        "3 %d %d 0\n", 
        npole_ID, 
        M-1
    );

    // south pole
    for (int jj = 0; jj < M-1; jj++) {
        const int ii = N-1;
        fprintf(
            fptr,
            "3 %d %d %d \n", 
            spole_ID, 
            (ii    )*M + (jj + 1), 
            (ii    )*M + (jj    )
        );
    }
    fprintf(
        fptr,
        "3 %d %d %d\n", 
        spole_ID, 
        (N-1)*M,
        (N-1)*M + M-1
    );

    fprintf(
        fptr, 
        "\nPOINT_DATA %d\n",
        num_pts
    );

    fclose(fptr);
    return 0;
}

int add_spherical_vtk_field(
    char * adr,
    char * field_name,
    double * field,
    int32_t i_max,
    int32_t j_max
) {
    FILE * fptr = fopen(adr, "a");

    if (fptr == NULL) {
        printf("Error opening file!\n");
        exit(1);
    }
    fprintf(
        fptr,
        "\nSCALARS %s float 1\nLOOKUP_TABLE default\n", 
        field_name
    );
    for (int i = 1; i < i_max-1; i++) {
        for (int j = 0; j < j_max-1; j++) {
            fprintf(
                fptr,
                "%e\n", 
                field[i*j_max + j]
            );
        }
    }
    fprintf(fptr, "%e\n", mean(field, i_max));
    fprintf(fptr, "%e\n", mean(field + (i_max-1)*j_max, i_max));

    fclose(fptr);

    return 0;
}
