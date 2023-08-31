#ifndef ENGINE_H
#define ENGINE_H

#include <string.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include "SDL2/include/SDL2/SDL.h"
#include "debug.h"

#define debug false
#define warnings false
#define DEFAULT_COLOUR ((color_t){.r=0, .g=255, .b=255, .a=SDL_ALPHA_OPAQUE})
#define MESH_INIT(mesh) mesh_t mesh; meshInitTris(&mesh)
#define MESH_INIT_CAPACITY 1
#define UNDEFINED -1
#define SUCCESS 0


/// @brief The equation for a plane (𝑎𝑥+𝑏𝑦+𝑐𝑧+𝑑=0)
typedef struct planeEquation_t
{
    double ax;
    double by;
    double cz;
    double d;
} planeEquation_t;

/// @brief A general purpose equation with 3 variables (𝑎+𝑏+𝑐=𝑑)
typedef struct equation3_t
{
    double a;
    double b;
    double c;
    double d;
} equation3_t;

/// @brief 3d-space coordinates (aliases: `point_t`, `rotation_t`)
struct xyz_t
{
    double x;
    double y;
    double z;
};

/// @brief A point in 3d space.
typedef struct xyz_t point_t;

/// @brief Rotation in 3d space (degrees).
typedef struct xyz_t rotation_t;

/// @brief A RGBA color
typedef struct color_t
{
    Uint8 r; //red
    Uint8 g; //green
    Uint8 b; //blue
    Uint8 a; //alpha
} color_t;

/// @brief A triangle in 3d-space
typedef struct triangle_t
{
    point_t p1;
    point_t p2;
    point_t p3;
    color_t color;
} triangle_t;

/// @brief A light in 3d-space
typedef struct light_t
{
    point_t position;
    double intensity;
} light_t;

/// @brief A camera in 3d-space
typedef struct camera_t
{
    point_t position;
    rotation_t rotation;
    int viewWidth;
    int viewHeight;
    double focalLength;
} camera_t;

/// @brief A collection of 3 variables
typedef struct abc_t
{
    double a;
    double b;
    double c;
} abc_t;

/// @brief A collection of triangles
typedef struct triangleCollection_t
{
    void **items;
    int capacity;
    int length;
} triangleCollection_t;

/// @brief An object
typedef struct mesh_t
{
    triangleCollection_t tris;
    point_t position;
    rotation_t rotation;
} mesh_t;


/// @brief Initialize the mesh by allocating all the momory
/// @param mesh the mesh to initialize
void meshInitTris(mesh_t *mesh)
{
    //initialize the capacity and allocate the memorys
    mesh->tris.capacity = MESH_INIT_CAPACITY;
    mesh->tris.length = 0;
    mesh->tris.items = malloc(sizeof(void *) * mesh->tris.capacity);
}

/// @brief resize the amount of triangles allocated for the mesh
/// @param mesh the mesh with the triangles
/// @param capacity the new amount of triangles the mesh should have space for
/// @return status: -1 = undefined | 0 = success
int meshResizeTris(mesh_t *mesh, int capacity)
{
    int status = UNDEFINED;
    if (mesh)
    {
        void **items = realloc(mesh->tris.items, sizeof(void *) * capacity);
        if (items)
        {
            mesh->tris.items = items;
            mesh->tris.capacity = capacity;
            status = SUCCESS;
        }
    }
    return status;
}

/// @brief release the allocated memory
/// @param mesh the mesh with the triangles
/// @return status: -1 = undefined | 0 = success
int meshFreeTris(mesh_t *mesh)
{
    int status = UNDEFINED;
    if (mesh)
    {
        free(mesh->tris.items);
        mesh->tris.items = NULL;
        status = SUCCESS;
    }
    return status;
}

/// @brief delete an item at an index and move all the other items to fill the space
/// @param mesh the mesh with the triangles
/// @param index the index of the triangle to delete
/// @return status: -1 = undefined | 0 = success
int meshDeleteTris(mesh_t *mesh, int index)
{
    int status = UNDEFINED;
    if (mesh)
    {
        if ((index < 0) || (index >= mesh->tris.length))
            return status;
        mesh->tris.items[index] = NULL;
        for (int i = index; (i < mesh->tris.length - 1); ++i)
        {
            mesh->tris.items[i] = mesh->tris.items[i + 1];
            mesh->tris.items[i + 1] = NULL;
        }
        mesh->tris.length--;
        /*if ((mesh->tris.length > 0) && ((mesh->tris.length) == (mesh->tris.capacity / 4)))
        {
            meshResizeTris(mesh, mesh->tris.capacity / 2);
        }*/
        meshResizeTris(mesh, mesh->tris.length);
        status = SUCCESS;
    }
    return status;
}

/// @brief gets the triangle of a mesh at a specific index
/// @param mesh the mesh with the triangles
/// @param index the index of the triangle to get
/// @return the item at the index
void *meshGetTris(mesh_t *mesh, int index)
{
    void *readData = NULL;
    if (mesh)
    {
        if ((index >= 0) && (index < mesh->tris.length))
        {
            readData = mesh->tris.items[index];
        }
    }
    return readData;
}

/// @brief set the triangle of a mesh at a specific index
/// @param mesh the mesh with the triangles
/// @param index the index of the triangle to set
/// @param item the triangle to set to
/// @return status: -1 = undefined | 0 = success
int meshSetTris(mesh_t *mesh, int index, void *item)
{
    int status = UNDEFINED;
    if (mesh)
    {
        if ((index >= 0) && (index < mesh->tris.length))
        {
            mesh->tris.items[index] = item;
            status = SUCCESS;
        }
    }
    return status;
}

/// @brief adds one triangle to the mesh at the back and allocates more space to the mesh if needed
/// @param mesh the mesh with the triangles
/// @param item the triangle to add
/// @return status: -1 = undefined | 0 = success
int meshAddTris(mesh_t *mesh, void *item)
{
    int status = UNDEFINED;
    if (mesh)
    {
        if (mesh->tris.capacity == mesh->tris.length)
        {
            status = meshResizeTris(mesh, mesh->tris.capacity + 1);
            if (status != UNDEFINED)
            {
                mesh->tris.items[mesh->tris.length++] = item;
            }
        }
        else
        {
            mesh->tris.items[mesh->tris.length++] = item;
            status = SUCCESS;
        }
    }
    return status;
}

/// @brief gets the amount of triangles in the mesh
/// @param mesh the mesh with the triangles
/// @return amount of triangles
int meshTotalTris(mesh_t *mesh)
{
    int total = UNDEFINED;
    if (mesh)
    {
        total = mesh->tris.length;
    }
    return total;
}

// timer (global variable)
time_t timer;

/// @brief starts the timer (stop the timer with "stopTimer()")
void startTimer()
{
    timer = time(NULL);
}

/// @brief stop the timer (start the timer with "startTimer()")
/// @return the seconds that passed since startTimer()
int stopTimer()
{
    return difftime(time(NULL), timer);
}

/// @brief Solves a system of three equations for 𝑎, 𝑏 and 𝑐. Time complexity: O(n³)
/// @param e1 equation one
/// @param e2 equation two
/// @param e3 equation three
/// @return 𝑎, 𝑏 and 𝑐
abc_t gaussianEliminationForThreeUnknowns(equation3_t e1, equation3_t e2, equation3_t e3)
{
    // coefficient matrix
    double A[3][3] = {
        { e1.a, e1.b, e1.c },
        { e2.a, e2.b, e2.c },
        { e3.a, e3.b, e3.c }
    };
    // right side matrix
    double B[3] = { e1.d, e2.d, e3.d };
    double tmpRow[3];
    double tmpVar;
    double factor;
    double sumAX;
    
    // Gaußsches Eliminationsverfahren:
    int n = 3;
    double X[3] = { 0, 0, 0 };

    //first loop specifys the fixed row
    for (int k = 0; k < n-1; k++)
    {
        if (abs(A[k][k]) < 1.0e-12)
        {
            for (int i = k + 1; i < n; i++)
            {
                if (abs(A[i][k]) > abs(A[k][k]))
                {
                    tmpRow[0] = A[k][0];
                    tmpRow[1] = A[k][1];
                    tmpRow[2] = A[k][2];
                    A[k][0] = A[i][0];
                    A[k][1] = A[i][1];
                    A[k][2] = A[i][2];
                    A[i][0] = tmpRow[0];
                    A[i][1] = tmpRow[1];
                    A[i][2] = tmpRow[2];
                    
                    tmpVar = B[k];
                    B[k] = B[i];
                    B[i] = tmpVar;

                    break;
                }
                
            }
            
        }

        //applies the elimination below the fixed row
        for (int i = k + 1; i < n; i++)
        {
            if (A[i][k] == 0)
            {
                continue;
            }

            factor = A[k][k] / A[i][k];
            for (int j = k; j < n; j++)
            {
                A[i][j] = A[k][j] - A[i][j] * factor;
            }

            //we also calculate the b vector of each row
            B[i] = B[k] - B[i] * factor;
        }
        
        
    }

    X[n - 1] = B[n - 1] / A[n - 1][n - 1];
    for (int i = n - 2; i > -1; i--)
    {
        sumAX = 0;

        for (int j = i + 1; j < n; j++)
        {
            sumAX += A[i][j] * X[j];
        }
        
        X[i] = (B[i] - sumAX) / A[i][i];
    }

    return (abc_t) {
        .a = X[0],
        .b = X[1],
        .c = X[2]
    };
}

/// @brief rotate point around the x axis using an origin
/// @param point the point to rotate
/// @param angle the angle of rotation
/// @param origin the point to rotate around
/// @return rotated point
point_t rotatePointX(point_t point, double angle, point_t origin)
{
    /* formulas: 
    x′ = x × cos(β) − y × sin(β)
    y′ = y × cos(β) + x × sin(β)        */
    double sinAngle = sin(angle);
    double cosAngle = cos(angle);

    // subract the origins transformation
    point.x -= origin.x;
    point.y -= origin.y;
    point.z -= origin.z;

    point_t rotatedPoint = {
        .x = point.x,
        .y = point.y * cosAngle - point.z * sinAngle,
        .z = point.z * cosAngle + point.y * sinAngle
    };

    // add the origins transformation
    rotatedPoint.x += origin.x;
    rotatedPoint.y += origin.y;
    rotatedPoint.z += origin.z;

    return rotatedPoint;
}

/// @brief rotate point around the y axis using an origin
/// @param point the point to rotate
/// @param angle the angle of rotation
/// @param origin the point to rotate around
/// @return rotated point
point_t rotatePointY(point_t point, double angle, point_t origin)
{
    /* formulas: 
    x′ = x × cos(β) − y × sin(β)
    y′ = y × cos(β) + x × sin(β)        */
    double sinAngle = sin(angle);
    double cosAngle = cos(angle);

    // subract the origins transformation
    point.x -= origin.x;
    point.y -= origin.y;
    point.z -= origin.z;

    point_t rotatedPoint = {
        .x = point.x * cosAngle - point.z * sinAngle,
        .y = point.y,
        .z = point.z * cosAngle + point.x * sinAngle
    };

    // add the origins transformation
    rotatedPoint.x += origin.x;
    rotatedPoint.y += origin.y;
    rotatedPoint.z += origin.z;

    return rotatedPoint;
}

/// @brief rotate point around the z axis using an origin
/// @param point the point to rotate
/// @param angle the angle of rotation
/// @param origin the point to rotate around
/// @return rotated point
point_t rotatePointZ(point_t point, double angle, point_t origin)
{
    /* formulas: 
    x′ = x × cos(β) − y × sin(β)
    y′ = y × cos(β) + x × sin(β)        */
    double sinAngle = sin(angle);
    double cosAngle = cos(angle);

    // subract the origins transformation
    point.x -= origin.x;
    point.y -= origin.y;
    point.z -= origin.z;

    point_t rotatedPoint = {
        .x = point.x * cosAngle - point.y * sinAngle,
        .y = point.y * cosAngle + point.x * sinAngle,
        .z = point.z
    };

    // add the origins transformation
    rotatedPoint.x += origin.x;
    rotatedPoint.y += origin.y;
    rotatedPoint.z += origin.z;

    return rotatedPoint;
}

/// @brief rotate a point around an origin point
/// @param point 
/// @param rotation 
/// @param origin 
/// @return rotated point
point_t rotatePoint(point_t point, rotation_t rotation, point_t origin)
{
    point_t rotatedPoint = rotatePointX(point, rotation.x, origin);
    rotatedPoint = rotatePointY(rotatedPoint, rotation.y, origin);
    rotatedPoint = rotatePointZ(rotatedPoint, rotation.z, origin);

    return rotatedPoint;
}

/// @brief Function for calculating the second binomial formula
/// @param a
/// @param b
/// @return result
static double binF(double a, double b)
{
    return pow(a, 2) - 2 * a * b + pow(b, 2);
}

/// @brief Function for calculating the distance between two points in 3d-space
/// @param p0 the first point
/// @param p1 the second point
/// @return distance
double calculateDistance(point_t p0, point_t p1)
{
    return sqrt(binF(p1.x, p0.x) + binF(p1.y, p0.y) + binF(p1.z, p0.z));
}

/// @brief Function for calculating the level of luminance based on the distance to the light source
/// @param p the point
/// @param l the light
/// @return luminance
Uint8 calculateLuminance(point_t p, light_t l)
{
    double distance = calculateDistance(p, l.position);
    if (distance == 0) { return 0; } // prevent div by 0

    double luminance = l.intensity / pow(distance, 2);
    
    if (debug)
    {
        printf("luminance: %f\n", luminance);
    }
    
    if (luminance > 255) { luminance = 255; }
    return round(luminance);
}

/// @brief Function for drawing a line between two points
/// @param rend sdl renderer
/// @param x0 x-coordinate of the first point
/// @param y0 y-coordinate of the first point
/// @param x1 x-coordinate of the second point
/// @param y1 y-coordinate of the second point
/// @param l0 the luminance level of the first point
/// @param l1 the luminance level of the second point
void drawLine(SDL_Renderer *rend, int x0, int y0, int x1, int y1, Uint8 l0, Uint8 l1)
{
    /* source: https://de.wikipedia.org/wiki/Bresenham-Algorithmus */

    int dx =  abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
    int dy = -abs(y1 - y0), sy = y0 < y1 ? 1 : -1;
    int err = dx + dy, e2; /* error value e_xy */
    float length = abs(x1 - x0) > abs(y1 - y0) ? abs(x1 - x0) : abs(y1 - y0);
    float increment = l0 < l1 ? -(l1 / length) : -(l0 / length);
    float l = l0 > l1 ? l0 : l1;

    while (1) {
        if (x0 == x1 && y0 == y1) break;
        e2 = 2 * err;
        if (e2 > dy) { err += dy; x0 += sx; } /* e_xy+e_x > 0 */
        if (e2 < dx) { err += dx; y0 += sy; } /* e_xy+e_y < 0 */
        l += increment;
        SDL_SetRenderDrawColor(rend, round(l), round(l), round(l), SDL_ALPHA_OPAQUE);
        SDL_RenderDrawPoint(rend, x0, y0);
    }
}

/// @brief calculates the cross product of two 3d vectors
/// @param a 3d vector
/// @param b 3d vector
/// @return 3d vector
point_t crossProduct(point_t a, point_t b)
{
    point_t s = {
        .x = a.y * b.z - a.z * b.y,
        .y = a.z * b.x - a.x * b.z,
        .z = a.x * b.y - a.y * b.x
    };
    
    return s;
}

/// @brief creates the equation for a plane from 3 non-linear points (𝑎𝑥+𝑏𝑦+𝑐𝑧+𝑑=0)
/// @param i the first point
/// @param j the second point
/// @param k the third point
/// @return plane equation
planeEquation_t createPlaneEquation(point_t i, point_t j, point_t k)
{
    // very helpful: https://mathinsight.org/forming_plane_examples

    point_t b = {
        .x = i.x - j.x,
        .y = i.y - j.y,
        .z = i.z - j.z
    }, c = {
        .x = i.x - k.x,
        .y = i.y - k.y,
        .z = i.z - k.z
    };

    // normal vector: n = b × c (cross product)
    point_t n = crossProduct(b, c);    

    // solve the equation: n · (𝑥 - ix, 𝑦 - iy, 𝑧 - iz) = 0
    planeEquation_t result = {
        .ax = n.x,
        .by = n.y,
        .cz = n.z,
        .d  = n.x * (-i.x) + n.y * (-i.y) + n.z * (-i.z)
    };

    return result;
}

/// @brief Calculates the intersection point I for a ray from C to P with a plane.
/// @param C the Camera as a point
/// @param P the Pixel as a point
/// @param E the Equation for the plane
/// @return the intersection I as a point of the plane and the ray.
point_t calculateIntersection(point_t C, point_t P, planeEquation_t E)
{
    double aPart = E.ax * (C.x + P.x);
    double bPart = E.by * (C.y + P.y);
    double cPart = E.cz * (C.z + P.z);
    double tCoefficient = aPart + bPart + cPart;
    E.d += C.x + C.y + C.z;
    if ( warnings && tCoefficient == 0 ) 
    { 
        printf("Division by zero in function 'calculateIntersection'! (engine.h)\n");
    }

    double t = (-E.d) / tCoefficient;
    point_t intersection = {
        .x = (1 - t) * C.x + t * P.x,
        .y = (1 - t) * C.y + t * P.y,
        .z = (1 - t) * C.z + t * P.z
    };
    return intersection;
}

/// @brief calculates if I ∈ △ABC (△ABC = T)
/// @details First the function computes the weights for the point I in the triangle to the 3 points of the triangle A, B and C. If one of the weights is negative then the point is not intersecting the triangle.
/// @param I point to calculate if inside the triangle
/// @param T the triangle
/// @return boolean value (true or false)
bool isPointInTriangle(point_t I, triangle_t T)
{
    point_t A = T.p1;
    point_t B = T.p2;
    point_t C = T.p3;

    // Ix = aAx + bBx + cCx
    if (debug) {
        printf("%fa + %fb + %fc = %f\n", A.x, B.x, C.x, I.x);
    }
    equation3_t equation1 = { .a = A.x, .b = B.x, .c = C.x, .d = I.x };

    // Iy = aAy + bBy + cCy
    if (debug) {
        printf("%fa + %fb + %fc = %f\n", A.y, B.y, C.y, I.y);
    }
    equation3_t equation2 = { .a = A.y, .b = B.y, .c = C.y, .d = I.y };

    // Iz = aAz + bBz + cCz
    if (debug) {
        printf("%fa + %fb + %fc = %f\n", A.z, B.z, C.z, I.z);
    }
    equation3_t equation3 = { .a = A.z, .b = B.z, .c = C.z, .d = I.z };

    // calculate the weights
    abc_t variables = gaussianEliminationForThreeUnknowns(equation1, equation2, equation3);
    double a = variables.a;
    double b = variables.b;
    double c = variables.c;

    if (debug)
    {
        printf("weights:\na: %f\nb: %f\nc: %f\n", a, b, c);
    }

    if (a >= 0 && b >= 0 && c >= 0) { 

        return true; 

    } else {
        
        return false;

    }

}

mesh_t parseObj()
{
    FILE *fp;
    char buff[255];
    char *portion;
    MESH_INIT(mesh);
    point_t point;
    int index = 0;
    int triangleLength = 0;
    int triangleCapacity = 0;
    long pointListLength = 0;
    long pointListCapacity = 0;
    /*triangle_t triangle;
    triangle_t *pTriangle = &triangle;*/

    if ((fp = fopen("F:\\Coding\\C\\ScreenProject\\GUI\\assets\\models\\cube.obj", "r")) == NULL)
    {
        printf("Error! opening file");
        exit(EXIT_FAILURE);
    }

    while (!feof(fp))
    {
        fgets(buff, 255, fp);
        switch (buff[0])
        {
        case 'v':
            pointListCapacity++;
            break;
        
        case 'f':
            triangleCapacity++;
            break;
        
        default:
            break;
        }
    }
    fseek(fp, 0, SEEK_SET);

    meshResizeTris(&mesh, triangleCapacity);
    point_t *pointList = (point_t *) calloc(pointListCapacity, sizeof(point_t));
    if (pointList == NULL)
    {
        printf("Unable to allocate memory\n");
        exit(EXIT_FAILURE);
    }

    triangle_t *triangleArray = (triangle_t *) calloc(triangleCapacity, sizeof(triangle_t));
    if (triangleArray == NULL)
    {
        printf("Unable to allocate memory\n");
        exit(EXIT_FAILURE);
    }

    while (!feof(fp))
    {
        fgets(buff, 255, fp);

        switch (buff[0])
        {
        case 'o':
            break;

        case 'v':
            portion = strtok(buff, " ");     // the operator part (o/v/s/f/...)
            portion = strtok(NULL, " ");     // the first parameter (x Coordinate)
            point.x = strtod(portion, NULL); // convert string to double
            portion = strtok(NULL, " ");     // the second parameter (y Coordinate)
            point.y = strtod(portion, NULL); // convert string to double
            portion = strtok(NULL, " ");     // the third parameter (z Coordinate)
            point.z = strtod(portion, NULL); // convert string to double
            portion = strtok(NULL, " ");     // set the pointer to NULL
            
            //printf("point %ld: (%lf, %lf, %lf)\n", pointListLength, point.x, point.y, point.z);
            pointList[pointListLength] = point;
            pointListLength++;

            break;

        case 's':
            break;

        case 'f':
            //printf("buff: %s\n", buff);
            //printf("triangle pointer: %p\n", pTriangle);
            portion = strtok(buff, " ");                        // the operator part (o/v/s/f/...)
            //printf("%s\n", portion);
            portion = strtok(NULL, " ");                        // the first parameter (point 1)
            if (portion == NULL)
            {
                break;
            }
            
            //printf("%s\n", portion);
            //printf("p1: %ld\n", strtol(portion, NULL, 10));
            index = strtol(portion, NULL, 10) - 1;
            //printf("p1 index: %d\n", index);
            //printf("x: %lf\n", pointList[index].x);
            //printf("y: %lf\n", pointList[index].y);
            //printf("z: %lf\n", pointList[index].z);
            triangleArray[triangleLength].p1 = pointList[index]; // put the point at index of the parameter into point 1 of the triangle
            //printf("P1\nx: %lf\ny: %lf\nz: %lf\n\n", triangle.p1.x, triangle.p1.y, triangle.p1.z);
            portion = strtok(NULL, " ");                        // the second parameter (point 2)
            //printf("%s\n", portion);
            //printf("p2: %ld\n", strtol(portion, NULL, 10));
            index = strtol(portion, NULL, 10) - 1;
            triangleArray[triangleLength].p2 = pointList[index]; // put the point at index of the parameter into point 2 of the triangle
            //printf("P2\nx: %lf\ny: %lf\nz: %lf\n\n", triangle.p2.x, triangle.p2.y, triangle.p2.z);
            portion = strtok(NULL, " ");                        // the third parameter (point 3)
            //printf("%s\n", portion);
            //printf("p3: %ld\n", strtol(portion, NULL, 10));
            index = strtol(portion, NULL, 10) - 1;
            triangleArray[triangleLength].p3 = pointList[index]; // put the point at index of the parameter into point 2 of the triangle
            //printf("P3\nx: %lf\ny: %lf\nz: %lf\n\n", triangle.p3.x, triangle.p2.y, triangle.p3.z);
            portion = strtok(NULL, " ");                        // set the pointer to NULL
            triangleArray[triangleLength].color = DEFAULT_COLOUR;

            triangle_t *ptrTris = &triangleArray[triangleLength];
            meshAddTris(&mesh, ptrTris);
            triangleLength++;

            break;
        
        default:
            break;
        }

    }

    //printf("total: %d\n", meshTotalTris(&mesh));
    //printf("x: %lf\n", pointList[47 - 1].x);
    //printf("y: %lf\n", pointList[47 - 1].y);
    //printf("z: %lf\n", pointList[47 - 1].z);
    //printf("x: %lf\n", ((triangle_t *)meshGetTris(&mesh, 0))->p1.x);
    //printf("y: %lf\n", ((triangle_t *)meshGetTris(&mesh, 0))->p1.y);
    //printf("z: %lf\n", ((triangle_t *)meshGetTris(&mesh, 0))->p1.z);
    /*for (int i = 0; i < meshTotalTris(&mesh); i++)
    {
        triangle_t *t = meshGetTris(&mesh, i);
        printf("triangle %d:\n", i);
        printf("(%lf, %lf, %lf)\n", t->p1.x, t->p1.y, t->p1.z);
        printf("(%lf, %lf, %lf)\n", t->p2.x, t->p2.y, t->p2.z);
        printf("(%lf, %lf, %lf)\n\n", t->p3.x, t->p3.y, t->p3.z);
    }*/
    
    free(pointList);
    fclose(fp);

    return mesh;
}

void renderScreen(SDL_Renderer *rend, camera_t camera, triangle_t triangle, light_t light)
{
    for (int y = 0; y < camera.viewHeight; y++)
    {
        for (int x = 0; x < camera.viewWidth; x++)
        {
            // set the pixel
            point_t pixel = {
                .x = camera.position.x + 0.01 * (x - 0.5 * camera.viewWidth),
                .y = camera.position.y + 0.01 * (y - 0.5 * camera.viewHeight),
                .z = camera.position.z + camera.focalLength
            };

            // rotate the pixel around the camera
            pixel = rotatePoint(pixel, camera.rotation, camera.position);

            // create the equation of the plane from the triangle
            planeEquation_t planeEquation = createPlaneEquation(triangle.p1, triangle.p2, triangle.p3);

            // calculate the intersection with the plane using a ray from the camera to the pixel
            point_t intersection = calculateIntersection(camera.position, pixel, planeEquation);
            
            if (debug)
            {
                printf("Px%iy%i = (%f, %f, %f)\n", x, y, pixel.x, pixel.y, pixel.z);
                printf("Plane: %fx + %fy + %fz + %f = 0\n", planeEquation.ax, planeEquation.by, planeEquation.cz, planeEquation.d);
                printf("Ix%iy%i = (%f, %f, %f)\n", x, y, pixel.x, pixel.y, pixel.z);
            }

            // check whether the intersection is an element of the triangle
            if (isPointInTriangle(intersection, triangle))
            {
                // calculate the luminance of the point based on the distance to the light source
                Uint8 luminance = calculateLuminance(intersection, light);
                color_t color = {
                    .r = triangle.color.r / 255.0 * luminance,
                    .g = triangle.color.g / 255.0 * luminance,
                    .b = triangle.color.b / 255.0 * luminance,
                    .a = triangle.color.a
                };
                
                // draw the point
                SDL_SetRenderDrawColor(rend, color.r, color.g, color.b, color.a);
                SDL_RenderDrawPoint(rend, x, y);

                if (debug)
                {
                    printf("drew point at %i %i (color: %i_%i_%i)\n", x, y, color.r, color.g, color.b);
                }

            }

        }
        
    }
    
}

/// @brief render a mesh onto the screen using a SDL-renderer. (CPU computed raytracing)
/// @param rend 
/// @param camera 
/// @param mesh 
/// @param light 
void renderMesh(SDL_Renderer *rend, camera_t camera, mesh_t mesh, light_t light)
{
    // prevent weirdness by stopping the function if there are no triangles in the mesh
    if (mesh.tris.length == 0)
    {
        return;
    }
    
    // calculate everything that stays constant
    int totalTris = meshTotalTris(&mesh);

    // allocate memory for a dynamic array of all the plane equations
    planeEquation_t *planeEquationList = (planeEquation_t *) calloc(totalTris, sizeof(planeEquation_t));
    if (planeEquationList == NULL)
    {
        if (warnings)
        {
            printf("Unable to allocate memory\n");
        }
        return;
    }

    // precalculate all the plane equations to not waste resources in the loop
    for (int i = 0; i < totalTris; i++)
    {
        triangle_t *pTriangle = (triangle_t *)meshGetTris(&mesh, i);
        planeEquation_t planeEquation = createPlaneEquation(pTriangle->p1, pTriangle->p2, pTriangle->p3);
        planeEquationList[i] = planeEquation;
    }

    for (int y = 0; y < camera.viewHeight; y++)
    {
        for (int x = 0; x < camera.viewWidth; x++)
        {
            // set the pixel
            point_t pixel = {
                .x = camera.position.x + 0.01 * (x - 0.5 * camera.viewWidth),
                .y = camera.position.y + 0.01 * (y - 0.5 * camera.viewHeight),
                .z = camera.position.z + camera.focalLength
            };

            // rotate the pixel around the camera
            pixel = rotatePoint(pixel, camera.rotation, camera.position);

            // FIXME: this for loop takes way too long
            int index;
            double distance = DBL_MAX;
            double temporaryDistance;
            point_t intersection, temporaryIntersection;
            for (int i = 0; i < totalTris; i++)
            {
                // calculate the intersection with the plane using a ray from the camera to the pixel
                temporaryIntersection = calculateIntersection(camera.position, pixel, planeEquationList[i]);
                temporaryDistance = calculateDistance(camera.position, temporaryIntersection);

                if (temporaryDistance < distance)
                {
                    // calculate the intersection with the plane using a ray from the camera to the pixel
                    intersection = temporaryIntersection;
                    distance = temporaryDistance;
                    index = i;
                }
                
            }

            // check whether the intersection is an element of the triangle
            triangle_t *triangle = (triangle_t *)meshGetTris(&mesh, index);
            if (isPointInTriangle(intersection, *triangle))
            {
                // calculate the luminance of the point based on the distance to the light source
                Uint8 luminance = calculateLuminance(intersection, light);
                color_t color = {
                    .r = triangle->color.r / 255.0 * luminance,
                    .g = triangle->color.g / 255.0 * luminance,
                    .b = triangle->color.b / 255.0 * luminance,
                    .a = triangle->color.a
                };
                
                // draw the point
                SDL_SetRenderDrawColor(rend, color.r, color.g, color.b, color.a);
                SDL_RenderDrawPoint(rend, x, y);

                if (debug)
                {
                    printf("drew point at %i %i (color: %i_%i_%i)\n", x, y, color.r, color.g, color.b);
                }

            }

        }
        
    }
    
    free(planeEquationList);
    return;
}

#endif