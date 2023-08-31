/// @brief Debug-Function for printing a 3x3 2d-array to standard output.
/// @param matrix 3x3 2d-array
void print3x3Matrix(double matrix[3][3])
{
    for (int y = 0; y < 3; y++)
    {
        printf("[ ");

        for (int x = 0; x < 3; x++)
        {
            printf("%f", matrix[y][x]);
            if (x != 3)
            {
                printf(",");
            }
            printf(" ");
        }

        printf("]");
        if (y != 2)
        {
            printf(",");
        }
        printf("\n");
    }
}

/*
/// @brief Debug-Function for printing all the triangles of a mesh.
/// @param mesh mesh
void printMesh(mesh_t mesh)
{
    for (int i = 0; i < meshTotalTris(&mesh); i++)
    {
        triangle_t *t = meshGetTris(&mesh, i);
        printf("triangle %d:\n", i);
        printf("(%lf, %lf, %lf)\n", t->p1.x, t->p1.y, t->p1.z);
        printf("(%lf, %lf, %lf)\n", t->p2.x, t->p2.y, t->p2.z);
        printf("(%lf, %lf, %lf)\n\n", t->p3.x, t->p3.y, t->p3.z);
    }
}
*/