#include "engine.h"

int main(int argc, char* argv[])
{
    if (SDL_Init(SDL_INIT_VIDEO) != 0)
    {
        printf("Error initializing SDL: %s\n", SDL_GetError());
        return 0;
    }
    printf("SDL successfully initialized!\n");

    SDL_Window* win = SDL_CreateWindow("Main Window",
                                       SDL_WINDOWPOS_CENTERED,
                                       SDL_WINDOWPOS_CENTERED,
                                       480, 360, 0);
    Uint32 render_flags = SDL_RENDERER_ACCELERATED;
    SDL_Renderer* rend = SDL_CreateRenderer(win, -1, render_flags);
    SDL_SetRelativeMouseMode(SDL_TRUE);
    SDL_bool done = SDL_FALSE;

    point_t point1T1 = {
        .x = 1, 
        .y = 1, 
        .z = 1
    };
    point_t point2T1 = {
        .x = 1, 
        .y = -1, 
        .z = 1
    };
    point_t point3T1 = {
        .x = 0, 
        .y = 1, 
        .z = 1
    };
    point_t point1T2 = {
        .x = 1, 
        .y = 1, 
        .z = 0
    };
    point_t point2T2 = {
        .x = 1, 
        .y = -1, 
        .z = 1.5
    };
    point_t point3T2 = {
        .x = 0, 
        .y = 1, 
        .z = 1.5
    };
    triangle_t triangle1 = {
        .p1 = point1T1,
        .p2 = point2T1,
        .p3 = point3T1,
        .color.r = 0,
        .color.g = 255,
        .color.b = 0,
        .color.a = SDL_ALPHA_OPAQUE,
    };
    triangle_t triangle2 = {
        .p1 = point1T2,
        .p2 = point2T2,
        .p3 = point3T2,
        .color.r = 0,
        .color.g = 0,
        .color.b = 255,
        .color.a = SDL_ALPHA_OPAQUE,
    };
    triangle_t triangle3 = {
        .p1 = (point_t){-1, -1, -1},
        .p2 = (point_t){-1, -1, 0.5},
        .p3 = (point_t){0.5, -1, -1},
        .color.r = 255,
        .color.g = 0,
        .color.b = 0,
        .color.a = SDL_ALPHA_OPAQUE,
    };
    light_t light = {
        .position.x = 0, 
        .position.y = 0, 
        .position.z = 1, 
        .intensity = 100
    };
    camera_t camera = {
        .position.x = 0, 
        .position.y = 0, 
        .position.z = 0, 
        .rotation.x = 90, 
        .rotation.y = 0, 
        .rotation.z = 0,
        .viewHeight = 360,
        .viewWidth = 480,
        .focalLength = 1
    };
    MESH_INIT(mesh);
    void *pTriangle1 = &triangle1;
    void *pTriangle2 = &triangle2;
    meshAddTris(&mesh, pTriangle1);
    meshAddTris(&mesh, pTriangle2);
    meshAddTris(&mesh, (void*)&triangle3);

    while (!done) 
    {
        SDL_Event event;

        SDL_SetRenderDrawColor(rend, 0, 0, 0, SDL_ALPHA_OPAQUE);
        SDL_RenderClear(rend);

        //drawLine(rend, 0, 0, 480, 360, 0, 122);
        //renderScreen(rend, camera, triangle, light);
        startTimer();
        renderMesh(rend, camera, mesh, light);
        printf("rendered mesh in %ds\n", stopTimer());

        /*SDL_SetRenderDrawColor(rend, 255, 255, 255, SDL_ALPHA_OPAQUE);
        SDL_RenderDrawLine(rend, 320, 200, 300, 240);
        SDL_RenderDrawLine(rend, 300, 240, 340, 240);
        SDL_RenderDrawLine(rend, 340, 240, 320, 200);*/
        SDL_RenderPresent(rend);

        while (SDL_PollEvent(&event)) 
        {
            switch (event.type)
            {
            case SDL_QUIT:
                done = SDL_TRUE;
                break;
            
            case SDL_KEYDOWN:
                switch (event.key.keysym.sym)
                {
                case SDLK_LEFT:
                    camera.rotation.y -= 0.1;
                    break;
                
                case SDLK_RIGHT:
                    camera.rotation.y += 0.1;
                    break;
                
                case SDLK_UP:
                    camera.rotation.x += 0.1;
                    break;
                
                case SDLK_DOWN:
                    camera.rotation.x -= 0.1;
                    break;

                case SDLK_w:
                    triangle1.p1.x -= .5;
                    break;

                case SDLK_s:
                    triangle1.p1.x += .5;
                    break;
                
                default:
                    break;
                }
                break;

            case SDL_MOUSEMOTION:
                camera.rotation.x -= event.motion.yrel * 0.001;
                camera.rotation.y += event.motion.xrel * 0.001;
                break;

            case SDL_MOUSEWHEEL:
                camera.focalLength += event.wheel.y;
                break;
            
            default:
                break;
            }            
        }
    }

    SDL_Quit();
    return 0;
}
