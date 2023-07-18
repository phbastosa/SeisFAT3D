# include "ray_tracing.hpp"

float Ray_Tracing::deg2rad(float degree)
{
    float pi = 4.0f*atanf(1.0f);

    return degree * pi / 180.0f;
}

void Ray_Tracing::set_specifications()
{
    V = new float[nPoints]();

    // import_binary_float(catch_parameter("vp_model_file", file), V, nPoints);

    for (int index = 0; index < nPoints; index++) V[index] = 1500.0f; 

    std::vector<std::string> splitted;

    // deg2rad()

    splitted = split(catch_parameter("plane_azimuth", file), ',');
    beg_azimuth = deg2rad(std::stof(splitted[0]));
    end_azimuth = deg2rad(std::stof(splitted[1]));
    azimuth_spacing = deg2rad(std::stof(splitted[2]));

    if (end_azimuth < beg_azimuth) 
        throw std::invalid_argument("Error: final plane azimuth angle is smaller than the initial.\n");

    splitted = split(catch_parameter("vertical_angle", file), ',');
    beg_vtangle = deg2rad(std::stof(splitted[0]));
    end_vtangle = deg2rad(std::stof(splitted[1]));    
    vtangle_spacing = deg2rad(std::stof(splitted[2]));

    if (end_azimuth < beg_azimuth) 
        throw std::invalid_argument("Error: final vertical angle is smaller than the initial.\n");

    std::vector<std::string>().swap(splitted);

    nbxl = 1; nbxr = 1;
    nbyl = 1; nbyr = 1;
    nbzu = 1; nbzd = 1;

    nxx = nx + nbxl + nbxr;
    nyy = ny + nbyl + nbyr;
    nzz = nz + nbzu + nbzd;

    volsize = nxx*nyy*nzz;

    S = new float[volsize]();

    expand_boundary(V, S);

    for (int index = 0; index < volsize; index++) S[index] = 1.0f / S[index];

    receiver_output_samples = total_nodes;
    receiver_output = new float[receiver_output_samples]();

    ray = new Ray[receiver_output_samples]();
}

void Ray_Tracing::get_first_arrivals()
{


}

void Ray_Tracing::get_ray_positions()
{

   
}

void Ray_Tracing::build_outputs()
{
// shot_index x y z ---> 4*ray[id]->x.size(); for all receivers

//     wavefield_output_samples = ;

//     wavefield_output = new float[wavefield_output_samples]();

    get_first_arrivals();
    get_ray_positions();
}


// # include <iostream>
// # include <stdlib.h>
// # include <math.h>
// # include <string.h>
// # include "rayTracing.h"

// using namespace std;

// int main(int argc, char **argv)
// {   
//     /* Model parameters */
//     float dx = 1.0;         /* Distance discretization parameter in meter */
//     float dz = 1.0;         /* Depth discretization parameter in meters */
//     int nx = 1000;          /* Total samples in x direction */
//     int nz = 500;           /* Total samples in z direction */

//     /* Ray parameters */
//     float pi = 4*atan(1);   /* Number pi = 3.14159265359*/
//     int nRay = 500;         /* Total samples in ray */
//     float dl = 8.0;         /* Length ray in meters */
//     int angles = 25;        /* Ray number per shot */
//     int edge = 0;           /* Parameter when ray reaches edge */    

//     /* Models */
//     float * vp = (float *) malloc(nx*nz*sizeof(float));    /* Velocities model */
//     float * rays = (float *) malloc(nx*nz*sizeof(float));  /* Matrix to trace rays for visualization */
//     float * s = (float *) malloc(nx*nz*sizeof(float));     /* Slowness model */

//     float * theta = (float *) malloc(angles*sizeof(float));       /* Angles per shot */ 
//     float * travelTimes = (float *) malloc(angles*sizeof(float)); /* Travel times per angles (per ray) */
//     float * sum = (float *) malloc(sizeof(float));                /* Parameter to compute the travel time */

//     Ray * ray = (Ray *) malloc(nRay*sizeof(Ray));    /* Ray object containing (x,z,px,pz,ds) */

//     const char * velocities = "model_smooth_plane_parallel_1000x500.bin"; 
//     importVector2D(vp,nx*nz,velocities);             /* Importing the velocities model */
//     readTextParameter(theta,"angles_NL.txt");        /* Importing angles to trace rays */
//     slowness(s,vp,nx*nz);                            /* Computing the slowness model */

//     // Shot loop

//     setShotPosition(ray,100,1,rays,vp,nx,nz);        /* Initializing the initial position per shot */

//     for (int ang = 0; ang < angles; ang++)
//     {    
//         iConditions(ray,nRay,ang,theta,vp,nx,sum);   /* Initializing PVI equation */

//         for (int path = 0; path < nRay; path++)
//         {        
//             computePx(ray,s,path,nx,dx,dl);          /* Computing x direction component */            
//             computePz(ray,s,path,nx,dz,dl);          /* Computing z direction component */

//             edge = computeX(ray,vp,path,nx,dl);      /* Computing x position point */
//             if(edge == 1) break;                     /* Testing if x is out of bounds */

//             edge = computeZ(ray,vp,path,nx,nz,dl);   /* Computing z position point */
//             if(edge == 1) break;                     /* Testing if z is out of bounds */
    
//             interval(ray,path,vp,nx,dx,dz,sum);      /* Computing distance and time travel per ray lenght */

//             reflectiveLayer(ray,path,vp,nx,6000);    /* Condition to reflect ray for high velocities */
//         }
//         interpolateRay(ray,rays,nx,nRay);            /* To interpolate each ray for visualization */

//         travelTimes[ang] = sum[0];                   /* Catching time travel per ray */
//         timeTravelMessage(ang,travelTimes);          /* Show in window the time calculated */
//     }
//     exportVector2D(rays,nx*nz,"RayModeled.bin");     /* Exporting the rays visualization */

//     return 0;
// }


