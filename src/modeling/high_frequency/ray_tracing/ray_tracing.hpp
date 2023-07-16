# ifndef RAY_TRACING_HPP
# define RAY_TRACING_HPP

# include "../../modeling.hpp"

class Ray_Tracing : public Modeling
{
private:


protected:

    float * V = nullptr;
    float * S = nullptr;

    void get_ray_positions();
    void get_first_arrivals();
    void set_specifications();

    virtual void set_model_boundaries() = 0;
    virtual void set_modeling_message() = 0;
    virtual void set_preconditioners() = 0;

public:

    void build_outputs();

    virtual void initial_setup() = 0;
    virtual void forward_solver() = 0;
    virtual void free_space() = 0;
};

# endif


// # ifndef RAYTRACING_H_INCLUDED
// # define RAYTRACING_H_INCLUDED

// /*  */
// typedef struct
// {
//     int x;
//     int z;
//     float px;
//     float pz;
//     float ds;
// } Ray;

// /* */
// void memSet(float * pointer, int n)
// {
//     for(int i = 0; i < n; i++) pointer[n] = 0.0f;
// }

// /* */
// void slowness(float * s,float * vp, int nPoints)
// {
//     for (int point = 0; point < nPoints; point++)
//         s[point] = 1/vp[point];
// }

// /* */
// void importVector2D(float * vector, int nPoints, const char * filename)
// {
//     FILE * read = fopen(filename,"rb");
//     fread(vector,sizeof(float),nPoints,read);
//     fclose(read);
// }

// /* */
// void readTextParameter(float * vector, const char * filename)
// {
//     FILE * rf = fopen(filename, "r");
//     if(rf != NULL) 
//     {
//         int ii = 0;
//         while(fscanf(rf,"%f",&vector[ii]) != EOF) 
//         {
//             ++ii;                    
//         }
//     } 
//     fclose(rf);
// }   

// /* */
// void exportVector2D(float * vector, int nPoints, const char * filename)
// {
//     FILE * write = fopen(filename, "wb");
//     fwrite(vector, sizeof(float), nPoints, write);
//     fclose(write);
// }

// /* */
// void setShotPosition(Ray * r, int x_i, int z_i, float * matRays,float * vp, int nx, int nz)
// {
//     r[0].x = x_i;
//     r[0].z = z_i;

//     for(int point = 0; point < nx*nz; point++) 
//         matRays[point] = vp[point] * 1e-4;    

//     matRays[r[0].z * nx + r[0].x] = 1.0f;
// }

// /* */
// void iConditions(Ray * r, int nRay, int actualAngle, float * theta, float * v, int nx, float * sum)
// {
//     int index;
//     for (index = 1; index < nRay; index++)
//     {
//         r[index].x = 0.0f;
//         r[index].z = 0.0f;
//         r[index].px = 0.0f;
//         r[index].pz = 0.0f;
//         r[index].ds = 0.0f;
//     }

//     r[0].px = cos(theta[actualAngle]) / v[r[0].z * nx + r[0].x]; 
//     r[0].pz = sin(theta[actualAngle]) / v[r[0].z * nx + r[0].x];
//     r[0].ds = 0.0f;
//     sum[0] = 0.0f;
// }

// /* */
// void computePx(Ray * r, float * s, int path, int nx, float dx, float dl)
// {
//     float sx_a = (s[r[path].z * nx + (r[path].x + 1)] - s[r[path].z * nx + (r[path].x - 1)])/(2*dx); 
//     float sx_half = (s[r[path].z * nx + (r[path].x + 1)] - s[r[path].z * nx + r[path].x])/dx;
//     float sx_b = (s[r[path].z * nx + (r[path].x + 2)] - s[r[path].z * nx + r[path].x])/(2*dx);

//     float k1 = sx_a;
//     float k2 = sx_half;
//     float k3 = sx_half;
//     float k4 = sx_b;
//     float k = (k1 + 2*k2 + 2*k3 + k4) / 6.0f;
    
//     r[path+1].px = r[path].px + k * dl;
// }

// /* */
// void computePz(Ray * r, float * s, int path, int nx, float dz, float dl)
// {
//     float sz_a = (s[(r[path].z + 1) * nx + r[path].x] - s[(r[path].z - 1) * nx + r[path].x])/(2*dz);
//     float sz_half = (s[(r[path].z + 1) * nx + r[path].x] - s[r[path].z * nx + r[path].x])/dz;
//     float sz_b = (s[(r[path].z + 2) * nx + r[path].x] - s[r[path].z * nx + r[path].x])/(2*dz);

//     float k1 = sz_a;
//     float k2 = sz_half;
//     float k3 = sz_half;
//     float k4 = sz_b;
//     float k = (k1 + 2*k2 + 2*k3 + k4) / 6.0f;

//     r[path+1].pz = r[path].pz + k * dl;
// }

// /* */
// int computeX(Ray * r, float * vp, int path, int nx, float dl)
// {
//     float vp_int_x = (vp[r[path].z * nx + (r[path].x + 1)] + vp[r[path].z * nx + r[path].x]) / 2.0f; 
//     float px_int = (r[path+1].px + r[path].px) / 2.0f;

//     float k1 = r[path].px * vp[r[path].z * nx + r[path].x];
//     float k2 = px_int * vp_int_x;
//     float k3 = px_int * vp_int_x;
//     float k4 = r[path+1].px * vp[r[path].z * nx + (r[path].x + 1)];
//     float k = (k1 + 2*k2 + 2*k3 + k4) / 6.0f;

//     r[path+1].x = r[path].x + k * dl;

//     if ((r[path+1].x <= 0) || (r[path+1].x >= nx-1)) 
//         return 1;
// }

// /* */
// int computeZ(Ray * r, float * vp, int path, int nx,int nz, float dl)
// {
//     float vp_int_z = (vp[(r[path].z + 1) * nx + r[path].x] + vp[r[path].z * nx + r[path].x]) / 2.0f;
//     float pz_int = (r[path+1].pz + r[path].pz) / 2.0f;

//     float k1 = r[path].pz * vp[r[path].z * nx + r[path].x];
//     float k2 = vp_int_z * pz_int;
//     float k3 = vp_int_z * pz_int;
//     float k4 = r[path+1].pz * vp[(r[path].z + 1) * nx + r[path].x];
//     float k = (k1 + 2*k2 + 2*k3 + k4) / 6.0f;

//     r[path+1].z = r[path].z + k * dl;

//     if ((r[path+1].z <= 0) || (r[path+1].z >= nz-1)) 
//         return 1;        
// }

// /* */
// void interval(Ray * r, int path, float * vp, int nx, float dx, float dz, float * sum)
// {
//     r[path+1].ds = sqrtf(powf(dx*(r[path+1].x - r[path].x),2.0f) + powf(dz*(r[path+1].z - r[path].z),2.0f));
//     sum[0] += r[path+1].ds / vp[r[path+1].z * nx + r[path+1].x];
// }

// /* */
// float lagrangeInterpolator(int point, Ray * r, int init, int final)
// {
//     int j, k;
//     float prod1, prod2;
//     float P, L;

//     P = 0.0;
//     for(k = init; k < final; k++)
//     {
//         L = 1.0; prod1 = 1; prod2 = 1;
//         for(j = init; j < final; j++)
//         {
//             if(k != j)
//             {
//                 prod1 *= (point - r[j].x);
//                 prod2 *= (r[k].x - r[j].x);
//             }            
//         }
//         L = prod1 / prod2;
//         P += r[k].z * L;
//     }
//     return P;
// }

// /* */
// void interpolateRay(Ray * r, float * matRays, int nx, int nRay)
// {
//     int len = 0;
//     for (int i = 0; i < nRay; i++)
//     {
//         if(r[i].x == 0) break;
//         ++len;
//     }

//     int total = r[len-1].x - r[0].x + 1;
//     int * X = (int *) malloc(total*sizeof(int));
//     int * Z = (int *) malloc(total*sizeof(int));
    
//     int point;
//     int par = 0;
//     for(int i = 0; i < total; i++)
//     {
//         if(r[par].x == r[0].x + i)
//         {
//             X[i] = r[par].x;
//             Z[i] = r[par].z; 
//             par++;
//         }
//         else
//         {
//             X[i] = r[0].x + i;
//             Z[i] = (int) lagrangeInterpolator(X[i],r,par-1,par+1);
//         }
//         matRays[Z[i] * nx + X[i]] = 1.0f;
//     }
// }

// /* */ 
// void reflectiveLayer(Ray * r, int path, float * vp, int nx, float refVel)
// {
//     if (vp[r[path+1].z * nx + r[path+1].x] > refVel)
//         r[path+1].pz = -r[path+1].pz;
// }

// /* */
// void timeTravelMessage(int ang, float * travelTimes)
// {
//     std::cout <<"Raio "<<ang+1<<": tempo = "<<travelTimes[ang]<<std::endl;
//     std::cout <<"---------------------"<<std::endl;
// }

// # endif