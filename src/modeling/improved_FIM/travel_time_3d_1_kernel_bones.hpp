#pragma once
#ifndef _TRAVEL_TIME_3D_KERNEL_BONES
#define _TRAVEL_TIME_3D_KERNEL_BONES
#include "travel_time_3d_0_bones.hpp"
namespace jarvis
{
    __device__ bool is_adjacent_node_converging(int tx, int ty, int tz, int idx, Frame vel_grid, NodeStatus *_mark);
    __device__ double diff3d_1(int ix, int iy, int iz, Frame vel_grid, float *_vel, double *_time);
    __device__ double diff3d_1diag(int ix, int iy, int iz, Frame vel_grid, float *_vel, double *_time);
    __device__ double diff3d_2(int ix, int iy, int iz, Frame vel_grid, float *_vel, double *_time);
    __device__ double diff3d_2diag(int ix, int iy, int iz, Frame vel_grid, float *_vel, double *_time);
    //
    __global__ void fim_kernel_0_cal_active_node(Frame vel_grid, float *_vel, double *_time, NodeStatus *_mark, int diff_order);
    __global__ void fim_kernel_1_mark_new_active_node(Frame vel_grid, float *_vel, double *_time, NodeStatus *_mark, int diff_order);
    __global__ void fim_kernel_2_set_as_converged_node(Frame vel_grid, NodeStatus *_mark);
    __global__ void fim_kernel_3_check_is_finishd(Frame vel_grid, NodeStatus *_mark, bool *endflag);
}
#endif