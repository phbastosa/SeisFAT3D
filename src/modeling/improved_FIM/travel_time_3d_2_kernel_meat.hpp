#pragma once
#ifndef _TRAVEL_TIME_3D_KERNEL_MEAT
#define _TRAVEL_TIME_3D_KERNEL_MEAT
#include "travel_time_3d_1_kernel_bones.hpp"
namespace jarvis
{
#define F_(a, i, j, k) (a)[(i) + (j)*n_rows + (k)*n_elem_slice]

#define SWAP(type, a, b) \
	{                    \
		type temp;       \
		temp = a;        \
		a = b;           \
		b = temp;        \
	}

	inline __device__ bool CudaNextPermutation(int *p, int n)
	{
		int last = n - 1;
		int i, j, k;
		i = last;
		while (i > 0 && p[i] < p[i - 1])
			i--;
		if (i == 0)
			return false;
		k = i;
		for (j = last; j >= i; j--)
			if (p[j] > p[i - 1] && p[j] < p[k])
				k = j;
		SWAP(int, p[k], p[i - 1]);
		for (j = last, k = i; j > k; j--, k++)
			SWAP(int, p[j], p[k]);
		return true;
	}

	inline __global__ void fim_kernel_0_cal_active_node(Frame vel_grid, float *_vel, double *_time, NodeStatus *_mark, int diff_order)
	{ // If this node is the active node, travel time is calculated
		set_cufield_3d_idx(vel_grid, idx, tx, ty, tz);
		double old_trav = LONG_MAX;
		double trav = LONG_MAX;
		if (idx < vel_grid.n_elem)
		{
			if (_mark[idx] == NodeStatus::active)
			{
				old_trav = _time[idx];
				if (diff_order == 11)
					trav = diff3d_1(tx, ty, tz, vel_grid, _vel, _time);
				if (diff_order == 12)
					trav = diff3d_1diag(tx, ty, tz, vel_grid, _vel, _time);
				if (diff_order == 21)
					trav = diff3d_2(tx, ty, tz, vel_grid, _vel, _time);
				if (diff_order == 22)
					trav = diff3d_2diag(tx, ty, tz, vel_grid, _vel, _time);
				//
				_time[idx] = trav;
				_mark[idx] = NodeStatus::converging; // Marked as converging node
			}
		}
	}

	inline __global__ void fim_kernel_1_mark_new_active_node(Frame vel_grid, float *_vel, double *_time, NodeStatus *_mark, int diff_order)
	{ // If the node is not active and has not converged, and the nearby node is to be converged, the node is calculated and marked as the active node
		set_cufield_3d_idx(vel_grid, idx, tx, ty, tz);
		double old_trav = UINT_MAX;
		double trav = UINT_MAX;
		if (idx < vel_grid.n_elem)
		{
			if (_mark[idx] != NodeStatus::converging && _mark[idx] != NodeStatus::not_active)
			{
				if (is_adjacent_node_converging(tx, ty, tz, idx, vel_grid, _mark) && _mark[idx] != NodeStatus::active)
				{
					old_trav = _time[idx];
					if (diff_order == 11)
						trav = diff3d_1(tx, ty, tz, vel_grid, _vel, _time);

					if (diff_order == 12)
						trav = diff3d_1diag(tx, ty, tz, vel_grid, _vel, _time);

					if (diff_order == 21)
						trav = diff3d_2(tx, ty, tz, vel_grid, _vel, _time);

					if (diff_order == 22)
						trav = diff3d_2diag(tx, ty, tz, vel_grid, _vel, _time);

					if (old_trav > trav)
					{
						_time[idx] = trav;
						_mark[idx] = NodeStatus::active;
					}
				}
			}
		}
	}

	inline __global__ void fim_kernel_2_set_as_converged_node(Frame vel_grid, NodeStatus *_mark)
	{ // Temporarily converted to converged node, which may become active node later in function "fim_kernel_1_mark_new_active_node"
		set_cufield_3d_idx(vel_grid, idx, tx, ty, tz);
		if (idx < vel_grid.n_elem)
		{
			if (_mark[idx] == NodeStatus::converging)
			{
				_mark[idx] = NodeStatus::converged;
			}
		}
	}

	inline __global__ void fim_kernel_3_check_is_finishd(Frame vel_grid, NodeStatus *_mark, bool *endflag)
	{ // Do not end as long as there is an active node
		set_cufield_3d_idx(vel_grid, idx, tx, ty, tz);
		if (idx < vel_grid.n_elem)
		{
			if (_mark[idx] == NodeStatus::active)
			{
				endflag[0] = false;
			}
		}
	}

	inline __device__ bool is_adjacent_node_converging(int tx, int ty, int tz, int idx, Frame vel_grid, NodeStatus *_mark)
	{
		set_frame_n(vel_grid);
		int lax = -1;
		int rax = 1;
		int lay = -1;
		int ray = 1;
		int laz = -1;
		int raz = 1;
		if (tx <= 0)
			lax = 1;
		if (tx >= n_rows - 1)
			rax = -1;
		if (ty <= 0)
			lay = 1;
		if (ty >= n_cols - 1)
			ray = -1;
		if (tz <= 0)
			laz = 1;
		if (tz >= n_slices - 1)
			raz = -1;
		if ((_mark[idx + lax] == NodeStatus::converging) || (_mark[idx + rax] == NodeStatus::converging) ||
			(_mark[idx + lay * n_rows] == NodeStatus::converging) || (_mark[idx + ray * n_rows] == NodeStatus::converging) ||
			(_mark[idx + laz * n_elem_slice] == NodeStatus::converging) || (_mark[idx + raz * n_elem_slice] == NodeStatus::converging))
		{
			return true;
		}
		else
		{
			return false;
		}
	}

	inline __device__ double diff3d_1(int ix, int iy, int iz, Frame vel_grid, float *_vel, double *_time)
	{
		set_frame_n(vel_grid);
		double d_rows = vel_grid.d_rows;
		double a, b, c, slown;
		double trav, travm;
		int xtmp, ytmp, ztmp;
		double ttnminx, ttnminy, ttnminz;
		int ax, ay, az;
		trav = UINT_MAX, travm = UINT_MAX;
		xtmp = ix;
		ytmp = iy;
		ztmp = iz;
		if (ix == n_rows)
			xtmp = n_rows - 1;
		if (iy == n_cols)
			ytmp = n_cols - 1;
		if (iz == n_slices)
			ztmp = n_slices - 1;
		slown = 1.0f / _vel[xtmp + ytmp * n_rows + ztmp * n_elem_slice];

		if ((ix >= 0 && ix < n_rows) && (iy >= 0 && iy < n_cols) && (iz >= 0 && iz < n_slices))
		{
			if (ix - 1 < 0)
				ax = 1;
			else if (ix + 1 >= n_rows)
				ax = -1;
			else if (F_(_time, (ix + 1), iy, iz) <= F_(_time, (ix - 1), iy, iz))
				ax = 1;
			else
				ax = -1;

			if (iy - 1 < 0)
				ay = 1;
			else if (iy + 1 >= n_cols)
				ay = -1;
			else if (F_(_time, ix, (iy + 1), iz) <= F_(_time, ix, (iy - 1), iz))
				ay = 1;
			else
				ay = -1;

			if (iz - 1 < 0)
				az = 1;
			else if (iz + 1 >= n_slices)
				az = -1;
			else if (F_(_time, ix, iy, (iz + 1)) <= F_(_time, ix, iy, (iz - 1)))
				az = 1;
			else
				az = -1;

			ttnminx = F_(_time, (ix + ax), iy, iz);
			ttnminy = F_(_time, ix, (iy + ay), iz);
			ttnminz = F_(_time, ix, iy, (iz + az));

			int p[3] = {0, 1, 2};
			double ttnmin[3] = {ttnminx, ttnminy, ttnminz};
			do
			{
				if (ttnmin[p[0]] <= ttnmin[p[1]] && ttnmin[p[1]] <= ttnmin[p[2]])
				{
					c = ttnmin[p[0]];
					b = ttnmin[p[1]];
					a = ttnmin[p[2]];
				}
			} while (CudaNextPermutation(p, 3));
			trav = c + d_rows * slown;
			if (trav > b)
			{
				trav = (b + c +
						sqrt(-b * b - c * c + 2.0f * b * c +
							 2.0f * (d_rows * d_rows * slown * slown))) /
					   2.0f;
				if (trav > a)
				{
					trav =
						(2.0f * (a + b + c) + sqrt(4.0f * (a + b + c) * (a + b + c) -
												   12.0f * (a * a + b * b + c * c -
															(d_rows * d_rows * slown * slown)))) /
						6.0f;
				}
			}
			travm = trav;
		}
		return travm;
	}

	inline __device__ double diff3d_1diag(int ix, int iy, int iz, Frame vel_grid, float *_vel, double *_time)
	{
		set_frame_n(vel_grid);
		double d_rows = vel_grid.d_rows;
		double a, b, c;
		double slown;
		double trav, travm;
		int xtmp, ytmp, ztmp;
		double ttnminx, ttnminy, ttnminz;
		int ax, ay, az;
		double trav_diag;
		trav = UINT_MAX, travm = UINT_MAX, trav_diag = UINT_MAX;

		xtmp = ix;
		ytmp = iy;
		ztmp = iz;
		if (ix == n_rows)
			xtmp = n_rows - 1;
		if (iy == n_cols)
			ytmp = n_cols - 1;
		if (iz == n_slices)
			ztmp = n_slices - 1;

		slown = 1.0f / F_(_vel, xtmp, ytmp, ztmp);

		if ((ix >= 0 && ix < n_rows) && (iy >= 0 && iy < n_cols) && (iz >= 0 && iz < n_slices))
		{
			if (ix - 1 < 0)
				ax = 1;
			else if (ix + 1 >= n_rows)
				ax = -1;
			else if (F_(_time, (ix + 1), iy, iz) <= F_(_time, (ix - 1), iy, iz))
				ax = 1;
			else
				ax = -1;

			if (iy - 1 < 0)
				ay = 1;
			else if (iy + 1 >= n_cols)
				ay = -1;
			else if (F_(_time, ix, (iy + 1), iz) <= F_(_time, ix, (iy - 1), iz))
				ay = 1;
			else
				ay = -1;

			if (iz - 1 < 0)
				az = 1;
			else if (iz + 1 >= n_slices)
				az = -1;
			else if (F_(_time, ix, iy, (iz + 1)) <= F_(_time, ix, iy, (iz - 1)))
				az = 1;
			else
				az = -1;

			ttnminx = F_(_time, (ix + ax), iy, iz);
			ttnminy = F_(_time, ix, (iy + ay), iz);
			ttnminz = F_(_time, ix, iy, (iz + az));

			int p[3] = {0, 1, 2}; //????123?????
			double ttnmin[3] = {ttnminx, ttnminy, ttnminz};
			do
			{
				if (ttnmin[p[0]] <= ttnmin[p[1]] && ttnmin[p[1]] <= ttnmin[p[2]])
				{
					c = ttnmin[p[0]];
					b = ttnmin[p[1]];
					a = ttnmin[p[2]];
				}
			} while (CudaNextPermutation(p, 3));

			trav = c + d_rows * slown;

			if (trav > b)
			{
				trav = (b + c +
						sqrt(-b * b - c * c + 2.0f * b * c +
							 2.0f * (d_rows * d_rows * slown * slown))) /
					   2.0f;

				if (trav > a)
				{
					trav =
						(2.0f * (a + b + c) + sqrt(4.0f * (a + b + c) * (a + b + c) -
												   12.0f * (a * a + b * b + c * c -
															(d_rows * d_rows * slown * slown)))) /
						6.0f;

					trav_diag =
						F_(_time, (ix + ax), (iy + ay), (iz + az)) + sqrt(3.0f) * d_rows * slown;

					if (trav_diag < trav)
					{
						trav = trav_diag;
					}
				}
			}
			travm = trav;
		}
		return travm;
	}

	inline __device__ double diff3d_2(int ix, int iy, int iz, Frame vel_grid, float *_vel, double *_time)
	{
		set_frame_n(vel_grid);
		double d_rows = vel_grid.d_rows;
		double d_cols = vel_grid.d_cols;
		double d_slices = vel_grid.d_slices;
		double a, b, c, slown;
		double a2, b2, c2;
		double aa, bb, cc;
		double aa1, bb1, cc1;
		double aa2, bb2, cc2;
		double aa3, bb3, cc3;
		double trav, trav1, trav2, trav3, trav4, trav5, trav6, travm, time_grad_vec_length;
		int xtmp, ytmp, ztmp;
		double ttnminx, ttnminy, ttnminz;
		double ttnminx2, ttnminy2, ttnminz2;
		ttnminx2 = UINT_MAX, ttnminy2 = UINT_MAX, ttnminz2 = UINT_MAX;
		int ax, ay, az;

		trav = 0.0f, trav1 = UINT_MAX, trav2 = UINT_MAX, trav3 = UINT_MAX,
		trav4 = UINT_MAX, trav5 = UINT_MAX, trav6 = UINT_MAX, travm = UINT_MAX, time_grad_vec_length = 0.0f;
		xtmp = ix;
		ytmp = iy;
		ztmp = iz;
		if (ix == n_rows)
			xtmp = n_rows - 1;
		if (iy == n_cols)
			ytmp = n_cols - 1;
		if (iz == n_slices)
			ztmp = n_slices - 1;
		slown = 1.0f / F_(_vel, xtmp, ytmp, ztmp);
		if ((ix >= 0 && ix < n_rows) && (iy >= 0 && iy < n_cols) && (iz >= 0 && iz < n_slices))
		{
			if (ix - 1 < 0)
				ax = 1;
			else if (ix + 1 >= n_rows)
				ax = -1;
			else if (F_(_time, (ix + 1), iy, iz) < F_(_time, (ix - 1), iy, iz))
				ax = 1;
			else
				ax = -1;

			if (iy - 1 < 0)
				ay = 1;
			else if (iy + 1 >= n_cols)
				ay = -1;
			else if (F_(_time, ix, (iy + 1), iz) < F_(_time, ix, (iy - 1), iz))
				ay = 1;
			else
				ay = -1;

			if (iz - 1 < 0)
				az = 1;
			else if (iz + 1 >= n_slices)
				az = -1;
			else if (F_(_time, ix, iy, (iz + 1)) < F_(_time, ix, iy, (iz - 1)))
				az = 1;
			else
				az = -1;

			ttnminx = F_(_time, (ix + ax), iy, iz);
			ttnminy = F_(_time, ix, (iy + ay), iz);
			ttnminz = F_(_time, ix, iy, (iz + az));

			ttnminx2 = F_(_time, (ix + 2 * ax), iy, iz);
			ttnminy2 = F_(_time, ix, (iy + 2 * ay), iz);
			ttnminz2 = F_(_time, ix, iy, (iz + 2 * az));

			int p[3] = {0, 1, 2};
			double ttnmin[3] = {ttnminx, ttnminy, ttnminz};
			double ttnmin2[3] = {ttnminx2, ttnminy2, ttnminz2};
			do
			{
				if (ttnmin[p[0]] <= ttnmin[p[1]] && ttnmin[p[1]] <= ttnmin[p[2]])
				{
					c = ttnmin[p[0]];
					b = ttnmin[p[1]];
					a = ttnmin[p[2]];
					c2 = ttnmin2[p[0]];
					b2 = ttnmin2[p[1]];
					a2 = ttnmin2[p[2]];
				}
			} while (CudaNextPermutation(p, 3));

			trav = c + d_rows * slown;

			if (trav > b)
			{
				trav1 = (b + c +
						 sqrt(-b * b - c * c + 2.0f * b * c +
							  2.0f * d_rows * d_rows * slown * slown)) /
						2.0f;
				if (trav1 < trav)
				{
					trav = trav1;
				}
				if (c2 < c)
				{
					trav2 =
						(3.0f * (4.0f * c - c2) + 4.0f * b +
						 2.0f * sqrt(13 * slown * slown * d_rows * d_rows -
									 (4.0f * c - c2 - 3.0f * b) * (4.0f * c - c2 - 3.0f * b))) /
						13.0f;
					if (trav2 < trav)
					{
						trav = trav2;
					}
					if (b2 < b)
					{
						trav3 = ((4.0f * c - c2) + (4.0f * b - b2) +
								 sqrt(8.0f * (slown * slown * d_rows * d_rows) -
									  ((4.0f * c - c2) - (4.0f * b - b2)) *
										  ((4.0f * c - c2) - (4.0f * b - b2)))) /
								6.0f;
					}
					if (trav3 < trav)
					{
						trav = trav3;
					}
				}
				if (trav > a)
				{
					trav =
						(2.0f * (a + b + c) + sqrt(4.0f * (a + b + c) * (a + b + c) -
												   12.0f * (a * a + b * b + c * c -
															d_rows * d_rows * (slown * slown)))) /
						6.0f;

					if (c2 < c)
					{
						aa1 = 9.0f / (4.0f * d_rows * d_rows);
						bb1 = (-24.0f * c + 6.0f * c2) / (4.0f * d_rows * d_rows);
						cc1 = ((4.0f * c - c2) * (4.0f * c - c2)) / (4.0f * d_rows * d_rows);

						aa2 = 1.0f / (d_cols * d_cols);
						bb2 = -2.0f * b / (d_cols * d_cols);
						cc2 = b * b / (d_cols * d_cols);

						aa3 = 1.0f / (d_slices * d_slices);
						bb3 = -2.0f * a / (d_slices * d_slices);
						cc3 = a * a / (d_slices * d_slices);

						aa = aa1 + aa2 + aa3;
						bb = bb1 + bb2 + bb3;
						cc = cc1 + cc2 + cc3 - (slown * slown);

						time_grad_vec_length = bb * bb - 4.0f * aa * cc;
						if (time_grad_vec_length < 0.0f)
							time_grad_vec_length = 0.0f;
						trav4 = (-bb + sqrt(time_grad_vec_length)) / (2.0f * aa);
						if (trav4 < trav)
						{
							trav = trav4;
						}

						if (b2 < b)
						{
							aa1 = 9.0f / (4.0f * d_rows * d_rows);
							bb1 = (-24.0f * c + 6.0f * c2) / (4.0f * d_rows * d_rows);
							cc1 = ((4.0f * c - c2) * (4.0f * c - c2)) / (4.0f * d_rows * d_rows);
							aa2 = 9.0f / (4.0f * d_cols * d_cols);
							bb2 = (-24.0f * b + 6.0f * b2) / (4.0f * d_cols * d_cols);
							cc2 = ((4.0f * b - b2) * (4.0f * b - b2)) / (4.0f * d_cols * d_cols);

							aa3 = 1.0f / (d_slices * d_slices);
							bb3 = -2.0f * a / (d_slices * d_slices);
							cc3 = a * a / (d_slices * d_slices);

							aa = aa1 + aa2 + aa3;
							bb = bb1 + bb2 + bb3;
							cc = cc1 + cc2 + cc3 - (slown * slown);
							time_grad_vec_length = bb * bb - 4.0f * aa * cc;
							if (time_grad_vec_length < 0.0f)
								time_grad_vec_length = 0.0f;
							trav5 = (-bb + sqrt(time_grad_vec_length)) / (2.0f * aa);
							if (trav5 < trav)
							{
								trav = trav5;
							}
							if (a2 < a)
							{
								aa1 = 9.0f / (4.0f * d_rows * d_rows);
								bb1 = (-24.0f * c + 6.0f * c2) / (4.0f * d_rows * d_rows);
								cc1 = ((4.0f * c - c2) * (4.0f * c - c2)) / (4.0f * d_rows * d_rows);
								//
								aa2 = 9.0f / (4.0f * d_cols * d_cols);
								bb2 = (-24.0f * b + 6.0f * b2) / (4.0f * d_cols * d_cols);
								cc2 = ((4.0f * b - b2) * (4.0f * b - b2)) / (4.0f * d_cols * d_cols);
								//
								aa3 = 9.0f / (4.0f * d_slices * d_slices);
								bb3 = (-24.0f * a + 6.0f * a2) / (4.0f * d_slices * d_slices);
								cc3 = ((4.0f * a - a2) * (4.0f * a - a2)) / (4.0f * d_slices * d_slices);
								//
								aa = aa1 + aa2 + aa3;
								bb = bb1 + bb2 + bb3;
								cc = cc1 + cc2 + cc3 - (slown * slown);
								time_grad_vec_length = bb * bb - 4.0f * aa * cc;
								if (time_grad_vec_length < 0.0f)
									time_grad_vec_length = 0.0f;
								trav6 = (-bb + sqrt(time_grad_vec_length)) / (2.0f * aa);
								if (trav6 < trav)
								{
									trav = trav6;
								}
							}
						}
					}
				}
			}
			travm = trav;
		}
		return travm;
	}

	inline __device__ double diff3d_2diag(int ix, int iy, int iz, Frame vel_grid, float *_vel, double *_time)
	{
		set_frame_nd(vel_grid);
		double a, b, c, slown;
		double a2, b2, c2;
		double aa, bb, cc;
		double aa1, bb1, cc1;
		double aa2, bb2, cc2;
		double aa3, bb3, cc3;
		double trav, trav1, trav2, trav3, trav4, trav5, trav6, trav7, trav_d, travm, time_grad_vec_length;
		int xtmp, ytmp, ztmp;
		double ttnminx, ttnminy, ttnminz;
		double ttnminx2, ttnminy2, ttnminz2;
		ttnminx2 = UINT_MAX, ttnminy2 = UINT_MAX, ttnminz2 = UINT_MAX;
		int ax, ay, az;

		double trav_diag = UINT_MAX;
		trav = 0.0f, trav1 = UINT_MAX, trav2 = UINT_MAX, trav3 = UINT_MAX,
		trav4 = UINT_MAX, trav5 = UINT_MAX, trav6 = UINT_MAX, trav7 = UINT_MAX, trav_d = UINT_MAX,
		travm = UINT_MAX, time_grad_vec_length = 0.0f;
		xtmp = ix;
		ytmp = iy;
		ztmp = iz;
		if (ix == n_rows)
			xtmp = n_rows - 1;
		if (iy == n_cols)
			ytmp = n_cols - 1;
		if (iz == n_slices)
			ztmp = n_slices - 1;
		slown = 1.0f / F_(_vel, xtmp, ytmp, ztmp);
		if ((ix >= 0 && ix < n_rows) && (iy >= 0 && iy < n_cols) && (iz >= 0 && iz < n_slices))
		{
			if (ix - 1 < 0)
				ax = 1;
			else if (ix + 1 >= n_rows)
				ax = -1;
			else if (F_(_time, (ix + 1), iy, iz) <= F_(_time, (ix - 1), iy, iz))
				ax = 1;
			else
				ax = -1;

			if (iy - 1 < 0)
				ay = 1;
			else if (iy + 1 >= n_cols)
				ay = -1;
			else if (F_(_time, ix, (iy + 1), iz) <= F_(_time, ix, (iy - 1), iz))
				ay = 1;
			else
				ay = -1;

			if (iz - 1 < 0)
				az = 1;
			else if (iz + 1 >= n_slices)
				az = -1;
			else if (F_(_time, ix, iy, (iz + 1)) <= F_(_time, ix, iy, (iz - 1)))
				az = 1;
			else
				az = -1;

			ttnminx = F_(_time, (ix + ax), iy, iz);
			ttnminy = F_(_time, ix, (iy + ay), iz);
			ttnminz = F_(_time, ix, iy, (iz + az));

			ttnminx2 = F_(_time, (ix + 2 * ax), iy, iz);
			ttnminy2 = F_(_time, ix, (iy + 2 * ay), iz);
			ttnminz2 = F_(_time, ix, iy, (iz + 2 * az));

			int p[3] = {0, 1, 2}; //????123?????
			double ttnmin[3] = {ttnminx, ttnminy, ttnminz};
			double ttnmin2[3] = {ttnminx2, ttnminy2, ttnminz2};
			do
			{
				if (ttnmin[p[0]] <= ttnmin[p[1]] && ttnmin[p[1]] <= ttnmin[p[2]])
				{
					c = ttnmin[p[0]];
					b = ttnmin[p[1]];
					a = ttnmin[p[2]];
					c2 = ttnmin2[p[0]];
					b2 = ttnmin2[p[1]];
					a2 = ttnmin2[p[2]];
				}
			} while (CudaNextPermutation(p, 3));

			trav = c + d_rows * slown;

			if (trav > b)
			{
				trav = (b + c +
						sqrt(-b * b - c * c + 2.0f * b * c +
							 2.0f * d_rows * d_rows * (slown * slown))) /
					   2.0f;

				if (c2 < c)
				{
					trav1 =
						(3.0f * (4.0f * c - c2) + 4.0f * b +
						 2.0f * sqrt(13 * slown * slown * d_rows * d_rows -
									 (4.0f * c - c2 - 3.0f * b) * (4.0f * c - c2 - 3.0f * b))) /
						13.0f;

					if (b2 < b)
					{
						trav2 = ((4.0f * c - c2) + (4.0f * b - b2) +
								 sqrt(8.0f * (slown * slown * d_rows * d_rows) -
									  ((4.0f * c - c2) - (4.0f * b - b2)) * ((4.0f * c - c2) - (4.0f * b - b2)))) /
								6.0f;
					}
				}
				//
				if (trav1 < trav)
				{
					trav = trav1;
				}
				if (trav2 < trav)
				{
					trav = trav2;
				}
				//
				if (trav > a)
				{
					trav =
						(2.0f * (a + b + c) + sqrt(4.0f * (a + b + c) * (a + b + c) -
												   12.0f * (a * a + b * b + c * c -
															d_rows * d_rows * (slown * slown)))) /
						6.0f;

					if (c2 < c)
					{
						aa1 = 9.0f / (4.0f * d_rows * d_rows);
						bb1 = (-24.0f * c + 6.0f * c2) / (4.0f * d_rows * d_rows);
						cc1 = ((4.0f * c - c2) * (4.0f * c - c2)) / (4.0f * d_rows * d_rows);

						aa2 = 1.0f / (d_cols * d_cols);
						bb2 = -2.0f * b / (d_cols * d_cols);
						cc2 = b * b / (d_cols * d_cols);

						aa3 = 1.0f / (d_slices * d_slices);
						bb3 = -2.0f * a / (d_slices * d_slices);
						cc3 = a * a / (d_slices * d_slices);

						aa = aa1 + aa2 + aa3;
						bb = bb1 + bb2 + bb3;
						cc = cc1 + cc2 + cc3 - (slown * slown);

						time_grad_vec_length = bb * bb - 4.0f * aa * cc;
						if (time_grad_vec_length < 0.0f)
							time_grad_vec_length = 0.0f;
						trav3 = (-bb + sqrt(time_grad_vec_length)) / (2.0f * aa);

						if (b2 < b)
						{
							aa1 = 9.0f / (4.0f * d_rows * d_rows);
							bb1 = (-24.0f * c + 6.0f * c2) / (4.0f * d_rows * d_rows);
							cc1 = ((4.0f * c - c2) * (4.0f * c - c2)) / (4.0f * d_rows * d_rows);
							aa2 = 9.0f / (4.0f * d_cols * d_cols);
							bb2 = (-24.0f * b + 6.0f * b2) / (4.0f * d_cols * d_cols);
							cc2 = ((4.0f * b - b2) * (4.0f * b - b2)) / (4.0f * d_cols * d_cols);

							aa3 = 1.0f / (d_slices * d_slices);
							bb3 = -2.0f * a / (d_slices * d_slices);
							cc3 = a * a / (d_slices * d_slices);

							aa = aa1 + aa2 + aa3;
							bb = bb1 + bb2 + bb3;
							cc = cc1 + cc2 + cc3 - (slown * slown);
							time_grad_vec_length = bb * bb - 4.0f * aa * cc;
							if (time_grad_vec_length < 0.0f)
								time_grad_vec_length = 0.0f;
							trav4 = (-bb + sqrt(time_grad_vec_length)) / (2.0f * aa);

							if (a2 < a)
							{
								aa1 = 9.0f / (4.0f * d_rows * d_rows);
								bb1 = (-24.0f * c + 6.0f * c2) / (4.0f * d_rows * d_rows);
								cc1 = ((4.0f * c - c2) * (4.0f * c - c2)) / (4.0f * d_rows * d_rows);
								aa2 = 9.0f / (4.0f * d_cols * d_cols);
								bb2 = (-24.0f * b + 6.0f * b2) / (4.0f * d_cols * d_cols);
								cc2 = ((4.0f * b - b2) * (4.0f * b - b2)) / (4.0f * d_cols * d_cols);
								aa3 = 9.0f / (4.0f * d_slices * d_slices);
								bb3 = (-24.0f * a + 6.0f * a2) / (4.0f * d_slices * d_slices);
								cc3 = ((4.0f * a - a2) * (4.0f * a - a2)) / (4.0f * d_slices * d_slices);
								aa = aa1 + aa2 + aa3;
								bb = bb1 + bb2 + bb3;
								cc = cc1 + cc2 + cc3 - (slown * slown);
								time_grad_vec_length = bb * bb - 4.0f * aa * cc;
								if (time_grad_vec_length < 0.0f)
									time_grad_vec_length = 0.0f;
								trav5 = (-bb + sqrt(time_grad_vec_length)) / (2.0f * aa);
								//
								trav_diag = F_(_time, ix + ax, iy + ay, iz + az) + sqrt(3.0f) * d_rows * slown;
							}
						}
					}
					//
					if (trav_diag < trav)
					{
						trav = trav_diag;
					}
					if (trav3 < trav)
					{
						trav = trav3;
					}
					if (trav4 < trav)
					{
						trav = trav4;
					}
					if (trav5 < trav)
					{
						trav = trav5;
					}
					//
				}
			}
			travm = trav;
		}
		return travm;
	}
}
#endif