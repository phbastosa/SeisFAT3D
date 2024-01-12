#pragma once
#ifndef TRAVEL_TIME_3D_MEAT
#define TRAVEL_TIME_3D_MEAT
#include "travel_time_3d_2_kernel_meat.hpp"
namespace jarvis
{
		inline void check_shot_line_is_out_of_bounds(Frame &model, tp3cuvec &shotline)
	{
		for (int i = 0; i < shotline.frame.n_elem; i++)
		{
			if (shotline(i).x < model.l_rows)
			{
				printf("ERROR:The shotline exceeds the model boundary!\n");
				printf("DETAILS:The %d shotpoint exceeds the model LEFT boundary!\n", i);
				std::abort();
			}
			if (shotline(i).x > model.r_rows)
			{
				printf("ERROR:The shotline exceeds the model boundary!\n");
				printf("DETAILS:The %d shotpoint exceeds the model RIGHT boundary!\n", i);
				std::abort();
			}
			if (shotline(i).y < model.l_cols)
			{
				printf("ERROR:The shotline exceeds the model boundary!\n");
				printf("DETAILS:The %d shotpoint exceeds the model FRONT boundary!\n", i);
				std::abort();
			}
			if (shotline(i).y > model.r_cols)
			{
				printf("ERROR:The shotline exceeds the model boundary!\n");
				printf("DETAILS:The %d shotpoint exceeds the model BEHIND boundary!\n", i);
				std::abort();
			}
			if (shotline(i).z < model.l_slices)
			{
				printf("ERROR:The shotline exceeds the model boundary!\n");
				printf("DETAILS:The %d shotpoint exceeds the model UP boundary!\n", i);
				std::abort();
			}
			if (shotline(i).z > model.r_slices)
			{
				printf("ERROR:The rshotline exceeds the model boundary!\n");
				printf("DETAILS:The %d shotpoint exceeds the model BOTTOM boundary!\n", i);
				std::abort();
			}
		}
	}

	inline void _get_time_at(dfield &time, tp3cuvec &position) //*Solve time by interpolation
	{
		set_frame_ndl(time.frame);
		int isx, isy, isz;
		float sx, sy, sz;
		float vecl;
		int idx;
		for (int ishot = 0; ishot < position.frame.n_elem; ishot++)
		{
			isx = int((position(ishot).x - l_rows) / d_rows);
			isy = int((position(ishot).y - l_cols) / d_cols);
			isz = int((position(ishot).z - l_slices) / d_slices);
			if (isx == n_rows)
				isx = n_rows - 1;
			if (isy == n_cols)
				isy = n_cols - 1;
			if (isz == n_slices)
				isz = n_slices - 1;
			if (isx >= n_rows - 1)
				isx = n_rows - 2;

			idx = isx + isy * time.frame.n_rows + isz * time.frame.n_elem_slice;
			float vec[3];
			if (isx <= time.frame.n_rows - 2 && isx >= 0)
			{
				vec[0] = time[idx + 1] - time[idx];
			}
			else
			{
				vec[0] = time[idx] - time[idx - 1];
			}
			//
			if (isy <= time.frame.n_cols - 2 && isy >= 0)
			{
				vec[1] = time[idx + time.frame.n_rows] - time[idx];
			}
			else
			{
				vec[1] = time[idx] - time[idx - time.frame.n_rows];
			}
			//
			if (isz <= time.frame.n_slices - 2 && isz >= 0)
			{
				vec[2] = time[idx + time.frame.n_elem_slice] - time[idx];
			}
			else
			{
				vec[2] = time[idx] - time[idx - time.frame.n_elem_slice];
			}
			sx = position(ishot).x - d_rows * isx;
			sy = position(ishot).y - d_cols * isy;
			sz = position(ishot).z - d_slices * isz;
			vecl = sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
			float tdis = d_rows / vecl * vec[0];
			float ptdis = (sx * vec[0] + sy * vec[1] + sz * vec[2]) / vecl;
			position(ishot).time = time[idx] + ptdis / tdis * vec[0];
		}
	}
	
	inline void travel_time_3d_module::module_init(fcufld *_vel_p, int _diff_order, int _extend_num, int _divide_num)
	{
#ifdef JARVIS_DEBUG
		if (_diff_order != 11 && _diff_order != 12 && _diff_order != 21 && _diff_order != 22)
		{
			cout << "travel_time_3d_module::Init():diff_order= 11 or 12 or 21 or 22!" << endl;
			std::abort();
		}
#endif
		vel_p = _vel_p;
		diff_order = _diff_order;
		extend_num = _extend_num;
		divide_num = _divide_num;
		vel_p->set_name("vel");
		time.set_name("time");
		mark.set_name("mark");
		time.cu_alloc(MemType::pin, vel_p->frame);
		mark.cu_alloc(MemType::pin, vel_p->frame);
		endflag.cu_alloc(MemType::pin, 1);
	}

	inline void travel_time_3d_module::cal_travel_time(tp3cuvec &_shotline, TravelType _travel_type)
	{
		if (_travel_type == TravelType::theoretical)
		{
			cal_true_time(_shotline);
		}
		if (extend_num > 1 && divide_num > 1)
		{
			is_source_refine = true;
		}
		//
		if (_travel_type == TravelType::normal)
		{
			if (is_source_refine == false)
			{
				cu_cal_time(_shotline);
			}
			else
			{
				printf("[ERROR]:cal_travel_time:is_source_refine is TRUE!");
				std::abort();
			}
		}
		//
		if (_travel_type == TravelType::refine)
		{
			if (is_source_refine == true)
			{
				cu_cal_refine_time(_shotline);
			}
			else
			{
				printf("[ERROR]:cal_travel_time:is_source_refine is FALSE!");
				std::abort();
			}
		}
	}

	inline dcufld travel_time_3d_module::GetTravelTimeField()
	{
#ifdef JARVIS_DEBUG
		if (&time[0])
		{
#endif
			return time;
#ifdef JARVIS_DEBUG
		}
		else if (&time[0] && !time[0])
		{
			std::cout << "travel_time_3d_module::GetTravelTimeField():\033[41;37m[ERROR]:\033[0mThe time of memory has been released by the copied object !" << std::endl;
			std::abort();
		}
		else
		{
			std::cout << "travel_time_3d_module::GetTravelTimeField():\033[41;37m[ERROR]:\033[0mtravel_time_3d_module Error!" << std::endl;
			std::abort();
		}
#endif
	}

	inline void travel_time_3d_module::get_device_time()
	{
		time.cu_stream_copy_d2h(jarvis_default_stream);
		cudaDeviceSynchronize();
	}

	inline void travel_time_3d_module::get_time_at(tp3cuvec &position)
	{
		_get_time_at(time, position);
	}

	inline void travel_time_3d_module::cal_true_time(tp3cuvec &_shotline)
	{
		set_frame_ndl(vel_p->frame);
		float sx, sy, sz, ds;
		for (int k = 0; k < n_slices; k++)
			for (int j = 0; j < n_cols; j++)
				for (int i = 0; i < n_rows; i++)
				{
					sx = l_rows + i * d_rows - _shotline(0).x;
					sy = l_cols + j * d_cols - _shotline(0).y;
					sz = l_slices + k * d_slices - _shotline(0).z;
					ds = sqrt(sx * sx + sy * sy + sz * sz);
					time(i, j, k) = ds / vel_p->mem_p()[0];
				}
	}

	inline void GetSubGridValue(const fcufld &a, fcufld &sub_a)
	{
		set_frame_ndl(a.frame);
		long sub_mnx = sub_a.frame.n_rows, sub_mnz = sub_a.frame.n_cols;
		double sub_mdx = sub_a.frame.d_rows, sub_mdz = sub_a.frame.d_cols;
		double sub_mlx = sub_a.frame.l_rows, sub_mlz = sub_a.frame.l_cols;
		int isx, isy, isz;
		double sx, sy, sz;
		for (int k = 0; k < sub_a.frame.n_slices; k++)
			for (int j = 0; j < sub_a.frame.n_cols; j++)
				for (int i = 0; i < sub_a.frame.n_rows; i++)
				{
					double posx = sub_a.frame.l_rows + sub_a.frame.d_rows * i;
					double posy = sub_a.frame.l_cols + sub_a.frame.d_cols * j;
					double posz = sub_a.frame.l_slices + sub_a.frame.d_slices * k;

					isx = int((posx - l_rows) / d_rows);
					isy = int((posy - l_cols) / d_cols);
					isz = int((posz - l_slices) / d_slices);
					if (isx == n_rows)
						isx = n_rows - 1;
					if (isx >= n_rows - 1)
						isx = n_rows - 2;
					//
					if (isy == n_cols)
						isy = n_cols - 1;
					if (isy >= n_cols - 1)
						isy = n_cols - 2;
					//
					if (isz == n_slices)
						isz = n_slices - 1;
					if (isz >= n_slices - 1)
						isz = n_slices - 2;

					field vec(3);
					vec(0) = double(a(isx + 1, isy, isz)) - a(isx, isy, isz);
					vec(1) = double(a(isx, isy + 1, isz)) - a(isx, isy, isz);
					vec(2) = double(a(isx, isy, isz + 1)) - a(isx, isy, isz);
					sx = posx - d_rows * isx;
					sy = posy - d_cols * isy;
					sz = posz - d_slices * isz;
					double dis = (sx * vec(0) + sy * vec(1) + sz * vec(2));
					sub_a(i, j, k) = float(dis / d_rows + a(isx, isz));
				}
	}

	inline void travel_time_3d_module::cu_cal_refine_time(tp3cuvec &_shotline)
	{
		TIC(0);
#ifdef JARVIS_DEBUG
		if (_shotline.frame.n_elem == 1 && extend_num > 0 && divide_num > 0)
		{
#endif
			double lbx = _shotline(0).x - extend_num * vel_p->frame.d_rows;
			double rbx = _shotline(0).x + extend_num * vel_p->frame.d_rows;
			double lby = _shotline(0).y - extend_num * vel_p->frame.d_cols;
			double rby = _shotline(0).y + extend_num * vel_p->frame.d_cols;
			double lbz = _shotline(0).z - extend_num * vel_p->frame.d_slices;
			double rbz = _shotline(0).z + extend_num * vel_p->frame.d_slices;

			if (lbx < vel_p->frame.l_rows)
				lbx = vel_p->frame.l_rows;
			if (rbx > vel_p->frame.r_rows)
				rbx = vel_p->frame.r_rows;
			//
			if (lby < vel_p->frame.l_cols)
				lby = vel_p->frame.l_cols;
			if (rby > vel_p->frame.r_cols)
				rby = vel_p->frame.r_cols;
			//
			if (lbz < vel_p->frame.l_slices)
				lbz = vel_p->frame.l_slices;
			if (rbz > vel_p->frame.r_slices)
				rbz = vel_p->frame.r_slices;
			//
			Frame refine_grid;
			refine_grid.l_rows = lbx;
			refine_grid.r_rows = rbx;
			refine_grid.l_cols = lby;
			refine_grid.r_cols = rby;
			refine_grid.l_slices = lbz;
			refine_grid.r_slices = rbz;
			//
			refine_grid.d_rows = vel_p->frame.d_rows / divide_num;
			refine_grid.d_cols = vel_p->frame.d_cols / divide_num;
			refine_grid.d_slices = vel_p->frame.d_slices / divide_num;
			refine_grid.n_rows = long((refine_grid.r_rows - refine_grid.l_rows) / refine_grid.d_rows + 1);
			refine_grid.n_cols = long((refine_grid.r_cols - refine_grid.l_cols) / refine_grid.d_cols + 1);
			refine_grid.n_slices = long((refine_grid.r_slices - refine_grid.l_slices) / refine_grid.d_slices + 1);
			//
			refine_grid.print_info("refine_grid:");
			cout << diff_order << endl;
			fcufld refinevel;
			refinevel.cu_alloc(MemType::pin, refine_grid);
			GetSubGridValue(*vel_p, refinevel);
			// refinevel.save("/home/calvin/WKSP/INCLUDE/jarvis/app/cpp/travel_time_3d/demo/out/refinevel" + refinevel.size_str(), SaveFormat::binary_raw);
			//
			travel_time_3d_module travelrefine;
			travelrefine.module_init(&refinevel, diff_order);
			travelrefine.cu_cal_time(_shotline);
			travelrefine.get_device_time();
			travelrefine.time.save("/home/calvin/WKSP/INCLUDE/jarvis/app/cpp/travel_time_3d/demo/out/refinetime" + travelrefine.time.size_str(), SaveFormat::binary_raw);
			//
			tp3cuvec refineline(MemType::pin, ((refine_grid.n_rows - 1) / divide_num + 1) *
												  ((refine_grid.n_cols - 1) / divide_num + 1) *
												  ((refine_grid.n_slices - 1) / divide_num + 1));
			//
			int linenum = 0;
			for (int k = 0; k < refine_grid.n_slices; k = k + divide_num)
				for (int j = 0; j < refine_grid.n_cols; j = j + divide_num)
					for (int i = 0; i < refine_grid.n_rows; i = i + divide_num)
					{
						refineline(linenum).x = refine_grid.l_rows + refine_grid.d_rows * i;
						refineline(linenum).y = refine_grid.l_cols + refine_grid.d_cols * j;
						refineline(linenum).z = refine_grid.l_slices + refine_grid.d_slices * k;
						refineline(linenum).time = travelrefine.time(i, j, k);
						linenum++;
					}
			cu_cal_time(refineline);
			//
#ifdef JARVIS_DEBUG
		}
		else if (_shotline.frame.n_elem > 1)
		{
			cout << "travel_time_3d_module::RefineTravel():\033[41;37m[ERROR]:\033[0mThis function is only applicable to point sources!";
			std::abort();
		}
		else if (extend_num <= 0)
		{
			cout << "travel_time_3d_module::RefineTravel:\033[41;37m[ERROR]:\033[0mError extend_num!" << endl;
			std::abort();
		}
		else if (divide_num <= 0)
		{
			cout << "travel_time_3d_module::RefineTravel:\033[41;37m[ERROR]:\033[0mError divide_num!" << endl;
			std::abort();
		}
		else
		{
			cout << "travel_time_3d_module::RefineTravel:\033[41;37m[ERROR]:\033[0mError !" << endl;
			std::abort();
		}
#endif
	}

	inline void travel_time_3d_module::cu_set_zero()
	{
		mark.cu_set_zero();
		time.cu_set_zero();
	}

	inline void travel_time_3d_module::del()
	{
		mark.del();
		time.del();
	}

	inline void travel_time_3d_module::set_source(tp3cuvec &_shotline)
	{
		set_frame_ndl(vel_p->frame);
		check_shot_line_is_out_of_bounds(vel_p->frame, _shotline);
		shotline = _shotline;
		mark.fill(NodeStatus::null);
		time.fill(LONG_MAX);
		float sx, sy, sz, ds;
		float ishotx, ishoty, ishotz;
		int isx, isy, isz;
		for (int ishot = 0; ishot < shotline.frame.n_elem; ishot++)
		{
			isx = int((shotline(ishot).x - l_rows) / d_rows);
			isy = int((shotline(ishot).y - l_cols) / d_cols);
			isz = int((shotline(ishot).z - l_slices) / d_slices);
			// ishotx = int((100 - l_rows) / d_rows);
			// ishoty = int((100 - l_cols) / d_cols);
			// ishotz = int((100 - l_slices) / d_slices);
			// if (abs(isx - ishotx) + abs(isy - ishoty) + abs(isz - ishotz) < 10)
			{
				mark(isx, isy, isz) = NodeStatus::converged; //加密网格标记
				time(isx, isy, isz) = shotline(ishot).time;
				//
				int lax = -1;
				int rax = 1;
				int lay = -1;
				int ray = 1;
				int laz = -1;
				int raz = 1;

				if (isx <= 0)
					lax = 1;
				if (isx >= n_rows - 1)
					rax = -1;

				if (isy <= 0)
					lay = 1;
				if (isy >= n_cols - 1)
					ray = -1;

				if (isz <= 0)
					laz = 1;
				if (isz >= n_slices - 1)
					raz = -1;

				if (mark(isx + lax, isy, isz) != NodeStatus::converged)
					mark(isx + lax, isy, isz) = NodeStatus::active;

				if (mark(isx + rax, isy, isz) != NodeStatus::converged)
					mark(isx + rax, isy, isz) = NodeStatus::active;

				if (mark(isx, isy + lay, isz) != NodeStatus::converged)
					mark(isx, isy + lay, isz) = NodeStatus::active;

				if (mark(isx, isy + ray, isz) != NodeStatus::converged)
					mark(isx, isy + ray, isz) = NodeStatus::active;

				if (mark(isx, isy, isz + laz) != NodeStatus::converged)
					mark(isx, isy, isz + laz) = NodeStatus::active;

				if (mark(isx, isy, isz + raz) != NodeStatus::converged)
					mark(isx, isy, isz + raz) = NodeStatus::active;
			}
		}
		time.cu_stream_copy_h2d(jarvis_default_stream);
		mark.cu_stream_copy_h2d(jarvis_default_stream);
		cudaDeviceSynchronize();
	}
}
#endif