#pragma once
#ifndef _TRAVEL_TIME_3D_BONES
#define _TRAVEL_TIME_3D_BONES
#include "include/.jarvis_1_field/.jarvis_2_device_field_header_out.h"
namespace jarvis
{
	struct TimePoint3D : Point3D
	{
		double time = 0;
#ifdef JARVIS_DEBUG
		void print()
		{
			printf("TimePoint3D:(%6.4f, %6.4f, %6.4f): %10.6f\n", x, y, z, time);
		}
#endif
		void set_pos(const Coord3D &a)
		{
			x = a.x;
			y = a.y;
			z = a.z;
		}
		Coord3D get_pos()
		{
			Coord3D a;
			a.x = x;
			a.y = y;
			a.z = z;
			return a;
		}
	};
	enum class TravelType
	{
		theoretical = 1,
		normal,
		refine
	};
	enum class NodeStatus
	{
		null = -1,
		not_active = 0,
		active = 1,
		converging = 2,
		converged = 3
	};
	typedef cuField<TimePoint3D> tp3cuvec;
	class travel_time_3d_module
	{
	public:
		fcufld *vel_p;
		int diff_order;
		int extend_num;
		int divide_num;
		bool is_source_refine = false;
		tp3cuvec shotline;
		cuField<NodeStatus> mark;
		dcufld time;
		cuField<bool> endflag;
		//
		void module_init(fcufld *_vel_p, int _diff_order, int _extend_num = 1, int _divide_num = 1);
		void cal_travel_time(tp3cuvec &_shotline, TravelType _travel_type);
		//*FUNC
		void get_time_at(tp3cuvec &_position);
		void get_device_time();
		dcufld GetTravelTimeField();
		void cu_set_zero();
		void del();

	protected:
		void set_source(tp3cuvec &_shotline);
		void cu_cal_time(tp3cuvec &_shotline);
		void cal_true_time(tp3cuvec &_shotline);
		void cu_cal_refine_time(tp3cuvec &_shotline);
	};
}
#endif