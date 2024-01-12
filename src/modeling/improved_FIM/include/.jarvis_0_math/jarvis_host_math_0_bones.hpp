#ifndef JARVIS_TFUN_BONES
#define JARVIS_TFUN_BONES
//**********************************Developer******************************************
// 2020.04.10 BY CAIWEI CALVIN CAI
//*************************************************************************************
#include ".jarvis_host_math_header_in.h"
//************************************ARRAY********************************************
namespace jarvis
{
	template <typename eT>
	eT min(eT a, eT b);
	template <typename eT>
	eT max(eT a, eT b);
	template <typename eT>
	eT min(eT a, eT b, eT c);
	template <typename eT>
	eT max(eT a, eT b, eT c);
	template <typename eT>
	eT min(eT *a, int l);
	template <typename eT>
	eT max(eT *a, int l);
	template <typename eT>
	eT min(eT **a, uint32_t n_rows, uint32_t n_cols);
	template <typename eT>
	eT max(eT **a, uint32_t n_rows, uint32_t n_cols);
	template <typename eT>
	eT min(eT ***a, uint32_t n_rows, uint32_t n_cols, uint32_t n_slices);
	template <typename eT>
	eT max(eT ***a, uint32_t n_rows, uint32_t n_cols, uint32_t n_slices);
	template <typename eT>
	eT abs_min(eT *a, uint64_t l);
	template <typename eT>
	eT abs_max(eT *a, uint64_t l);
	template <typename eT>
	eT abs_min(eT **a, uint32_t n_rows, uint32_t n_cols);
	template <typename eT>
	eT abs_max(eT **a, uint32_t n_rows, uint32_t n_cols);
	template <typename eT>
	eT mean(eT *a, uint64_t l);
	template <typename eT>
	eT abs_mean(eT *a, uint64_t l);
	template <typename eT>
	void min(eT *a, uint64_t l, eT &min_val, uint64_t &pos);
	template <typename eT>
	void max(eT *a, uint64_t l, eT &max_val, uint64_t &pos);
	template <typename eT>
	void min(eT **a, int n_rows, uint32_t n_cols, eT *min_val, uint32_t &pos_x, uint32_t &pos_y);
	template <typename eT>
	void max(eT **a, uint32_t n_rows, uint32_t n_cols, eT *max_val, uint32_t &pos_x, uint32_t &pos_y);
	template <typename eT>
	void limit(eT &a, double lb, double rb);
	// Full Permutation
	template <typename T>
	bool next_permutation(T *p, uint32_t n);
	// Show Permutation
	template <typename T>
	void show_permutation(T *p, uint32_t n);
	//*********************************************************************************
	template <typename T>
	class Component2
	{
	public:
		T x;
		T y;
		enum class c
		{
			x = 0,
			y,
		};
	};
	struct Coord2D : public Component2<double>
	{
		void print();
		double GetDistanceTo(Coord2D a);
	};
	//
	struct Point2D : Coord2D
	{
		uint32_t ind;
		void print();
	};
	//
	struct uvec2 : public Component2<double>
	{
		bool is_unit()
		{
			if (x * x + y * y > 0.999999)
				return true;
			else
				return false;
		}
	};
	//
	template <typename T>
	class Component3
	{
	public:
		T x;
		T y;
		T z;
		enum class c
		{
			x = 0,
			y,
			z
		};
	};
	//
	typedef Component3<float> fcomp3;
	struct Coord3D : public Component3<float>
	{
		void print();
		double GetDistanceTo(Coord3D a);
	};
	//
	struct Point3D : Coord3D
	{
		uint32_t ind;
		void print();
		void set_pos(Coord3D a);
		Coord3D get_pos();
	};
	//
	struct uvec3 : public Component3<float>
	{
		bool is_unit()
		{
			if (x * x + y * y + z * z > 0.999)
				return true;
			else
				return false;
		}
	};
}
#endif