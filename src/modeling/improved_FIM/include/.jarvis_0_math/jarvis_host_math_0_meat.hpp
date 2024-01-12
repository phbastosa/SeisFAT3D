#ifndef JARVIS_TFUN_MEAT
#define JARVIS_TFUN_MEAT
//**********************************Developer******************************************
// 2020.04.10 BY CAIWEI CALVIN CAI
//*************************************************************************************
#include "jarvis_host_math_0_bones.hpp"
//************************************ARRAY********************************************
namespace jarvis
{
	template <typename eT>
	eT min(eT a, eT b)
	{
		if (a > b)
			return b;
		else
			return a;
	}
	template <typename eT>
	eT min(eT a, eT b, eT c)
	{
		eT tmp;
		tmp = min<eT>(a, b);
		tmp = min<eT>(tmp, c);
		return tmp;
	}
	template <typename eT>
	eT max(eT a, eT b)
	{
		if (a > b)
			return a;
		else
			return b;
	};
	template <typename eT>
	eT max(eT a, eT b, eT c)
	{
		eT tmp;
		tmp = max<eT>(a, b);
		tmp = max<eT>(tmp, c);
		return tmp;
	};
	template <typename eT>
	eT min(eT *a, uint64_t l)
	{
		eT min;
		uint64_t i;
		min = a[0];
		for (i = 0; i < l; i++)
		{
			if (min > a[i])
			{
				min = a[i];
			}
		}
		return min;
	}
	template <typename eT>
	eT min(eT *a, int l)
	{
		eT min;
		uint64_t i;
		min = a[0];
		for (i = 0; i < l; i++)
		{
			if (min > a[i])
			{
				min = a[i];
			}
		}
		return min;
	}
	template <typename eT>
	eT max(eT *a, uint64_t l)
	{
		eT max;
		uint64_t i;
		max = a[0];
		for (i = 0; i < l; i++)
		{
			if (max < a[i])
			{
				max = a[i];
			}
		}
		return max;
	}
	template <typename eT>
	eT max(eT *a, int l)
	{
		eT max;
		int i;
		max = a[0];
		for (i = 0; i < l; i++)
		{
			if (max < a[i])
			{
				max = a[i];
			}
		}
		return max;
	}
	template <typename eT>
	eT min(eT **a, uint32_t n_rows, uint32_t n_cols)
	{
		eT a_min;
		uint32_t i, j;
		a_min = a[0][0];
		for (i = 0; i < n_rows; i++)
		{
			for (j = 0; j < n_cols; j++)
			{
				if (a_min > a[i][j])
				{
					a_min = a[i][j];
				}
			}
		}
		return a_min;
	}
	template <typename eT>
	eT max(eT **a, uint32_t n_rows, uint32_t n_cols)
	{
		eT max;
		uint32_t i, j;
		max = a[0][0];
		for (i = 0; i < n_rows; i++)
		{
			for (j = 0; j < n_cols; j++)
			{
				if (max < a[i][j])
				{
					max = a[i][j];
				}
			}
		}
		return max;
	}
	template <typename eT>
	eT min(eT ***a, uint32_t n_rows, uint32_t n_cols, uint32_t n_slices)
	{
		eT min;
		uint32_t i, j, k;
		min = a[0][0][0];
		for (i = 0; i < n_rows; i++)
			for (j = 0; j < n_cols; j++)
				for (k = 0; k < n_slices; k++)
				{
					if (min > a[i][j][k])
					{
						min = a[i][j][k];
					}
				}
		return min;
	}
	template <typename eT>
	eT max(eT ***a, uint32_t n_rows, uint32_t n_cols, uint32_t n_slices)
	{
		eT _max;
		_max = a[0][0][0];
		uint32_t i, j, k;
		for (i = 0; i < n_rows; i++)
			for (j = 0; j < n_cols; j++)
				for (k = 0; k < n_slices; k++)
				{
					if (_max < a[i][j][k])
					{
						_max = a[i][j][k];
					}
				}
		return _max;
	}
	template <typename eT>
	void min(eT *a, uint64_t l, eT &min_val, uint64_t &pos)
	{
		eT min;
		min = a[0];
		for (uint64_t i = 0; i < l; i++)
		{
			if (min > a[i])
			{
				min = a[i];
				pos = i;
			}
		}
		min_val = min;
	}
	template <typename eT>
	void max(eT *a, uint64_t l, eT &max_val, uint64_t &pos)
	{
		eT max;
		max = a[0];
		uint64_t i;
		for (i = 0; i < l; i++)
		{
			if (max < a[i])
			{
				max = a[i];
				pos = i;
			}
		}
		max_val = max;
	}
	template <typename eT>
	void min(eT **a, uint32_t n_rows, uint32_t n_cols, eT *min_val, uint32_t &pos_x, uint32_t &pos_y)
	{
		eT min;
		uint32_t i, j;
		min = a[0][0];
		for (i = 0; i < n_rows; i++)
			for (j = 0; j < n_cols; j++)
			{
				if (min > a[i][j])
				{
					min = a[i][j];
					pos_x = i;
					pos_y = j;
				}
			}
		min_val = min;
	}
	template <typename eT>
	void max(eT **a, uint32_t n_rows, uint32_t n_cols, eT *max_val, uint32_t &pos_x, uint32_t &pos_y)
	{
		eT max;
		uint32_t i, j;
		max = a[0][0];
		for (i = 0; i < n_rows; i++)
			for (j = 0; j < n_cols; j++)
			{
				if (max < a[i][j])
				{
					max = a[i][j];
					pos_x = i;
					pos_y = j;
				}
			}
		max_val = max;
	}
	template <typename eT>
	eT abs_min(eT *a, uint64_t l)
	{
		uint64_t i;
		eT temp = (eT)fabs((double)a[0]);
		for (i = 1; i < l; i++)
		{
			if (temp > a[i])
				temp = (eT)fabs((double)a[i]);
		}
		return temp;
	}
	template <typename eT>
	eT abs_max(eT *a, uint64_t l)
	{
		uint64_t i;
		eT temp = (eT)fabs((double)a[0]);
		for (i = 1; i < l; i++)
		{
			if (temp < a[i])
				temp = (eT)fabs((double)a[i]);
		}
		return temp;
	}
	template <typename eT>
	eT abs_mean(eT *a, uint64_t l)
	{
		uint64_t i;
		eT Sum = 0;
		for (i = 0; i < l; i++)
		{
			Sum += a[i] * a[i];
		}
		return sqrt(Sum / l);
	}
	template <typename eT>
	eT mean(eT *a, uint64_t l)
	{
		uint64_t i;
		eT _sum = 0;
		for (i = 0; i < l; i++)
		{
			_sum += a[i];
		}
		return _sum / l;
	}
	template <typename eT>
	eT abs_min(eT **a, uint32_t n_rows, uint32_t n_cols)
	{
		uint32_t i, j;
		eT valueout;
		eT temp = (eT)(fabs((double)(a[0][0])));
		for (i = 0; i < n_rows; i++)
		{
			for (j = 0; j < n_cols; j++)
			{
				valueout = (eT)(fabs((double)(a[i][j])));
				if (temp > valueout)
					temp = valueout;
			}
		}
		return temp;
	}
	template <typename eT>
	eT abs_max(eT **a, uint32_t n_rows, uint32_t n_cols)
	{
		uint32_t i, j;
		eT valueout;
		eT temp = (eT)(fabs((double)(a[0][0])));
		for (i = 0; i < n_rows; i++)
		{
			for (j = 0; j < n_cols; j++)
			{
				valueout = (eT)(fabs((double)(a[i][j])));
				if (temp < valueout)
					temp = valueout;
			}
		}
		return temp;
	}
	template <typename eT>
	eT abs_sum(eT *a, int l)
	{
		eT abs_sum = 0;
		int i;
		for (i = 0; i < l; i++)
		{
			abs_sum += abs(a[i]);
		}
		return abs_sum;
	}
	template <typename eT>
	void limit(eT &a, double lb, double rb)
	{
		if (a < lb)
			a = lb;
		else if (a > rb)
			a = rb;
	}
	template <typename eT>
	void UnitVec(eT &a, eT &b)
	{
		eT _a = a;
		eT _b = b;
		a = _a / sqrt(pow(_a, 2) + pow(_b, 2));
		b = _b / sqrt(pow(_a, 2) + pow(_b, 2));
	}
	template <typename eT>
	void UnitVec(eT &a, eT &b, eT &c)
	{
		eT _a = a;
		eT _b = b;
		eT _c = c;
		a = _a / sqrt(pow(_a, 2) + pow(_b, 2) + pow(_c, 2));
		b = _b / sqrt(pow(_a, 2) + pow(_b, 2) + pow(_c, 2));
		c = _c / sqrt(pow(_a, 2) + pow(_b, 2) + pow(_c, 2));
	}
	template <typename T>
	void Sort(T *&a, uint64_t l)
	{
		T t;
		uint64_t i, j, k;
		for (i = 0; i < l - 1; i++)
		{
			k = i;
			for (j = i + 1; j < l; j++)
				if (a[j] < a[k])
					k = j;
			t = a[k];
			a[k] = a[i];
			a[i] = t;
		}
	}
	template <typename T>
	void Diff(T *&a_i, uint64_t l, T di, T *a_o)
	{
		uint64_t i;
		for (i = 0; i < l - 1; i++)
		{
			a_o[i] = (a_i[i + 1] - a_i[i]) / di;
		}
		if (l > 0)
		{
			a_o[l] = a_o[l - 1];
		}
	}

	template <typename T>
	void Swap(T &a, T &b)
	{
		T temp;
		temp = a;
		a = b;
		b = temp;
	}
	// Full Permutation
	template <typename T>
	bool next_permutation(T *p, uint32_t n)
	{
		uint32_t last = n - 1;
		uint32_t i, j, k;
		static uint32_t Count = 0;
		if (Count == 0)
		{
			Count++;
			return true;
		}
		i = last;
		while (i > 0 && p[i] < p[i - 1])
			i--;
		if (i == 0)
			return false;
		k = i;
		for (j = last; j >= i; j--)
			if (p[j] > p[i - 1] && p[j] < p[k])
				k = j;
		Swap<T>(p[k], p[i - 1]);
		for (j = last, k = i; j > k; j--, k++)
			Swap<T>(p[j], p[k]);
		Count++;
		return true;
	}
	// Show Permutation
	template <typename T>
	void show_permutation(T *p, uint32_t n)
	{
		uint32_t i;
		for (i = 0; i < n; i++)
			std::cout << p[i];
		std::cout << std::endl;
	}

	inline void fzero(float *A, uint32_t u_max)
	{
		uint32_t u = 0;
		while (++u <= u_max)
			*(A++) = 0.0;
		// memset(A, 0, sizeof(float)*u_max);
	}
	//*********************************************************************************
}
#endif