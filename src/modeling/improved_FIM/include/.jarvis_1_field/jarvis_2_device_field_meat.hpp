#ifndef JARVIS_DEVICE_TCUFIELD_MEAT
#define JARVIS_DEVICE_TCUFIELD_MEAT
//**********************************Developer*************************************
// 2021.05.14 BY CAIWEI CALVIN CAI
//********************************************************************************
#include "jarvis_2_device_field_bones.hpp"
//********************************CLASS_TEMPLATE**********************************
namespace jarvis
{
	template <typename eT>
	inline cuField<eT>::cuField()
	{
		this->mem = 0;
		cu_mem = 0;
	}
	//
	template <typename eT>
	inline cuField<eT>::cuField(string _name)
	{
		this->mem = 0;
		this->name = _name;
		cu_mem = 0;
	}
	//
	template <typename eT>
	inline cuField<eT>::cuField(eT *_mem_ptr)
	{
		this->mem = _mem_ptr;
	}
	//
	template <typename eT>
	inline cuField<eT>::cuField(MemType _mem_type, uint32_t _n_rows, uint32_t _n_cols, uint32_t _n_slices)
	{
		cu_alloc(_mem_type, _n_rows, _n_cols, _n_slices);
	}
	//
	template <typename eT>
	inline cuField<eT>::cuField(MemType _mem_type, const Frame &_frame)
	{
		cu_alloc(_mem_type, _frame);
	}
	//
	template <typename eT>
	inline void cuField<eT>::cu_alloc(MemType _mem_type, uint32_t _n_rows, uint32_t _n_cols, uint32_t _n_slices)
	{
		this->frame.set_n(_n_rows, _n_cols, _n_slices);
		mem_type = _mem_type;
		if (mem_type == MemType::pin)
		{
			JARVIS_HANDLE_ERROR(cudaHostAlloc((void **)&this->mem, this->frame.n_elem * sizeof(eT), cudaHostAllocDefault));
		}
		else
			this->mem = (eT *)calloc(this->frame.n_elem, sizeof(eT));
		JARVIS_HANDLE_ERROR(cudaMalloc((void **)&(cu_mem), this->frame.n_elem * sizeof(eT)));
		cudaMemset(cu_mem, 0, this->frame.n_elem * sizeof(eT));
#ifdef JARVIS_DEBUG
		this->mem_state = MemState::allocated;
		cu_mem_state = MemState::allocated;
#endif
	}
	//
	template <typename eT>
	inline void cuField<eT>::cu_alloc(MemType _mem_type, const Frame &_frame)
	{
		cu_alloc_only_host(_mem_type, _frame);
		JARVIS_HANDLE_ERROR(cudaMalloc((void **)&(cu_mem), this->frame.n_elem * sizeof(eT)));
		cudaMemset(cu_mem, 0, this->frame.n_elem * sizeof(eT));
#ifdef JARVIS_DEBUG
		cu_mem_state = MemState::allocated;
#endif
	}
	//
	template <typename eT>
	inline void cuField<eT>::cu_alloc_by_frame(MemType _mem_type)
	{
#ifdef JARVIS_DEBUG
		if (this->mem_state == MemState::inital)
		{
#endif
			mem_type = _mem_type;
			if (mem_type == MemType::pin)
			{
				JARVIS_HANDLE_ERROR(cudaHostAlloc((void **)&this->mem, this->frame.n_elem * sizeof(eT), cudaHostAllocDefault));
			}
			else
			{
				this->mem = (eT *)calloc(this->frame.n_elem, sizeof(eT));
			}
			JARVIS_HANDLE_ERROR(cudaMalloc((void **)&(cu_mem), this->frame.n_elem * sizeof(eT)));
			cudaMemset(cu_mem, 0, this->frame.n_elem * sizeof(eT));
#ifdef JARVIS_DEBUG
			this->mem_state = MemState::allocated;
			cu_mem_state = MemState::allocated;
		}
#ifndef JARVIS_NO_WARNING
		else if (this->mem_state == MemState::allocated)
		{
			printf("alloc():[WARNING]:cuField host memory has been allocated!\n");
		}
#endif
#endif
	}
	//
	template <typename eT>
	inline void cuField<eT>::cu_alloc_only_host(MemType _mem_type, uint32_t _n_rows, uint32_t _n_cols, uint32_t _n_slices)
	{
#ifdef JARVIS_DEBUG
		if (this->mem_state == MemState::inital && cu_mem_state == MemState::inital)
		{
#endif
			this->frame.set_n(_n_rows, _n_cols, _n_slices);
			mem_type = _mem_type;
			if (mem_type == MemType::pin)
			{
				JARVIS_HANDLE_ERROR(cudaHostAlloc((void **)&this->mem, this->frame.n_elem * sizeof(eT), cudaHostAllocDefault));
			}
			else
				this->mem = (eT *)calloc(this->frame.n_elem, sizeof(eT));
#ifdef JARVIS_DEBUG
			this->mem_state = MemState::allocated;
		}
#ifndef JARVIS_NO_WARNING
		else if (this->mem_state == MemState::allocated)
		{
			printf("alloc():[WARNING]:cuField host memory has been allocated!\n");
		}
		else if (cu_mem_state == MemState::allocated)
		{
			printf("alloc():[WARNING]:cuField device memory has been allocated!\n");
		}
#endif
#endif
	}
	//
	template <typename eT>
	inline void cuField<eT>::cu_alloc_only_host(MemType _mem_type, const Frame &_frame)
	{
#ifdef JARVIS_DEBUG
		if (this->mem_state == MemState::inital)
		{
#endif
			this->frame.copy(_frame);
			mem_type = _mem_type;
			if (mem_type == MemType::pin)
			{
				JARVIS_HANDLE_ERROR(cudaHostAlloc((void **)&this->mem, this->frame.n_elem * sizeof(eT), cudaHostAllocDefault));
			}
			else
			{
				this->mem = (eT *)calloc(this->frame.n_elem, sizeof(eT));
			}
#ifdef JARVIS_DEBUG
			this->mem_state = MemState::allocated;
		}
#ifndef JARVIS_NO_WARNING
		else if (this->mem_state == MemState::allocated)
		{
			printf("alloc():[WARNING]:cuField host memory has been allocated!\n");
		}
#endif
#endif
	}
	//
	template <typename eT>
	inline void cuField<eT>::set_frame_for_alloc(uint32_t _n_rows, uint32_t _n_cols, uint32_t _n_slices)
	{
#ifdef JARVIS_DEBUG
		if (this->mem_state == MemState::inital && cu_mem_state == MemState::inital)
		{
#endif
			this->frame.set_n(_n_rows, _n_cols, _n_slices);
			mem_type = MemType::pin;
#ifdef JARVIS_DEBUG
			this->mem_state = MemState::allocated;
			cu_mem_state = MemState::allocated;
		}
#ifndef JARVIS_NO_WARNING
		else if (this->mem_state == MemState::allocated)
		{
			printf("alloc():[WARNING]:cuField host memory has been allocated!\n");
		}
		else if (cu_mem_state == MemState::allocated)
		{
			printf("alloc():[WARNING]:cuField device memory has been allocated!\n");
		}
#endif
#endif
	}
	template <typename eT>
	inline void cuField<eT>::set_frame_for_alloc(const Frame &_frame)
	{
#ifdef JARVIS_DEBUG
		if (this->mem_state == MemState::inital)
		{
#endif
			this->frame.copy(_frame);
			mem_type = MemType::pin;
#ifdef JARVIS_DEBUG
			this->mem_state = MemState::allocated;
			cu_mem_state = MemState::allocated;
		}
#ifndef JARVIS_NO_WARNING
		else if (this->mem_state == MemState::allocated)
		{
			printf("alloc():[WARNING]:cuField host memory has been allocated!\n");
		}
#endif
#endif
	}
	//
	template <typename eT>
	inline void cuField<eT>::set_frame_for_alloc(uint32_t _n_rows, uint32_t _n_cols, uint32_t _n_slices, float _d_rows, float _d_cols, float _d_slices, float _l_rows, float _l_cols, float _l_slices)
	{
#ifdef JARVIS_DEBUG
		if (this->mem_state == MemState::inital)
		{
#endif
			this->frame.set_ndl(_n_rows, _n_cols, _n_slices, _d_rows, _d_cols, _d_slices, _l_rows, _l_cols, _l_slices);
			mem_type = MemType::pin;
#ifdef JARVIS_DEBUG
			this->mem_state = MemState::allocated;
			cu_mem_state = MemState::allocated;
		}
#ifndef JARVIS_NO_WARNING
		else if (this->mem_state == MemState::allocated)
		{
			printf("alloc():[WARNING]:cuField host memory has been allocated!\n");
		}
#endif
#endif
	}

	//
	template <typename eT>
	inline void cuField<eT>::cu_alloc_only_device(uint32_t _n_rows, uint32_t _n_cols, uint32_t _n_slices)
	{
		this->frame.set_n(_n_rows, _n_cols, _n_slices);
#ifdef JARVIS_DEBUG
		cu_mem_state = MemState::allocated;
#endif
		JARVIS_HANDLE_ERROR(cudaMalloc((void **)&(cu_mem), this->frame.n_elem * sizeof(eT)));
		cudaMemset(cu_mem, 0, this->frame.n_elem * sizeof(eT));
	}
	//
	template <typename eT>
	inline void cuField<eT>::cu_alloc_only_device(const Frame &_frame)
	{
		this->frame.copy(_frame);
#ifdef JARVIS_DEBUG
		cu_mem_state = MemState::allocated;
#endif
		JARVIS_HANDLE_ERROR(cudaMalloc((void **)&(cu_mem), this->frame.n_elem * sizeof(eT)));
		cudaMemset(cu_mem, 0, this->frame.n_elem * sizeof(eT));
	}
	//
	template <typename eT>
	inline void cuField<eT>::cu_stream_copy_h2d(cudaStream_t s)
	{
#ifdef JARVIS_DEBUG
		if (this->mem_state == MemState::allocated && cu_mem_state == MemState::allocated && mem_type == MemType::pin)
		{
#endif
			cudaMemcpyAsync(cu_mem, this->mem, this->frame.n_elem * sizeof(eT), cudaMemcpyHostToDevice, s);
#ifdef JARVIS_DEBUG
		}
		else if (mem_type == MemType::npin)
		{
			printf("cu_stream_copy_h2d:\033[41;37m[ERROR]:\033[0mMem is not Pinned!");
			std::abort();
		}
		else if (this->mem_state != MemState::allocated || cu_mem_state != MemState::allocated)
		{
			printf("cu_stream_copy_h2d:\033[41;37m[ERROR]:\033[0mHost and Device Mem not allocated!");
			std::abort();
		}
#endif
	}
	//
	template <typename eT>
	inline void cuField<eT>::cu_stream_copy_h2d(cuField<eT> &a_dst, cudaStream_t s)
	{
#ifdef JARVIS_DEBUG
		if (this->mem_state == MemState::allocated && mem_type == MemType::pin)
		{
#endif
			cudaMemcpyAsync(a_dst.cu_mem, this->mem, this->frame.n_elem * sizeof(eT), cudaMemcpyHostToDevice, s);
#ifdef JARVIS_DEBUG
		}
		else if (mem_type == MemType::npin)
		{
			printf("cu_stream_copy_h2d:\033[41;37m[ERROR]:\033[0mMem is not Pinned!");
			std::abort();
		}
		else if (this->mem_state != MemState::allocated)
		{
			printf("cu_stream_copy_h2d:\033[41;37m[ERROR]:\033[0mHost Mem not allocated!");
			std::abort();
		}
#endif
	}
	//
	template <typename eT>
	inline void cuField<eT>::cu_copy_h2d()
	{
#ifdef JARVIS_DEBUG
		if (this->mem_state == MemState::allocated && cu_mem_state == MemState::allocated && mem_type != MemType::pin)
		{
#endif
			cudaMemcpy(cu_mem, this->mem, this->frame.n_elem * sizeof(eT), cudaMemcpyHostToDevice);
#ifdef JARVIS_DEBUG
		}
		else if (mem_type == MemType::pin)
		{
			printf("cu_copy_h2d:\033[41;37m[ERROR]:\033[0mMem is Pinned! plz use the [stream_copy_h2d]");
			std::abort();
		}
		else if (this->mem_state != MemState::allocated || cu_mem_state != MemState::allocated)
		{
			printf("cu_copy_h2d:\033[41;37m[ERROR]:\033[0mHost and Device Mem not allocated!");
			std::abort();
		}
#endif
	}
	//
	template <typename eT>
	inline void cuField<eT>::cu_copy_h2d(cuField<eT> &a_dst)
	{
		cudaMemcpy(a_dst.cu_mem, this->mem, this->frame.n_elem * sizeof(eT), cudaMemcpyHostToDevice);
	}
	//
	template <typename eT>
	inline void cuField<eT>::cu_stream_copy_d2h(cudaStream_t s)
	{
		cudaMemcpyAsync(this->mem, cu_mem, this->frame.n_elem * sizeof(eT), cudaMemcpyDeviceToHost, s);
	}
	//
	template <typename eT>
	inline void cuField<eT>::cu_stream_copy_d2h(cuField<eT> &a_dst, cudaStream_t s)
	{
		cudaMemcpyAsync(a_dst.mem, cu_mem, this->frame.n_elem * sizeof(eT), cudaMemcpyDeviceToHost, s);
	}
	//
	template <typename eT>
	inline void cuField<eT>::cu_copy_d2h()
	{
		cudaMemcpy(this->mem, cu_mem, this->frame.n_elem * sizeof(eT), cudaMemcpyDeviceToHost);
	}
	//
	template <typename eT>
	inline void cuField<eT>::cu_copy_d2h(cuField<eT> &a_dst)
	{
		cudaMemcpy(a_dst.mem, cu_mem, this->frame.n_elem * sizeof(eT), cudaMemcpyDeviceToHost);
	}
	//
	template <typename eT>
	inline void cuField<eT>::cu_set_zero()
	{
		cudaMemset(cu_mem, 0, this->frame.n_elem * sizeof(eT));
	}
	//
	template <typename eT>
	inline void cuField<eT>::cu_save(string file_path, SaveFormat _format)
	{
		if (this->frame.n_elem == 0)
		{
			printf("\033[41;37m[ERROR]:\033[0m cuField cu_save() it is empty!!\n");
			std::abort();
		}
		if (mem_type == MemType::npin)
		{

			this->cu_copy_d2h();
			cudaDeviceSynchronize();
			this->save(file_path, _format);
		}
		else if (mem_type == MemType::pin)
		{
			this->cu_stream_copy_d2h(jarvis_default_stream);
			cudaDeviceSynchronize();
			this->save(file_path, _format);
		}
		else
		{
			printf("cu_save error!");
			std::abort();
		}
	}
	//
	template <typename eT>
	void cuField<eT>::to_cuField(MemType _mem_type, Field<eT> &obj)
	{
		this->frame.copy(obj.frame);
		cu_alloc_only_device(this->frame);
		this->mem = obj.get_mem_p();
		this->mem_state = MemState::allocated;
		if (_mem_type == MemType::pin)
		{
			cudaHostRegister((void *)&this->mem, this->frame.n_elem * sizeof(eT), cudaHostAllocDefault);
			mem_type = _mem_type;
		}
	}
//
#ifdef JARVIS_DEBUG
#ifdef JARVIS_MEMORY_RELEASE_MANAGEMENT
	template <typename eT>
	cuField<eT>::~cuField()
	{
		if (this->mem_state == MemState::allocated)
		{
			printf("cuField [%s] HOST Memory has not been released!", this->name.c_str());
			std::abort();
		}
		if (cu_mem_state == MemState::allocated)
		{
			printf("cuField [%s] DEVICE Memory has not been released!", this->name.c_str());
			std::abort();
		}
	}
#endif
	template <typename eT>
	inline eT cuField<eT>::cu_get_value(uint32_t i, uint32_t j, uint32_t k)
	{
		if (mem_type == MemType::pin)
			cudaMemcpyAsync(this->mem, cu_mem, this->frame.n_elem * sizeof(eT), cudaMemcpyDeviceToHost, 0);
		else
			cudaMemcpy(this->mem, cu_mem, this->frame.n_elem * sizeof(eT), cudaMemcpyDeviceToHost);
		cudaDeviceSynchronize();
		return this->mem[i + j * this->n_rows + k * this->n_elem_slice];
	}
#endif
	//
	template <typename eT>
	inline MemState cuField<eT>::cu_get_mem_state()
	{
		return cu_mem_state;
	}
	//
	template <typename eT>
	inline void cuField<eT>::host_del()
	{
#ifdef JARVIS_DEBUG
		if (this->mem_state == MemState::allocated)
		{
#endif
			if (mem_type == MemType::npin)
			{
				free(this->mem);
				this->mem = 0;
			}
			else if (mem_type == MemType::pin)
			{
				cudaFreeHost(this->mem);
			}
#ifdef JARVIS_DEBUG
			this->mem_state = MemState::released;
		}
#ifndef JARVIS_NO_WARNING
		else if (this->mem_state == MemState::inital)
		{
			printf("host_del()  :[WARNING]:cuField [%s] HOST memory is not allocated and does not need to be released!\n", this->name.c_str());
		}
		else if (this->mem_state == MemState::released)
		{
			printf("host_del()  :[WARNING]:cuField [%s] HOST does not need to be released again!\n", this->name.c_str());
		}
#endif
#endif
	}
	template <typename eT>
	inline void cuField<eT>::del()
	{
		host_del();
		device_del();
	}
	//
	template <typename eT>
	inline void cuField<eT>::device_del()
	{
#ifdef JARVIS_DEBUG
		if (cu_mem_state == MemState::allocated)
		{
#endif
			cudaFree(cu_mem);
#ifdef JARVIS_DEBUG
			cu_mem_state = MemState::released;
		}
#ifndef JARVIS_NO_WARNING
		else if (cu_mem_state == MemState::inital)
		{
			printf("device_del():[WARNING]:cuField [%s] DEVICE memory is not allocated and does not need to be released!\n", this->name.c_str());
		}
		else if (cu_mem_state == MemState::released)
		{
			printf("device_del():[WARNING]:cuField [%s] DEVICE memory does not need to be released again!\n", this->name.c_str());
		}
#endif
#endif
	}
}
#endif
