#ifndef JARVIS_DEVICE_TCUFIELD_BONES
#define JARVIS_DEVICE_TCUFIELD_BONES
//**********************************Developer*************************************
// 2021.05.14 BY CAIWEI CALVIN CAI
//********************************************************************************
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "jarvis_1_host_field_meat.hpp"
//********************************CLASS_TEMPLATE**********************************
namespace jarvis
{
	enum class MemType
	{
		npin = 0,
		pin
	};

	template <typename eT>
	class cuField : public Field<eT>
	{
	protected:
		MemState cu_mem_state = MemState::inital;
		MemType mem_type = MemType::npin;

	public:
		eT *cu_mem = 0;

	public:
		cuField();
		cuField(string _name);
		cuField(eT *_mem_ptr);
		cuField(MemType _mem_type, uint32_t _n_rows, uint32_t _n_cols = 1, uint32_t _n_slices = 1);
		cuField(MemType _mem_type, const Frame &_frame);
		void cu_alloc(MemType _mem_type, uint32_t _n_rows, uint32_t _n_cols = 1, uint32_t _n_slices = 1);
		void cu_alloc(MemType _mem_type, const Frame &_frame);
		void cu_alloc_only_host(MemType _mem_type, uint32_t _n_rows, uint32_t _n_cols = 1, uint32_t _n_slices = 1);
		void cu_alloc_only_host(MemType _mem_type, const Frame &_frame);
		void cu_alloc_only_device(uint32_t _n_rows, uint32_t _n_cols = 1, uint32_t _n_slices = 1);
		void cu_alloc_only_device(const Frame &_frame);
		void cu_copy_h2d();
		void cu_copy_d2h();
		void cu_copy_h2d(cuField<eT> &_dst);
		void cu_copy_d2h(cuField<eT> &_dst);
		void cu_stream_copy_h2d(cudaStream_t s);
		void cu_stream_copy_d2h(cudaStream_t s);
		void cu_stream_copy_h2d(cuField<eT> &_dst, cudaStream_t s);
		void cu_stream_copy_d2h(cuField<eT> &_dst, cudaStream_t s);
		void cu_save(string file_path, SaveFormat _format);
		void cu_set_zero();
		void set_frame_for_alloc(uint32_t _n_rows, uint32_t _n_cols = 1, uint32_t _n_slices = 1);
		void set_frame_for_alloc(const Frame &_frame);
		void set_frame_for_alloc(uint32_t _n_rows, uint32_t _n_cols, uint32_t _n_slices, float _d_rows, float _d_cols, float _d_slices, float _l_rows, float _l_cols, float _l_slices);
		void cu_alloc_by_frame(MemType _mem_type);
		MemState cu_get_mem_state();
		void to_cuField(MemType _mem_type, Field<eT> &obj);

		void host_del();
		void device_del();
		void del();
		eT cu_get_value(uint32_t i, uint32_t j = 0, uint32_t k = 0);
#ifdef JARVIS_DEBUG
	public:
#ifdef JARVIS_MEMORY_RELEASE_MANAGEMENT
		~cuField();
#endif
#endif
	};
	//******************************************************************
	typedef cuField<double> cufld;
	typedef cuField<double> dcufld;
	typedef cuField<float> fcufld;
	typedef cuField<int> icufld;
	//******************************************************************
#define set_cufield_2d_idx(mesh, idx, i, j)               \
	int n_rows = (mesh).n_rows, n_cols = (mesh).n_cols;   \
	float d_rows = (mesh).d_rows, d_cols = (mesh).d_cols; \
	int idx = threadIdx.x + blockDim.x * blockIdx.x;      \
	int i = idx % n_rows;                                 \
	int j = idx / n_rows

#define set_cufield_3d_idx(grid, idx, i, j, k)                                        \
	int n_rows = (grid).n_rows, n_cols = (grid).n_cols, n_slices = (grid).n_slices;   \
	float d_rows = (grid).d_rows, d_cols = (grid).d_cols, d_slices = (grid).d_slices; \
	int n_elem_slice = (grid).n_elem_slice;                                           \
	int idx = threadIdx.x + blockDim.x * blockIdx.x;                                  \
	int idx_n_elem_slice = (n_rows) * (n_cols);                                       \
	int idx_n_elem = idx_n_elem_slice * (n_slices);                                   \
	int i = (idx % (n_elem_slice)) % n_rows;                                          \
	int j = (idx % (n_elem_slice)) / n_rows;                                          \
	int k = idx / (n_elem_slice)

#define set_cuda_cube_idx(n_rows, n_cols, n_slices, idx, i, j, k) \
	int idx = threadIdx.x + blockDim.x * blockIdx.x;              \
	int idx_n_elem_slice = (n_rows) * (n_cols);                   \
	int idx_n_elem = idx_n_elem_slice * (n_slices);               \
	int i = idx % (idx_n_elem_slice) % (n_rows);                  \
	int j = idx % (idx_n_elem_slice) / (n_rows);                  \
	int k = idx / (idx_n_elem_slice)

#ifndef jarvis_default_stream
#define jarvis_default_stream (cudaStream_t)0
#endif

	static void HandleError(cudaError_t err, const char *file, int line)
	{
		if (err != cudaSuccess)
		{
			printf("%s in %s at line %d\n", cudaGetErrorString(err), file, line);
			exit(EXIT_FAILURE);
		}
	}
#define JARVIS_HANDLE_ERROR(err) (HandleError(err, __FILE__, __LINE__))

#define CUDA_TIC(a)                       \
	cudaEvent_t startTime##a, endTime##a; \
	cudaEventCreate(&startTime##a);       \
	cudaEventCreate(&endTime##a);         \
	cudaEventRecord(startTime##a, 0)

#define CUDA_TOC(a, str)                                      \
	cudaEventRecord(endTime##a, 0);                           \
	cudaEventSynchronize(startTime##a);                       \
	cudaEventSynchronize(endTime##a);                         \
	float time##a;                                            \
	cudaEventElapsedTime(&time##a, startTime##a, endTime##a); \
	printf(" %s : %f ms \n", str, time##a)

#define jarvis_cuda_kernel_size(cuda_threads_num) \
	(cuda_threads_num + jarvis_const::block_size - 1) / jarvis_const::block_size, jarvis_const::block_size

#ifndef GPU_NUM_PER_NODE
#define GPU_NUM_PER_NODE 4
#endif
	//******************************************************************
	struct jarvisCudaStream
	{
		cudaStream_t cal_stream;
		cudaStream_t d2h_stream;
		cudaStream_t h2d_stream;
		jarvisCudaStream()
		{
			cudaStreamCreate(&cal_stream);
			cudaStreamCreate(&d2h_stream);
			cudaStreamCreate(&h2d_stream);
		}
		void sync_cal_stream()
		{
			cudaStreamSynchronize(cal_stream);
		}
		void sync_d2h_stream()
		{
			cudaStreamSynchronize(d2h_stream);
		}
		void sync_h2d_stream()
		{
			cudaStreamSynchronize(h2d_stream);
		}
	};

	static void jarvisGetGPUMemory(int i_device)
	{
		int deviceCount = 0;
		cudaError_t error_id = cudaGetDeviceCount(&deviceCount);
		if (deviceCount == 0)
		{
			std::cout << "\033[41;37mThe current PC does not support CUDA devices!\033[0m" << std::endl;
		}

		size_t gpu_total_size;
		size_t gpu_free_size;

		cudaError_t cuda_status = cudaMemGetInfo(&gpu_free_size, &gpu_total_size);

		if (cudaSuccess != cuda_status)
		{
			std::cout << "Error: cudaMemGetInfo fails : " << cudaGetErrorString(cuda_status) << std::endl;
			exit(1);
		}

		double total_memory = double(gpu_total_size) / (1024.0 * 1024.0);
		double free_memory = double(gpu_free_size) / (1024.0 * 1024.0);
		double used_memory = total_memory - free_memory;

		string free_memory_str = "[free\tcuda memory]:\t";
		if (free_memory < 10 * 1024.0)
		{
			std::cout << "[device\tnumber]:\t" << i_device << "\n"
					  << "[total\tcuda memory]:\t" << total_memory << " MB\n"
					  << "[used\tcuda memory]:\t" << used_memory << " MB\n"
					  << "\033[41;37m[free\tcuda memory]:\t" << free_memory << " MB\033[0m\n"
					  << std::endl;
		}
		else
		{
			std::cout << "[device\tnumber]:\t" << i_device << "\n"
					  << "[total\tcuda memory]:\t" << total_memory << " MB\n"
					  << "[used\tcuda memory]:\t" << used_memory << " MB\n"
					  << "[free\tcuda memory]:\t" << free_memory << " MB\n"
					  << std::endl;
		}

		if (free_memory < 10 * 1024.0)
		{
			printf("\033[41;37m[WARNING]: free memory is not enough!\033[0m\n");
			std::abort();
		}
	}
	//
	inline void CheckCudaInfo()
	{
		cudaDeviceProp deviceProp;
		int deviceCount;
		cudaGetDeviceCount(&deviceCount);
		for (int i = 0; i < deviceCount; i++)
		{
			cudaGetDeviceProperties(&deviceProp, i);
			cout << "Device " << i + 1 << ":" << endl;
			cout << "Graphics card model:" << deviceProp.name << endl;
			cout << "Total device global memory in MB:"
				 << deviceProp.totalGlobalMem / 1024 / 1024 << endl;
			cout << "The maximum available shared memory (in KB) in a thread block on the device:"
				 << deviceProp.sharedMemPerBlock / 1024 << endl;
			cout << "Number of 32-bit registers available for a thread block on the device:"
				 << deviceProp.regsPerBlock << endl;
			cout << "The maximum number of threads that a thread block on a device can contain:"
				 << deviceProp.maxThreadsPerBlock << endl;
			cout << "Version number of the device's compute capability:"
				 << deviceProp.major << "." << deviceProp.minor << endl;
			cout << "Number of multiprocessors on the device:" << deviceProp.multiProcessorCount << endl;
			cout << "canMapHostMemory:" << deviceProp.canMapHostMemory << endl;
			int canAccessPeer;

			cudaDeviceCanAccessPeer(&canAccessPeer, 0, 4);
			cout << "support peer to peer==%d" << canAccessPeer << endl;
		}
	}

	inline __device__ bool cuJacobi(float *matrix, int dim, float *eigenvectors, float *eigenvalues, float precision, int max)
	{
		for (int i = 0; i < dim; i++)
		{
			eigenvectors[i * dim + i] = 1.0f;
			for (int j = 0; j < dim; j++)
			{
				if (i != j)
					eigenvectors[i * dim + j] = 0.0f;
			}
		}

		int nCount = 0; // current iteration
		while (1)
		{
			// find the largest element on the off-diagonal line of the matrix
			double dbMax = matrix[1];
			int nRow = 0;
			int nCol = 1;
			for (int i = 0; i < dim; i++)
			{ // row
				for (int j = 0; j < dim; j++)
				{ // column
					double d = fabs(matrix[i * dim + j]);
					if ((i != j) && (d > dbMax))
					{
						dbMax = d;
						nRow = i;
						nCol = j;
					}
				}
			}

			if (dbMax < precision) // precision check
				break;
			if (nCount > max) // iterations check
				break;
			nCount++;

			double dbApp = matrix[nRow * dim + nRow];
			double dbApq = matrix[nRow * dim + nCol];
			double dbAqq = matrix[nCol * dim + nCol];
			// compute rotate angle
			double dbAngle = 0.5 * atan2(-2 * dbApq, dbAqq - dbApp);
			double dbSinTheta = sin(dbAngle);
			double dbCosTheta = cos(dbAngle);
			double dbSin2Theta = sin(2 * dbAngle);
			double dbCos2Theta = cos(2 * dbAngle);
			matrix[nRow * dim + nRow] = dbApp * dbCosTheta * dbCosTheta +
										dbAqq * dbSinTheta * dbSinTheta + 2 * dbApq * dbCosTheta * dbSinTheta;
			matrix[nCol * dim + nCol] = dbApp * dbSinTheta * dbSinTheta +
										dbAqq * dbCosTheta * dbCosTheta - 2 * dbApq * dbCosTheta * dbSinTheta;
			matrix[nRow * dim + nCol] = 0.5 * (dbAqq - dbApp) * dbSin2Theta + dbApq * dbCos2Theta;
			matrix[nCol * dim + nRow] = matrix[nRow * dim + nCol];

			for (int i = 0; i < dim; i++)
			{
				if ((i != nCol) && (i != nRow))
				{
					int u = i * dim + nRow; // p
					int w = i * dim + nCol; // q
					dbMax = matrix[u];
					matrix[u] = matrix[w] * dbSinTheta + dbMax * dbCosTheta;
					matrix[w] = matrix[w] * dbCosTheta - dbMax * dbSinTheta;
				}
			}

			for (int j = 0; j < dim; j++)
			{
				if ((j != nCol) && (j != nRow))
				{
					int u = nRow * dim + j; // p
					int w = nCol * dim + j; // q
					dbMax = matrix[u];
					matrix[u] = matrix[w] * dbSinTheta + dbMax * dbCosTheta;
					matrix[w] = matrix[w] * dbCosTheta - dbMax * dbSinTheta;
				}
			}

			// compute eigenvector
			for (int i = 0; i < dim; i++)
			{
				int u = i * dim + nRow; // p
				int w = i * dim + nCol; // q
				dbMax = eigenvectors[u];
				eigenvectors[u] = eigenvectors[w] * dbSinTheta + dbMax * dbCosTheta;
				eigenvectors[w] = eigenvectors[w] * dbCosTheta - dbMax * dbSinTheta;
			}
		}

		for (int i = 0; i < dim; i++)
		{
			eigenvalues[i] = matrix[i * dim + i];
		}
		return true;
	}
}
#endif