#ifndef JARVIS_HOST_TFIELD_BONES
#define JARVIS_HOST_TFIELD_BONES
//**********************************Developer*************************************
// 2021.04.9 BY CAIWEI CALVIN CAI
//********************************************************************************
#include "jarvis_0_host_frame_meat.hpp"
//********************************CLASS_TEMPLATE**********************************
namespace jarvis
{
	enum class MemState
	{
		released = -1,
		inital = 0,
		allocated = 1
	};

	enum class SaveFormat
	{
		ascii_txt,
		ascii_grd,
		ascii_xyz,
		binary_raw,
		binary_fld
	};

	template <typename T>
	class Field
	{
	protected:
		T *mem = 0;

	public:
		string name = "un_named";
		Frame frame;
		//*construct func
		Field();
		Field(string _name);
		Field(T *_mem_ptr); //*(for fortran mixed coding)
		Field(uint32_t _n_rows, uint32_t _n_cols = 1, uint32_t _n_slices = 1);
		Field(const Frame &_frame);
		//* operator func
		T &operator[](uint64_t i);
		T &operator()(uint32_t i, uint32_t j = 0, uint32_t k = 0);
		T &operator()(uint32_t i, uint32_t j = 0, uint32_t k = 0) const;
		// memory func
		void alloc(uint32_t _n_rows, uint32_t _n_cols = 1, uint32_t _n_slices = 1);
		void alloc(const Frame &_frame);
		void fill(T a);
		void zeros(uint32_t _n_rows, uint32_t _n_cols = 1, uint32_t _n_slices = 1); // alloc()+fill(0)
		void copy(const Field<T> &obj);
		void del();
		//* name
		void set_name(string name);
		string get_name();
		//* math func
		T max();
		T min();
		T abs_sum();
		T abs_max();
		T mean();
		T absmean();
		T *get_mem_p();
		T *&mem_p();
		//* io func
		void save(string file_path, SaveFormat _format);
		void read_binary(string file_path);
		string size_str();

	protected:
		void output_binary_raw(string file_path);
		void output_binary_fld(string file_path);
		void output_ascii_xyz(string file_path);
		void output_ascii_grd(string file_path);
		void output_sgy(string file_path);
		void output_ascii_txt_2d(string file_path);
		void output_ascii_txt_3d(string file_path);
		//
#ifdef JARVIS_DEBUG
	protected:
		MemState mem_state = MemState::inital;

	public:
#ifdef JARVIS_MEMORY_RELEASE_MANAGEMENT
		~Field();
#endif
		MemState get_mem_state() const;
		void set_mem_state(MemState _mem_state);
		void print_frame_info();
#endif
	};
	//******************************************************************
	typedef Field<complex<double>> cx_field;
	typedef Field<complex<double>> cx_dfield;
	typedef Field<complex<float>> cx_ffield;
	typedef Field<double> field;
	typedef Field<double> dfield;
	typedef Field<float> ffield;
	typedef Field<int> ifield;
	//******************************************************************
#define field_for(frame)                         \
	for (int k = 0; k < (frame).n_slices; k++)   \
		for (int j = 0; j < (frame).n_cols; j++) \
			for (int i = 0; i < (frame).n_rows; i++)

	//* time
#define TIC(s)                          \
	clock_t clockBegin##s, clockEnd##s; \
	clockBegin##s = clock();

#define TOC(s, string)                                                               \
	clockEnd##s = clock();                                                           \
	double elapsed_time##s = (double)(clockEnd##s - clockBegin##s) / CLOCKS_PER_SEC; \
	cout << string << elapsed_time##s << "s" << endl;
	//******************************************************************
}
#endif