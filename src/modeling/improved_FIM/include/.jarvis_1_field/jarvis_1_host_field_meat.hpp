#ifndef JARVIS_HOST_TFIELD_MEAT
#define JARVIS_HOST_TFIELD_MEAT
//**********************************Developer*************************************
// 2021.04.9 BY CAIWEI CALVIN CAI
//********************************************************************************
#include "jarvis_1_host_field_bones.hpp"
//********************************CLASS_TEMPLATE**********************************
namespace jarvis
{
	template <typename T>
	Field<T>::Field()
	{
		mem = 0;
	}

	template <typename T>
	Field<T>::Field(string _name)
	{
		mem = 0;
		name = _name;
	}

	template <typename T>
	Field<T>::Field(T *_mem_ptr)
	{
		mem = _mem_ptr;
	}

	template <typename T>
	Field<T>::Field(uint32_t _n_rows, uint32_t _n_cols, uint32_t _n_slices)
	{
		frame.set_n(_n_rows, _n_cols, _n_slices);
		mem = (T *)calloc(frame.n_elem, sizeof(T));
#ifdef JARVIS_DEBUG
		mem_state = MemState::allocated;
#endif
	}

	template <typename T>
	inline Field<T>::Field(const Frame &_frame)
	{
#ifdef JARVIS_DEBUG
		if (mem_state == MemState::inital)
		{
#endif
			frame.copy(_frame);
			mem = (T *)calloc(frame.n_elem, sizeof(T));
#ifdef JARVIS_DEBUG
			mem_state = MemState::allocated;
		}
#ifndef JARVIS_NO_WARNING
		else if (mem_state == MemState::allocated)
		{
			printf("Field(frame):[WARNING]:Field [%s] memory has been allocated!\n", name.c_str());
		}
#endif
#endif
	}

	template <typename T>
	void Field<T>::alloc(uint32_t _n_rows, uint32_t _n_cols, uint32_t _n_slices)
	{
#ifdef JARVIS_DEBUG
		if (mem_state == MemState::inital)
		{
#endif
			frame.set_n(_n_rows, _n_cols, _n_slices);
			mem = (T *)calloc(frame.n_elem, sizeof(T));
#ifdef JARVIS_DEBUG
			mem_state = MemState::allocated;
		}
#ifndef JARVIS_NO_WARNING
		else if (mem_state == MemState::allocated)
		{
			printf("alloc():[WARNING]:Field [%s] memory has been allocated!\n", name.c_str());
		}
#endif
#endif
	}

	template <typename T>
	inline void Field<T>::alloc(const Frame &_frame)
	{
#ifdef JARVIS_DEBUG
		if (mem_state == MemState::inital)
		{
#endif
			frame.copy(_frame);
			mem = (T *)calloc(frame.n_elem, sizeof(T));
#ifdef JARVIS_DEBUG
			mem_state = MemState::allocated;
		}
#ifndef JARVIS_NO_WARNING
		else if (mem_state == MemState::allocated)
		{
			printf("Field(frame):[WARNING]:Field [%s] memory has been allocated!\n", name.c_str());
		}
#endif
#endif
	}

	template <typename T>
	void Field<T>::fill(T a)
	{
#ifdef JARVIS_DEBUG
		if (mem_state == MemState::allocated)
		{
#endif
			for (long i = 0; i < frame.n_elem; i++)
				mem[i] = a;
#ifdef JARVIS_DEBUG
		}
		else if (mem_state == MemState::inital)
		{
			printf("fill():\033[41;37m[ERROR]:\033[0mField [%s] memory not allocated!\n", name.c_str());
			std::abort();
		}
		else if (mem_state == MemState::released)
		{
			printf("fill():\033[41;37m[ERROR]:\033[0mField [%s] memory has been released!\n", name.c_str());
			std::abort();
		}
#endif
	}

	template <typename T>
	inline void Field<T>::zeros(uint32_t _n_rows, uint32_t _n_cols, uint32_t _n_slices)
	{
#ifdef JARVIS_DEBUG
		if (mem_state == MemState::inital)
		{
#endif
			uint32_t i;
			frame.set_n(_n_rows, _n_cols, _n_slices);
			mem = (T *)calloc(frame.n_elem, sizeof(T));
			for (i = 0; i < frame.n_elem; i++)
				mem[i] = 0;
#ifdef JARVIS_DEBUG
			mem_state = MemState::allocated;
		}
#ifndef JARVIS_NO_WARNING
		else if (mem_state == MemState::allocated)
		{
			printf("zeros():[WARNING]:Field [%s] memory has been allocated!\n", name.c_str());
			std::abort();
		}
#endif
#endif
	}

	template <typename T>
	inline T *Field<T>::get_mem_p()
	{
		return mem;
	}

	template <typename T>
	inline T *&Field<T>::mem_p()
	{
		return mem;
	}

	template <typename T>
	T Field<T>::max()
	{
		return jarvis::max(mem, frame.n_elem);
	}

	template <typename T>
	T Field<T>::min()
	{
		return jarvis::min(mem, frame.n_elem);
	}

	template <typename T>
	T Field<T>::abs_sum()
	{
		return jarvis::abs_sum(mem, frame.n_elem);
	}

	template <typename T>
	T Field<T>::abs_max()
	{
		return jarvis::abs_max(mem, frame.n_elem);
	}

	template <typename T>
	T Field<T>::absmean()
	{
		return jarvis::abs_mean(mem, frame.n_elem);
	}

	template <typename T>
	T Field<T>::mean()
	{
		return jarvis::mean(mem, frame.n_elem);
	}

	template <typename T>
	void Field<T>::copy(const Field<T> &obj)
	{
#ifdef JARVIS_DEBUG
		mem_state = obj.get_mem_state();
#endif
		frame.frame_type = obj.frame.frame_type;
		frame.copy(obj.frame);
		frame.n_rows = obj.frame.n_rows;
		frame.n_cols = obj.frame.n_cols;
		frame.n_slices = obj.frame.n_slices;
		frame.n_elem_slice = obj.frame.n_elem_slice;
		frame.n_elem = obj.frame.n_elem;
#ifdef JARVIS_DEBUG
		if (mem_state == MemState::allocated)
		{
#endif
			uint32_t i;
			for (i = 0; i < frame.n_elem; i++)
				mem[i] = obj.mem[i];
#ifdef JARVIS_DEBUG
		}
#endif
	}

	template <typename T>
	T &Field<T>::operator[](uint64_t i)
	{
#ifdef JARVIS_DEBUG
		if (mem_state == MemState::allocated && (i >= 0 && i < frame.n_elem))
		{
#endif
			return mem[i];
#ifdef JARVIS_DEBUG
		}
		else if (mem_state == MemState::inital)
		{
			printf("operator[index]:\033[41;37m[ERROR]:\033[0mField [%s] memory not allocated!\n", name.c_str());
			std::abort();
		}
		else if (mem_state == MemState::released)
		{
			printf("operator[index]:\033[41;37m[ERROR]:\033[0mField [%s] memory has been released!\n", name.c_str());
			std::abort();
		}
		else if (!(i >= 0 && i < frame.n_elem))
		{
			printf("operator[index]:\033[41;37m[ERROR]:\033[0mField [%s] index exceeded!\n", name.c_str());
			std::abort();
		}
		else
		{
			printf("operator[index]:\033[41;37m[ERROR]:\033[0mField [%s]  ERROR!\n", name.c_str());
			std::abort();
		}
#endif
	}

	template <typename T>
	T &Field<T>::operator()(uint32_t i, uint32_t j, uint32_t k)
	{
#ifdef JARVIS_DEBUG
		if (mem_state == MemState::allocated && (i >= 0 && i < frame.n_rows) && (j >= 0 && j < frame.n_cols) && (k >= 0 && k < frame.n_slices))
		{
#endif
			return mem[i + j * frame.n_rows + k * frame.n_elem_slice];
#ifdef JARVIS_DEBUG
		}
		else if (mem_state == MemState::inital)
		{
			printf("operator(index):\033[41;37m[ERROR]:\033[0mField [%s] memory not allocated!\n", name.c_str());
			std::abort();
		}
		else if (mem_state == MemState::released)
		{
			printf("operator(index):\033[41;37m[ERROR]:\033[0mField [%s] memory has been released!\n", name.c_str());
			std::abort();
		}
		else if (!(i >= 0 && i < frame.n_rows))
		{
			printf("operator(index):\033[41;37m[ERROR]:\033[0mField [%s] row index exceeded!\n", name.c_str());
			std::abort();
		}
		else if (!(j >= 0 && j < frame.n_cols))
		{
			printf("operator(index):\033[41;37m[ERROR]:\033[0m:Field [%s] col index exceeded!\n", name.c_str());
			std::abort();
		}
		else if (!(k >= 0 && k < frame.n_slices))
		{
			printf("operator(index):\033[41;37m[ERROR]:\033[0m:Field [%s] slice index exceeded!\n", name.c_str());
			std::abort();
		}
		else
		{
			printf("operator(index):\033[41;37m[ERROR]:\033[0mField [%s] ERROR!\n", name.c_str());
			std::abort();
		}
#endif
	}

	template <typename T>
	T &Field<T>::operator()(uint32_t i, uint32_t j, uint32_t k) const
	{
#ifdef JARVIS_DEBUG
		if (mem_state == MemState::allocated && (i >= 0 && i < frame.n_rows) && (j >= 0 && j < frame.n_cols) && (k >= 0 && k < frame.n_slices))
		{
#endif
			return mem[i + j * frame.n_rows + k * frame.n_elem_slice];
#ifdef JARVIS_DEBUG
		}
		else if (mem_state == MemState::inital)
		{
			printf("operator(index):\033[41;37m[ERROR]:\033[0mField [%s] memory not allocated!\n", name.c_str());
			std::abort();
		}
		else if (mem_state == MemState::released)
		{
			printf("operator(index):\033[41;37m[ERROR]:\033[0mField [%s] memory has been released!\n", name.c_str());
			std::abort();
		}
		else if (!(i >= 0 && i < frame.n_rows))
		{
			printf("operator(index):\033[41;37m[ERROR]:\033[0mField [%s] row index exceeded!\n", name.c_str());
			std::abort();
		}
		else if (!(j >= 0 && j < frame.n_cols))
		{
			printf("operator(index):\033[41;37m[ERROR]:\033[0m:Field [%s] col index exceeded!\n", name.c_str());
			std::abort();
		}
		else if (!(k >= 0 && k < frame.n_slices))
		{
			printf("operator(index):\033[41;37m[ERROR]:\033[0m:Field [%s] slice index exceeded!\n", name.c_str());
			std::abort();
		}
		else
		{
			printf("operator(index):\033[41;37m[ERROR]:\033[0mField [%s] ERROR!\n", name.c_str());
			std::abort();
		}
#endif
	}

	template <typename T>
	void Field<T>::save(string file_path, SaveFormat _format)
	{
		if (_format == SaveFormat::binary_raw)
		{
			file_path += size_str() + ".raw";
			output_binary_raw(file_path);
		}
		else if (_format == SaveFormat::binary_fld)
		{
			file_path += ".fld";
			output_binary_fld(file_path);
		}
		else if (_format == SaveFormat::ascii_xyz)
		{
			if (frame.frame_type == FrameType::grid)
			{
				file_path += ".txt";
				output_ascii_xyz(file_path);
			}
			else
			{
				printf("save():\033[41;37m[ERROR]:\033[0m [%s] Format Error! check frame_type!", name.c_str());
				std::abort();
			}
		}
		else if (_format == SaveFormat::ascii_grd)
		{
			if (frame.frame_type == FrameType::mesh)
			{
				file_path += ".grd";
				output_ascii_grd(file_path);
			}
			else
			{
				printf("save():\033[41;37m[ERROR]:\033[0m [%s] Format Error! check frame_type!", name.c_str());
				std::abort();
			}
		}
		else if (_format == SaveFormat::ascii_txt)
		{
			if (frame.frame_type == FrameType::vec || frame.frame_type == FrameType::line || frame.frame_type == FrameType::mat || frame.frame_type == FrameType::mesh)
			{
				file_path += ".txt";
				output_ascii_txt_2d(file_path);
			}
			else if (frame.frame_type == FrameType::cube || frame.frame_type == FrameType::grid)
			{
				file_path += ".txt";
				output_ascii_txt_3d(file_path);
			}
			else
			{
				printf("save():\033[41;37m[ERROR]:\033[0m [%s] Format Error! check frame_type!", name.c_str());
				std::abort();
			}
		}
		else
		{
			printf("save():\033[41;37m[ERROR]:\033[0m [%s] Format Error!", name.c_str());
			std::abort();
		}
	}

	template <typename T>
	void Field<T>::read_binary(string file_path)
	{
		string FileEndStr = file_path.substr(file_path.length() - 4, file_path.length());
		string EndStr = ".fld";
		if (FileEndStr == EndStr)
		{
			FILE *fp;
			fp = fopen(file_path.c_str(), "rb");
			if (fp == 0)
			{
				printf("read_binary():\033[41;37m[ERROR]:\033[0m File [%s] open error!", name.c_str());
				std::abort();
			}
			fread(&frame, sizeof(Frame), 1, fp);
			mem = (T *)calloc(frame.n_elem, sizeof(T));
			fread(mem, frame.n_elem * sizeof(T), 1, fp);
			fclose(fp);
		}
		else
		{
			printf("read_binary:\033[41;37m[ERROR]:\033[0m [%s] Can't read field-raw file %s ,exit\n", file_path.c_str());
			std::abort();
		}
	}

	template <typename T>
	void Field<T>::output_binary_raw(string file_path)
	{
		FILE *fp = fopen(file_path.c_str(), "wb");
		fwrite(mem, frame.n_elem * sizeof(T), 1, fp);
		fclose(fp);
		printf("\033[45;30m[save Field-raw File]:\033[0m%s\n", file_path.c_str());
	}

	template <typename T>
	void Field<T>::output_binary_fld(string file_path)
	{
		FILE *fp = fopen(file_path.c_str(), "wb");
		fwrite(&frame, sizeof(Frame), 1, fp);
		fwrite(mem, frame.n_elem * sizeof(T), 1, fp);
		printf("\033[45;30m[save Field-fld File]:\033[0m%s\n", file_path.c_str());
	}

	template <typename T>
	void Field<T>::output_ascii_xyz(string file_path)
	{
		FILE *fptime;
		const char *FilePathChar = file_path.c_str();
		printf("save X-Y-Z-V File:%s\n", file_path.c_str());
		if (FilePathChar != NULL)
		{
			if ((fptime = fopen(FilePathChar, "w")) == NULL)
			{
				printf("grid_output:\033[41;37m[ERROR]:\033[0mCan't create output file %s ,exit\n", FilePathChar);
				std::abort();
			}

			if (typeid(mem[0]).name()[0] == 'd' || typeid(mem[0]).name()[0] == 'f')
			{
				double x, y, z;
				for (int idx = 0; idx < frame.n_elem; idx++)
				{
					int i = idx % frame.n_elem_slice % frame.n_rows;
					int j = idx % frame.n_elem_slice / frame.n_rows;
					int k = idx / frame.n_elem_slice;
					x = frame.l_rows + i * frame.d_rows;
					y = frame.l_cols + j * frame.d_cols;
					z = frame.l_slices + k * frame.d_slices;
					fprintf(fptime, "%7.3f,%7.3f,%7.3f,%8.6e\n", x, y, z, mem[idx]);
				}
			}
			else if (typeid(mem[0]).name()[0] == 'i')
			{
				double x, y, z;
				for (int idx = 0; idx < frame.n_elem; idx++)
				{
					int i = idx % frame.n_elem_slice % frame.n_rows;
					int j = idx % frame.n_elem_slice / frame.n_rows;
					int k = idx / frame.n_elem_slice;
					x = frame.l_rows + i * frame.d_rows;
					y = frame.l_cols + j * frame.d_cols;
					z = frame.l_slices + k * frame.d_slices;
					fprintf(fptime, "%7.3f,%7.3f,%7.3f,%d\n", x, y, z, mem[idx]);
				}
			}
		}
		fclose(fptime);
	}

	template <typename T>
	void Field<T>::output_ascii_grd(string file_path)
	{
		string EndStr = file_path.substr(file_path.length() - 4, file_path.length());
		string EndStr1 = ".grd";
		if (EndStr == EndStr1)
		{
			cout << "save .grd File:" << file_path << endl;
			ofstream fout(file_path);
			if (fout.is_open())
			{
				fout << "DSAA" << endl;
				fout << frame.n_rows << " " << frame.n_cols << endl;
				fout << frame.l_rows << " " << frame.r_rows << endl;
				fout << frame.l_cols << " " << frame.r_cols << endl;
				fout << jarvis::min(mem, frame.n_elem) << " " << jarvis::max(mem, frame.n_elem) << endl;
				for (int j = 0; j < frame.n_cols; j++)
				{
					for (int i = 0; i < frame.n_rows; i++)
					{
						fout << mem[i + j * frame.n_rows] << " ";
					}
					fout << endl;
				}
				fout.close();
			}
			else
			{
				printf("MeshOutPut():\033[41;37m[ERROR]:\033[0mCan't create GRD output file %s ,exit\n", file_path.c_str());
				std::abort();
			}
		}
		else
		{
			printf("save X-Y-V File:%s\n", file_path.c_str());
			ofstream fout(file_path);
			if (fout.is_open())
			{
				for (int i = 0; i < frame.n_rows; i++)
				{
					for (int j = 0; j < frame.n_cols; j++)
					{
						fout << frame.l_rows + i * frame.d_rows << " " << frame.l_cols + j * frame.d_cols << " " << mem[i + j * frame.n_rows] << endl;
					}
				}
				fout.close();
			}
			else
			{
				printf("MeshOutPut():\033[41;37m[ERROR]:\033[0mCan't create XYV output file %s ,exit\n", file_path.c_str());
				std::abort();
			}
		}
	}

	template <typename T>
	void Field<T>::output_ascii_txt_2d(string file_path)
	{
		printf("save TXT File:%s\n", file_path.c_str());
		ofstream fout(file_path);
		if (fout.is_open())
		{
			for (int i = 0; i < frame.n_rows; i++)
			{
				for (int j = 0; j < frame.n_cols; j++)
				{
					fout << mem[i + j * frame.n_rows] << " ";
				}
				fout << endl;
			}
			fout.close();
		}
		else
		{
			printf("MeshOutPut():\033[41;37m[ERROR]:\033[0mCan't create TXT output file %s ,exit\n", file_path.c_str());
			std::abort();
		}
	}

	template <typename T>
	void Field<T>::output_ascii_txt_3d(string file_path)
	{
		printf("save TXT File:%s\n", file_path.c_str());
		FILE *fp;
		if (fp = fopen(file_path.c_str(), "wt"))
		{
			uint32_t i, j, k;
			for (k = 0; k < frame.n_slices; k++)
			{
				for (j = 0; j < frame.n_cols; j++)
				{
					for (i = 0; i < frame.n_rows; i++)
					{
						fprintf(fp, "%d,%d,%d,%10.6e\n", i, j, k, mem[i + j * frame.n_rows + k * frame.n_elem_slice]);
					}
				}
			}
			fclose(fp);
		}
		else
		{
			printf("CubeOutPut():\033[41;37m[ERROR]:\033[0mCan't create TXT output file %s ,exit\n", file_path.c_str());
			std::abort();
		}
	}

	template <typename T>
	string Field<T>::size_str()
	{
		return "_" + to_string(frame.n_rows) + "_" + to_string(frame.n_cols) + "_" + to_string(frame.n_slices) +
			   "_" + to_string(frame.l_rows) + "_" + to_string(frame.l_cols) + "_" + to_string(frame.l_slices);
	}

	template <typename T>
	void Field<T>::output_sgy(string file_path)
	{
	}

	template <typename T>
	void Field<T>::set_name(string _name)
	{
		name = _name;
	}
	template <typename T>
	string Field<T>::get_name()
	{
		return name;
	}

#ifdef JARVIS_DEBUG
#ifdef JARVIS_MEMORY_RELEASE_MANAGEMENT
	template <typename T>
	Field<T>::~Field()
	{
		if (mem_state == MemState::allocated)
		{
			printf("Field [%s] Memory has not been released!", name.c_str());
			std::abort();
		}
	}
#endif
	template <typename T>
	MemState Field<T>::get_mem_state() const
	{
		return mem_state;
	}

	template <typename T>
	void Field<T>::set_mem_state(MemState _mem_state)
	{
		mem_state = _mem_state;
	}

	template <typename T>
	inline void Field<T>::print_frame_info()
	{
		cout << "FrameType:" << uint32_t(frame.frame_type) << endl;
		cout << "Frame_X_Num:" << setw(5) << frame.n_rows << "; INTERVAL:" << setw(6) << frame.d_rows << "; RANGE:"
			 << "[" << frame.l_rows << ", " << frame.r_rows << "]" << endl;
		cout << "Frame_Y_Num:" << setw(5) << frame.n_cols << "; INTERVAL:" << setw(6) << frame.d_cols << "; RANGE:"
			 << "[" << frame.l_cols << ", " << frame.r_cols << "]" << endl;
		cout << "Frame_Z_Num:" << setw(5) << frame.n_slices << "; INTERVAL:" << setw(6) << frame.d_slices << "; RANGE:"
			 << "[" << frame.l_slices << ", " << frame.r_slices << "]" << endl;
	}
#endif
	//
	template <typename T>
	void Field<T>::del()
	{
#ifdef JARVIS_DEBUG
		if (mem_state == MemState::allocated)
		{
#endif
			free(mem);
			mem = 0;
#ifdef JARVIS_DEBUG
			mem_state = MemState::released;
		}
#ifndef JARVIS_NO_WARNING
		else if (mem_state == MemState::inital)
		{
			printf("del():[WARNING]:Field [%s] memory is not allocated and does not need to be released!\n", name.c_str());
		}
		else if (mem_state == MemState::released)
		{
			printf("del():[WARNING]:Field [%s] does not need to be released again!\n", name.c_str());
		}
#endif
#endif
	}
}
#endif