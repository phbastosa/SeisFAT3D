#ifndef JARVIS_HOST_FRAME_BONES
#define JARVIS_HOST_FRAME_BONES
//**********************************Developer*************************************
// 2021.04.9 BY CAIWEI CALVIN CAI
//********************************************************************************
#include ".jarvis_0_field_header_in.h"
namespace jarvis
{
    enum class FrameType
    {
        null = 0,
        vec = 1,
        mat = 2,
        cube = 3,
        line = 10,
        mesh = 20,
        grid = 30
    };
    struct Frame
    {
        FrameType frame_type = FrameType::null;
        uint32_t n_rows, n_cols, n_slices;
        float d_rows, d_cols, d_slices;
        float l_rows, l_cols, l_slices;
        float r_rows, r_cols, r_slices;
        uint32_t n_elem_slice;
        uint64_t n_elem;
        //
        Frame();
        Frame sparse(uint32_t a);
        void set_n(uint32_t _n_rows, uint32_t _n_cols, uint32_t _n_slices);
        void set_nd(uint32_t _n_rows, float _d_rows);
        void set_nd(uint32_t _n_rows, uint32_t _n_cols, float _d_rows, float _d_cols);
        void set_nd(uint32_t _n_rows, uint32_t _n_cols, uint32_t _n_slices, float _d_rows, float _d_cols, float _d_slices);
        void set_ndl(uint32_t _n_rows, uint32_t _n_cols, uint32_t _n_slices, float _d_rows, float _d_cols, float _d_slices, float _l_rows, float _l_cols, float _l_slices);
        void copy(const Frame &_frame);
        void re_setndl();
        void set_model_type();
        FrameType get_model_type();
        bool check();
        bool is_grid();
        bool is_mesh();
        void print_info(string s = "");
    };
    inline bool frame_compare(Frame a, Frame b);

#define set_frame_n(frame)                                                             \
    int n_rows = (frame).n_rows, n_cols = (frame).n_cols, n_slices = (frame).n_slices; \
    int n_elem_slice = (frame).n_elem_slice

#define set_frame_nd(frame)                                                              \
    int n_rows = (frame).n_rows, n_cols = (frame).n_cols, n_slices = (frame).n_slices;   \
    float d_rows = (frame).d_rows, d_cols = (frame).d_cols, d_slices = (frame).d_slices; \
    int n_elem_slice = (frame).n_elem_slice

#define set_frame_ndl(frame)                                                             \
    int n_rows = (frame).n_rows, n_cols = (frame).n_cols, n_slices = (frame).n_slices;   \
    float d_rows = (frame).d_rows, d_cols = (frame).d_cols, d_slices = (frame).d_slices; \
    float l_rows = (frame).l_rows, l_cols = (frame).l_cols, l_slices = (frame).l_slices; \
    int n_elem_slice = (frame).n_elem_slice

#define set_frame_ndlr(frame)                                                            \
    int n_rows = (frame).n_rows, n_cols = (frame).n_cols, n_slices = (frame).n_slices;   \
    float d_rows = (frame).d_rows, d_cols = (frame).d_cols, d_slices = (frame).d_slices; \
    float l_rows = (frame).l_rows, l_cols = (frame).l_cols, l_slices = (frame).l_slices; \
    float r_rows = (frame).r_rows, r_cols = (frame).r_cols, r_slices = (frame).r_slices; \
    int n_elem_slice = (frame).n_elem_slice

#define set_frame_nd_pre(pre, frame)       \
    int pre##_n_rows = (frame).n_rows,     \
        pre##_n_cols = (frame).n_cols,     \
        pre##_n_slices = (frame).n_slices, \
        pre##_n_elem_slice = (frame).n_elem_slice
}
#endif