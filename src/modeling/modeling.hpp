# ifndef MODELING_HPP
# define MODELING_HPP

# include <chrono>

# include "geometry/geometry.hpp"
# include "geometry/regular/regular.hpp"
# include "geometry/circular/circular.hpp"
# include "geometry/streamer/streamer.hpp"

class Modeling
{
private:

    std::chrono::_V2::system_clock::time_point ti, tf;

protected:

    int receiver_output_samples;
    int wavefield_output_samples;

    bool export_receiver_output;
    bool export_wavefield_output;

    std::string receiver_output_file;
    std::string wavefield_output_file;
    std::string receiver_output_folder;
    std::string wavefield_output_folder;

    float * receiver_output = nullptr;
    float * wavefield_output = nullptr;

    int nx, ny, nz;
    float dx, dy, dz;

    float * S = nullptr;
    float * V = nullptr;

public: 

    int total_shots;
    int total_nodes;

    virtual void set_parameters(std::string file); 

    // virtual void set_components() = 0;

    // virtual void initial_setup() = 0;
    // virtual void forward_solver() = 0;

    void set_runtime();
    void get_runtime();

    void export_outputs();
};

# endif