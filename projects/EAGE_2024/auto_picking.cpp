# include <cmath>
# include <chrono>
# include <string>
# include <fstream>
# include <iostream>
# include <algorithm>

# include "../../src/io/io.hpp"

int main(int argc, char**argv)
{
    auto ti = std::chrono::system_clock::now();

    int nt = 6001;
    int total_shots = 10 * 10; 
    int total_nodes = 51 * 51;

    float dt = 1e-3f;
    
    float window = 0.01f;
    float amp_cut = 1e-15f;

    std::string data_folder = "../inputs/data/";    

    int iw = (int)(window/dt) + 1;

    float * picks = new float[total_nodes]();    
    float * gather = new float[total_nodes * nt]();

    for (int shot = 0; shot < total_shots; shot++)
    {
        import_binary_float(data_folder + "synthetic_seismogram_" + std::to_string(nt) + "x" + std::to_string(total_nodes) + "_shot_" + std::to_string(shot+1) + ".bin", gather, total_nodes * nt);

        float gatherMaxAmp = *std::max_element(gather, gather + total_nodes * nt);

        for (int i = 0; i < total_nodes * nt; i++)
            gather[i] *= 1.0f / gatherMaxAmp;

        std::cout<<"Running shot " + std::to_string(shot+1) + " of " + std::to_string(total_shots) << std::endl;

        for (int trace = 0; trace < total_nodes; trace++)
        {
            if (trace % (total_nodes / 10) == 0)
            {
                std::cout<<"    Running trace " + std::to_string(trace) + " of " + std::to_string(total_nodes) << std::endl;
            }

            float * A = new float[nt]();
            float * B = new float[nt]();
            float * S = new float[nt](); 

            for (int time = iw; time < nt - iw; time++)
            {
                for (int k = 0; k < iw; k++)
                {
                    A[time] += gather[(time - k) + trace*nt] + 1e-15f;
                    B[time] += gather[(time + k) + trace*nt] + 1e-15f;
                }
                
                S[time] = fabsf(B[time] / A[time] * (B[time] - A[time]));
                S[time] *= (B[time] + A[time]) * (B[time] - A[time]) / powf(time + iw, 2.0f);
            }

            float ampMax = *std::max_element(S, S+nt);

            for (int k = 0; k < nt; k++)
            {
                S[k] *= 1.0f / ampMax;
                
                if (S[k] > amp_cut)
                {
                    picks[trace] = k * dt + 0.045f;            
                    break;
                }
            }    

            delete[] A;
            delete[] B;
            delete[] S;
        }

        export_binary_float(data_folder + "rawPicks_" + std::to_string(total_nodes) + "_shot_"+ std::to_string(shot+1) + ".bin", picks, total_nodes);
    }

    delete[] picks;
    delete[] gather;

    auto tf = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = tf - ti;

    std::cout << "\nRun time: " << elapsed_seconds.count() << " s\n";

    return 0;
}