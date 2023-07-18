# include "shooting.hpp"

void Shooting::set_modeling_message()
{
    std::cout<<"Running:\n";
    std::cout<<"[0] - Shooting Ray Tracing: Initial Value Problem\n"; 
    std::cout<<"    Cervenỳ (2001)\n\n";
}

void Shooting::initial_setup()
{
    for (int index = 0; index < receiver_output_samples; index++)
    {
        std::vector<float>().swap(ray[index].x);
        std::vector<float>().swap(ray[index].y);
        std::vector<float>().swap(ray[index].z);

        std::vector<float>().swap(ray[index].px);
        std::vector<float>().swap(ray[index].py);
        std::vector<float>().swap(ray[index].pz);

        ray[index].total_time = 0.0f;
    }

    source_id = ((int)(geometry->shots.z[shot_id]/dz) + nbzu) + 
                ((int)(geometry->shots.x[shot_id]/dx) + nbxl)*nzz + 
                ((int)(geometry->shots.y[shot_id]/dy) + nbyl)*nxx*nzz; 
}

void Shooting::forward_solver()
{
    for (float angle = beg_vtangle; angle <= end_vtangle; angle += vtangle_spacing)
    {
        for (float azimuth = beg_azimuth; azimuth <= end_azimuth; azimuth += azimuth_spacing)
        {
            // auxRay.px.push_back(cosf(theta[angle]) * S[]); 
            // auxRay.pz.push_back(sinf(theta[angle]) / v[r[0].z * nx + r[0].x]);    

            // while (true)
            // {


            // }

        }
    }



}

void Shooting::free_space()
{


}