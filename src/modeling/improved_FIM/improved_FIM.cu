# include "improved_FIM.cuh"

void Improved_FIM::set_specific_boundary()
{
    int pdx = (BLOCK_LENGTH - nx % BLOCK_LENGTH) % BLOCK_LENGTH;
    int pdy = (BLOCK_LENGTH - ny % BLOCK_LENGTH) % BLOCK_LENGTH;
    int pdz = (BLOCK_LENGTH - nz % BLOCK_LENGTH) % BLOCK_LENGTH;    
    
    if (pdx == 0)      { nbxl = 2; nbxr = 2; }
    else if (pdx == 1) { nbxl = 3; nbxr = 2; }
    else if (pdx == 2) { nbxl = 3; nbxr = 3; }
    else if (pdx == 3) { nbxl = 4; nbxr = 3; }

    if (pdy == 0)      { nbyl = 2; nbyr = 2; }
    else if (pdy == 1) { nbyl = 3; nbyr = 2; }
    else if (pdy == 2) { nbyl = 3; nbyr = 3; }
    else if (pdy == 3) { nbyl = 4; nbyr = 3; }

    if (pdz == 0)      { nbzu = 2; nbzd = 2; }
    else if (pdz == 1) { nbzu = 3; nbzd = 2; }
    else if (pdz == 2) { nbzu = 3; nbzd = 3; }
    else if (pdz == 3) { nbzu = 4; nbzd = 3; }
}

void Improved_FIM::set_eikonal_volumes()
{
    eikonal_method = std::string("hd_fim");
	eikonal_message = std::string("[3] - Improved FIM (Cai, Zhu & Li, 2023)");
    
    T = new float[volsize](); 
}

void Improved_FIM::initialization()
{

}

void Improved_FIM::forward_solver()
{
    Frame vel_grid;

    vel_grid.set_nd(nxx, nyy, nzz, 25, 25, 25);

    fcufld Vel(MemType::npin, vel_grid);

    for (int k = 0; k < nzz; k++)
        for (int j = 0; j < nyy; j++)
            for (int i = 0; i < nxx; i++)
                Vel(i,j,k) = 1.0f / S[k + i*nzz + j*nxx*nzz];    

    Vel.cu_copy_h2d();

    tp3cuvec source(MemType::pin, 1);

    source(0).x = geometry->shots.x[shot_id] + dx*nbxl; 
    source(0).y = geometry->shots.y[shot_id] + dy*nbyl;
    source(0).z = geometry->shots.z[shot_id] + dz*nbzu;

    source(0).time = S[source_id] * sqrtf(powf((sidx-nbxl)*dx - geometry->shots.x[shot_id], 2.0f) + powf((sidy-nbyl)*dy - geometry->shots.y[shot_id], 2.0f) + powf((sidz-nbzu)*dz - geometry->shots.z[shot_id], 2.0f));    

    travel_time_3d_2rd_diag.module_init(&Vel, 22); // second-order difference with diagonal node

    travel_time_3d_2rd_diag.cal_travel_time(source, TravelType::normal);
    // travel_time_3d_2rd_diag.cal_travel_time(source, TravelType::refine); It's not working!

    travel_time_3d_2rd_diag.get_device_time();
    
    for (int index = 0; index < travel_time_3d_2rd_diag.time.frame.n_elem; index++)
    {
        int k = (int) (index / (nxx*nzz));         // y direction
        int j = (int) (index - k*nxx*nzz) / nzz;   // x direction
        int i = (int) (index - j*nzz - k*nxx*nzz); // z direction        
    
        T[i + j*nzz + k*nxx*nzz] = travel_time_3d_2rd_diag.time[j + k*nxx + i*nxx*nyy];
    }
}

namespace jarvis
{
	void travel_time_3d_module::cu_cal_time(tp3cuvec &_shotline)
	{
		set_source(_shotline);
		do
		{
			endflag[0] = true;
			endflag.cu_stream_copy_h2d(jarvis_default_stream);
			fim_kernel_0_cal_active_node<<<jarvis_cuda_kernel_size(vel_p->frame.n_elem)>>>(vel_p->frame, vel_p->cu_mem, time.cu_mem, mark.cu_mem, diff_order);
			fim_kernel_1_mark_new_active_node<<<jarvis_cuda_kernel_size(vel_p->frame.n_elem)>>>(vel_p->frame, vel_p->cu_mem, time.cu_mem, mark.cu_mem, diff_order);
			fim_kernel_2_set_as_converged_node<<<jarvis_cuda_kernel_size(vel_p->frame.n_elem)>>>(vel_p->frame, mark.cu_mem);
			fim_kernel_3_check_is_finishd<<<jarvis_cuda_kernel_size(vel_p->frame.n_elem)>>>(vel_p->frame, mark.cu_mem, endflag.cu_mem);
			endflag.cu_stream_copy_d2h(jarvis_default_stream);
			cudaStreamSynchronize(jarvis_default_stream);
		} while (endflag[0] == false);
		endflag[0] == true;
	}
}

void Improved_FIM::free_space()
{

}