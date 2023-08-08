# include "inversion.hpp"

// void Inversion::set_parameters()
// {
//     Modeling::file = file;
//     Modeling::set_parameters();

//     iteration = 0;

//     n_model = nPoints;
//     n_data = total_nodes * total_shots;

//     tolerance = std::stof(catch_parameter("tolerance", file));
//     max_iteration = std::stof(catch_parameter("max_iteration", file));

//     obs_data_folder = catch_parameter("obs_data_folder", file);
//     obs_data_prefix = catch_parameter("obs_data_prefix", file);

//     dobs = new float[n_data]();    
//     dcal = new float[n_data]();    

//     source = new float[volsize]();
//     adjoint = new float[volsize]();
    
//     gradient = new float[n_model]();

//     cudaMalloc((void**)&(d_source), volsize*sizeof(float));
//     cudaMalloc((void**)&(d_adjoint), volsize*sizeof(float));
// }

// void Inversion::import_obs_data()
// {    
//     int ptr = 0; 

//     float * data = new float[total_nodes]();

//     for (int shot = 0; shot < total_shots; shot++)
//     {
//         import_binary_float(obs_data_folder + obs_data_prefix + std::to_string(shot+1) + ".bin", data, total_nodes);

//         for (int node = ptr; node < ptr + total_nodes; node++) 
//             dobs[node] = data[node - ptr];

//         ptr += total_nodes;        
//     }

//     delete[] data;

//     set_runtime();
// }

// void Inversion::forward_modeling()
// {
//     set_slowness();

//     for (int shot = 0; shot < total_shots; shot++)
//     {
//         shot_id = shot;

//         info_message();

//         initial_setup();

//         forward_solver();

//         build_outputs();

//         fill_calculated_data();

//         adjoint_state_solver();



//     }


// }

// void Inversion::fill_calculated_data()
// {
//     for (int node = 0; node < total_nodes; node++)
//     {
//         int index = node + shot_id * total_nodes;

//         dcal[index] = receiver_output[node]; 
//     }
// }

// void Inversion::adjoint_state_solver()
// {
//     float cell_volume = dx * dy * dz;

//     for (int index = 0; index < volsize; index++) 
//     {
//         source[index] = 0.0f;    
//         adjoint[index] = 1e6f;

//         int k = (int) (index / (nxx*nzz));        
//         int j = (int) (index - k*nxx*nzz) / nzz;    
//         int i = (int) (index - j*nzz - k*nxx*nzz);  

//         if ((i == 0) || (i == nzz-1) || 
//             (j == 0) || (j == nxx-1) || 
//             (k == 0) || (k == nyy-1))  
//         {    
//             adjoint[index] = 0.0f;        
//         }
//     }

//     for (int node = 0; node < total_nodes; node++)
//     {
//         int node_id = node + shot_id*total_nodes;

//         float tmp = dcal[node_id] + t0 - dobs[node_id]; 

//         int i = (int)(geometry->nodes.z[node] / dz);
//         int j = (int)(geometry->nodes.x[node] / dx);
//         int k = (int)(geometry->nodes.y[node] / dy);

//         int index = i + j*nzz + k*nzz*nxx;

//         source[index] += tmp / cell_volume; 
//         source[index + 1] += tmp / cell_volume; 
//         source[index + nzz] += tmp / cell_volume; 
//         source[index + 1 + nzz] += tmp / cell_volume; 
//         source[index + nxx*nzz] += tmp / cell_volume; 
//         source[index + 1 + nxx*nzz] += tmp / cell_volume; 
//         source[index + nzz + nxx*nzz] += tmp / cell_volume; 
//         source[index + 1 + nzz + nxx*nzz] += tmp / cell_volume; 
//     }

// 	cudaMemcpy(d_source, source, volsize*sizeof(float), cudaMemcpyHostToDevice);
// 	cudaMemcpy(d_adjoint, adjoint, volsize*sizeof(float), cudaMemcpyHostToDevice);

//     for (int sweep = 0; sweep < nSweeps; sweep++)
// 	{ 
// 		int start = (sweep == 3 || sweep == 5 || sweep == 6 || sweep == 7) ? totalLevels : meshDim;
// 		int end = (start == meshDim) ? totalLevels + 1 : meshDim - 1;
// 		int incr = (start == meshDim) ? true : false;

// 		int xSweepOff = (sweep == 3 || sweep == 4) ? nxx : 0;
// 		int ySweepOff = (sweep == 2 || sweep == 5) ? nyy : 0;
// 		int zSweepOff = (sweep == 1 || sweep == 6) ? nzz : 0;
		
// 		for (int level = start; level != end; level = (incr) ? level + 1 : level - 1)
// 		{			
// 			int xs = max(1, level - (nyy + nzz));	
// 			int ys = max(1, level - (nxx + nzz));

// 			int xe = min(nxx, level - (meshDim - 1));
// 			int ye = min(nyy, level - (meshDim - 1));	
		
// 			int xr = xe - xs + 1;
// 			int yr = ye - ys + 1;

// 			int nThreads = xr * yr;
				
// 			dim3 bs(16, 16, 1);

// 			if (nThreads < threadsPerBlock) { bs.x = xr; bs.y = yr; } 

// 			dim3 gs(iDivUp(xr, bs.x), iDivUp(yr , bs.y), 1);

//             adjoint_state_kernel<<<gs,bs>>>(d_adjoint, d_source, d_T, level, xs, ys, xSweepOff, 
//                                             ySweepOff, zSweepOff, nxx, nyy, nzz, dx, dy, dz);

//             cudaDeviceSynchronize();
// 		}
// 	}

//     cudaMemcpy(adjoint, d_adjoint, volsize*sizeof(float), cudaMemcpyDeviceToHost);

//     for (int index = 0; index < n_model; index++) 
//     {
//         int k = (int) (index / (nx*nz));        
//         int j = (int) (index - k*nx*nz) / nz;    
//         int i = (int) (index - j*nz - k*nx*nz);  

//         int indp = (i+padb) + (j+padb)*nzz + (k+padb)*nxx*nzz;

//         if (T[indp] <= 0.3) 
//             gradient[index] += 0.01*T[indp]*adjoint[indp]*S[indp]*S[indp]*cell_volume;
//         else 
//             gradient[index] += adjoint[indp]*S[indp]*S[indp]*cell_volume;

//         if ((i == 0) || (i == nz-1) || 
//             (j == 0) || (j == nx-1) || 
//             (k == 0) || (k == ny-1))  
//         {    
//             gradient[index] = 0.0f;        
//         }
//     }
// }

// void Inversion::compute_residuals()
// {
//     float r = 0.0f;
//     for (int i = 0; i < n_data; i++)
//         r += powf(dobs[i] - dcal[i], 2.0f);

//     residuo.push_back(sqrtf(r));
// }

// void Inversion::check_convergence()
// {
//     if ((iteration >= max_iteration) || (residuo.back() <= tolerance))
//     {
//         std::cout<<"\nFinal residuo: "<<residuo.back()<<std::endl;
//         converged = true;
//     }
//     else
//     {
//         iteration += 1;
//         converged = false;
//     }
// }

// void Inversion::optimization()
// {


// }

// void Inversion::model_update()
// {



// }

// void Inversion::export_results()
// {
//     get_runtime();

//     export_binary_float("gradient.bin", gradient, n_model);
// }

// __global__ void adjoint_state_kernel(float * adjoint, float * source, float * T, int level, int xOffset, int yOffset, 
//                                      int xSweepOffset, int ySweepOffset, int zSweepOffset, int nxx, int nyy, int nzz, 
//                                      float dx, float dy, float dz)
// {
// 	int x = (blockIdx.x * blockDim.x + threadIdx.x) + xOffset;
// 	int y = (blockIdx.y * blockDim.y + threadIdx.y) + yOffset;

// 	if ((x <= nxx) && (y <= nyy)) 
// 	{
// 		int z = level - (x + y);
		
// 		if ((z > 0) && (z <= nzz))	
// 		{
// 			int i = abs(z - zSweepOffset);
// 			int j = abs(x - xSweepOffset);
// 			int k = abs(y - ySweepOffset);

// 			if ((i > 0) && (i < nzz-1) && (j > 0) && (j < nxx-1) && (k > 0) && (k < nyy-1))
// 			{		
//                 float a1 = -1.0f*(T[i + j*nzz + k*nxx*nzz] - T[i + (j-1)*nzz + k*nxx*nzz]) / dx;
//                 float ap1 = (a1 + abs(a1)) / 2.0f;
//                 float am1 = (a1 - abs(a1)) / 2.0f;

//                 float a2 = -1.0f*(T[i + (j+1)*nzz + k*nxx*nzz] - T[i + j*nzz + k*nxx*nzz]) / dx;
//                 float ap2 = (a2 + abs(a2)) / 2.0f;
//                 float am2 = (a2 - abs(a2)) / 2.0f;

//                 float b1 = -1.0f*(T[i + j*nzz + k*nxx*nzz] - T[i + j*nzz + (k-1)*nxx*nzz]) / dy;
//                 float bp1 = (b1 + abs(b1)) / 2.0f;
//                 float bm1 = (b1 - abs(b1)) / 2.0f;

//                 float b2 = -1.0f*(T[i + j*nzz + (k+1)*nxx*nzz] - T[i + j*nzz + k*nxx*nzz]) / dy;
//                 float bp2 = (b2 + abs(b2)) / 2.0f;
//                 float bm2 = (b2 - abs(b2)) / 2.0f;

//                 float c1 = -1.0f*(T[i + j*nzz + k*nxx*nzz] - T[(i-1) + j*nzz + k*nxx*nzz]) / dz;
//                 float cp1 = (c1 + abs(c1)) / 2.0f;
//                 float cm1 = (c1 - abs(c1)) / 2.0f;

//                 float c2 = -1.0f*(T[(i+1) + j*nzz + k*nxx*nzz] - T[i + j*nzz + k*nxx*nzz]) / dz;
//                 float cp2 = (c2 + abs(c2)) / 2.0f;
//                 float cm2 = (c2 - abs(c2)) / 2.0f;

//                 float d = (ap2 - am1)/dx + (bp2 - bm1)/dy + (cp2 - cm1)/dz;

//                 if (d < 1e-6f)
//                 {
//                     adjoint[i + j*nzz + k*nxx*nzz] = 0.0f;    
//                 }
//                 else
//                 {
//                     float e = (ap1*adjoint[i + (j-1)*nzz + k*nxx*nzz] - am2*adjoint[i + (j+1)*nzz + k*nxx*nzz]) / dx +
//                               (bp1*adjoint[i + j*nzz + (k-1)*nxx*nzz] - bm2*adjoint[i + j*nzz + (k+1)*nxx*nzz]) / dy +
//                               (cp1*adjoint[(i-1) + j*nzz + k*nxx*nzz] - cm2*adjoint[(i+1) + j*nzz + k*nxx*nzz]) / dz;

//                     float f = (e + source[i + j*nzz + k*nxx*nzz]) / d;
            
//                     adjoint[i + j*nzz + k*nxx*nzz] = min(adjoint[i + j*nzz + k*nxx*nzz], f);
//                 }
//             }
//         }
//     }
// }