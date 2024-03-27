# include "eikonal.hpp"

// void Eikonal::set_parameters()
// {
//     set_general_parameters();

//     set_acquisition_geometry();

//     set_velocity_model();
    
//     set_specific_boundary();         
//     set_boundaries();
    
//     set_slowness_model();
//     set_outputs();

//     set_eikonal_volumes();
// }

// void Eikonal::set_runtime()
// {
//     ti = std::chrono::system_clock::now();
// }

// void Eikonal::info_message()
// {
//     get_RAM_usage();
//     get_GPU_usage();

//     auto clear = system("clear");
        
//     std::cout<<"Model dimensions (z = "<<(nz-1)*dz<<", x = "<<(nx-1)*dx<<", y = "<<(ny-1)*dy<<") m\n\n";

//     std::cout<<"Shot "<<shot_id+1<<" of "<<geometry->shots.total;

//     std::cout<<" at position (z = "<<geometry->shots.z[shot_id]<<", x = " 
//                                    <<geometry->shots.x[shot_id]<<", y = " 
//                                    <<geometry->shots.y[shot_id]<<") m\n\n";

//     std::cout<<"Memory usage: \n";
//     std::cout<<"RAM = "<<RAM<<" Mb\n";
//     std::cout<<"GPU = "<<vRAM<<" Mb\n\n";

//     std::cout<<"Eikonal solver:\n";
//     std::cout<<eikonal_message<<"\n\n";
// }

// void Eikonal::initial_setup()
// {
//     sidx = (int)(geometry->shots.x[shot_id] / dx) + nbxl;
//     sidy = (int)(geometry->shots.y[shot_id] / dy) + nbyl;
//     sidz = (int)(geometry->shots.z[shot_id] / dz) + nbzu;

//     source_id = sidz + sidx*nzz + sidy*nxx*nzz;

//     initialization();
// }

// void Eikonal::build_outputs()
// {
//     get_travel_times();
//     get_first_arrivals();
// }

// void Eikonal::export_outputs()
// {
//     if (export_receiver_output) 
//         export_binary_float(receiver_output_file, receiver_output, receiver_output_samples);
    
//     if (export_wavefield_output) 
//         export_binary_float(wavefield_output_file, wavefield_output, wavefield_output_samples);
// }

// void Eikonal::get_runtime()
// {
//     tf = std::chrono::system_clock::now();

//     std::chrono::duration<double> elapsed_seconds = tf - ti;

//     std::ofstream runTimeFile("elapsedTime.txt",std::ios::in | std::ios::app);
//     runTimeFile << "# Eikonal " + eikonal_method + " - Spacing: " + std::to_string(dx) + " m\n";
//     runTimeFile << "# Run Time [s]; RAM usage [MB]; GPU memory usage [MB]\n";
//     runTimeFile << std::to_string(elapsed_seconds.count()) + ";" + std::to_string(RAM) + ";" + std::to_string(vRAM) + "\n";
//     runTimeFile.close();

//     std::cout<<"\nRun time: "<<elapsed_seconds.count()<<" s."<<std::endl;
// }

// void Eikonal::expand_boundary(float * input, float * output)
// {
//     for (int y = nbyl; y < nyy - nbyr; y++)
//     {
//         for (int x = nbxl; x < nxx - nbyr; x++)
//         {
//             for (int z = nbzu; z < nzz - nbzd; z++)
//             {
//                 output[z + x*nzz + y*nxx*nzz] = input[(z - nbzu) + (x - nbxl)*nz + (y - nbyl)*nx*nz];       
//             }
//         }
//     }

//     for (int y = nbyl; y < nyy - nbyr; y++)
//     {
//         for (int x = nbxl; x < nxx - nbyr; x++)
//         {
//             for (int z = 0; z < nbzu; z++)
//             {
//                 output[z + x*nzz + y*nxx*nzz] = input[0 + (x - nbxl)*nz + (y - nbyl)*nx*nz];
//             }
//         }
//     }

//     for (int y = nbyl; y < nyy - nbyr; y++)
//     {
//         for (int x = nbxl; x < nxx - nbyr; x++)
//         {
//             for (int z = nzz - nbzd; z < nzz; z++)
//             {
//                 output[z + x*nzz + y*nxx*nzz] = input[(nz - 1) + (x - nbxl)*nz + (y - nbyl)*nx*nz];
//             }
//         }
//     }

//     for (int y = nbyl; y < nyy - nbyr; y++)
//     {
//         for (int x = nbxl; x < nxx - nbyr; x++)
//         {
//             for (int z = nzz - nbzd; z < nzz; z++)
//             {
//                 output[z + x*nzz + y*nxx*nzz] = input[(nz - 1) + (x - nbxl)*nz + (y - nbyl)*nx*nz];
//             }
//         }
//     }

//     for (int y = nbyl; y < nyy - nbyr; y++)
//     {
//         for (int x = 0; x < nbxl; x++)
//         {
//             for (int z = 0; z < nzz; z++)
//             {
//                 output[z + x*nzz + y*nxx*nzz] = output[z + nbxl*nzz + y*nxx*nzz];
//             }
//         }
//     }

//     for (int y = nbyl; y < nyy - nbyr; y++)
//     {
//         for (int x = nxx-nbxr; x < nxx; x++)
//         {
//             for (int z = 0; z < nzz; z++)
//             {
//                 output[z + x*nzz + y*nxx*nzz] = output[z + (nxx - nbxr - 1)*nzz + y*nxx*nzz];
//             }
//         }
//     }

//     for (int y = 0; y < nbyl; y++)
//     {
//         for (int x = 0; x < nxx; x++)
//         {
//             for (int z = 0; z < nzz; z++)
//             {
//                 output[z + x*nzz + y*nxx*nzz] = output[z + x*nzz + nbyl*nxx*nzz];
//             }
//         }
//     }

//     for (int y = nyy - nbyr; y < nyy; y++)
//     {
//         for (int x = 0; x < nxx; x++)
//         {
//             for (int z = 0; z < nzz; z++)
//             {
//                 output[z + x*nzz + y*nxx*nzz] = output[z + x*nzz + (nyy - nbyr - 1)*nxx*nzz];
//             }
//         }
//     }
// }

// void Eikonal::reduce_boundary(float * input, float * output)
// {
//     for (int index = 0; index < nPoints; index++)
//     {
//         int y = (int) (index / (nx*nz));         
//         int x = (int) (index - y*nx*nz) / nz;    
//         int z = (int) (index - x*nz - y*nx*nz);  

//         output[z + x*nz + y*nx*nz] = input[(z + nbzu) + (x + nbxl)*nzz + (y + nbyl)*nxx*nzz];
//     }
// }

// void Eikonal::set_general_parameters()
// {
//     get_GPU_initMem();

//     threadsPerBlock = 256;

//     nx = std::stoi(catch_parameter("x_samples", file));
//     ny = std::stoi(catch_parameter("y_samples", file));
//     nz = std::stoi(catch_parameter("z_samples", file));
    
//     dx = std::stof(catch_parameter("x_spacing", file));
//     dy = std::stof(catch_parameter("y_spacing", file));
//     dz = std::stof(catch_parameter("z_spacing", file));

//     nPoints = nx*ny*nz;

//     export_receiver_output = str2bool(catch_parameter("export_receiver_output", file));
//     export_wavefield_output = str2bool(catch_parameter("export_wavefield_output", file));

//     receiver_output_folder = catch_parameter("receiver_output_folder", file); 
//     wavefield_output_folder = catch_parameter("wavefield_output_folder", file);
// }

// void Eikonal::set_acquisition_geometry()
// {
//     std::vector<Geometry *> possibilities = 
//     {
//         new Regular(), 
//         new Circular()
//     };

//     auto type = std::stoi(catch_parameter("geometry_type", file));

//     geometry = possibilities[type];

//     geometry->file = file;

//     geometry->set_geometry();

//     total_shots = geometry->shots.total;
//     total_nodes = geometry->nodes.total;

//     check_geometry_overflow();

//     std::vector<Geometry *>().swap(possibilities); 
// }

// void Eikonal::set_velocity_model()
// {
//     V = new float[nPoints]();

//     import_binary_float(catch_parameter("input_model_file", file), V, nPoints);
// }

// void Eikonal::set_boundaries()
// {
//     nxx = nx + nbxl + nbxr;
//     nyy = ny + nbyl + nbyr;
//     nzz = nz + nbzu + nbzd;

//     volsize = nxx*nyy*nzz;
// }

// void Eikonal::set_slowness_model()
// {
//     S = new float[volsize]();

//     set_velocity_model();

//     expand_boundary(V, S);

//     for (int index = 0; index < volsize; index++) 
//         S[index] = 1.0f / S[index];

//     delete[] V;
// }

// void Eikonal::set_outputs()
// {
//     wavefield_output_samples = nPoints;
//     receiver_output_samples = total_nodes;

//     receiver_output = new float[receiver_output_samples]();
//     wavefield_output = new float[wavefield_output_samples]();
// }

// void Eikonal::get_RAM_usage()
// {
//     struct rusage usage;
//     getrusage(RUSAGE_SELF, &usage);
//     RAM = (int) (usage.ru_maxrss / 1024);
// }

// void Eikonal::get_GPU_usage()
// {
// 	size_t freeMem, totalMem;
// 	cudaMemGetInfo(&freeMem, &totalMem);
//     vRAM = (int) ((totalMem - freeMem) / (1024 * 1024));
//     vRAM -= ivRAM;
// }

// void Eikonal::get_GPU_initMem()
// {
// 	size_t freeMem, totalMem;
// 	cudaMemGetInfo(&freeMem, &totalMem);
//     ivRAM = (int) ((totalMem - freeMem) / (1024 * 1024));
// }

// void Eikonal::check_geometry_overflow()
// {
//     for (int shot = 0; shot < total_shots; shot++)
//     {
//         if ((geometry->shots.x[shot] < 0) || (geometry->shots.x[shot] > (nx-1)*dx) || 
//             (geometry->shots.y[shot] < 0) || (geometry->shots.y[shot] > (ny-1)*dy) ||
//             (geometry->shots.z[shot] < 0) || (geometry->shots.z[shot] > (nz-1)*dz))       
//         throw std::invalid_argument("\033[31mError: shots geometry overflow!\033[0;0m");
//     }

//     for (int node = 0; node < total_nodes; node++)
//     {
//         if ((geometry->nodes.x[node] < 0) || (geometry->nodes.x[node] > (nx-1)*dx) || 
//             (geometry->nodes.y[node] < 0) || (geometry->nodes.y[node] > (ny-1)*dy) ||
//             (geometry->nodes.z[node] < 0) || (geometry->nodes.z[node] > (nz-1)*dz))       
//         throw std::invalid_argument("\033[31mError: nodes geometry overflow!\033[0;0m");
//     }
// }

// void Eikonal::get_travel_times()
// {
//     reduce_boundary(T, wavefield_output);

//     wavefield_output_file = wavefield_output_folder + eikonal_method + "_time_volume_" + std::to_string(nz) + "x" + std::to_string(nx) + "x" + std::to_string(ny) + "_shot_" + std::to_string(shot_id+1) + ".bin";
// }

// void Eikonal::get_first_arrivals()
// {
//     for (int r = 0; r < total_nodes; r++)
//     {
//         float x = geometry->nodes.x[r];
//         float y = geometry->nodes.y[r];
//         float z = geometry->nodes.z[r];

//         float x0 = floorf(x / dx) * dx;
//         float y0 = floorf(y / dy) * dy;
//         float z0 = floorf(z / dz) * dz;

//         float x1 = floorf(x / dx) * dx + dx;
//         float y1 = floorf(y / dy) * dy + dy;
//         float z1 = floorf(z / dz) * dz + dz;

//         int id = ((int)(z / dz)) + ((int)(x / dx))*nz + ((int)(y / dy))*nx*nz;

//         float c000 = wavefield_output[id];
//         float c001 = wavefield_output[id + 1];
//         float c100 = wavefield_output[id + nz]; 
//         float c101 = wavefield_output[id + 1 + nz]; 
//         float c010 = wavefield_output[id + nx*nz]; 
//         float c011 = wavefield_output[id + 1 + nx*nz]; 
//         float c110 = wavefield_output[id + nz + nx*nz]; 
//         float c111 = wavefield_output[id + 1 + nz + nx*nz];

//         float xd = (x - x0) / (x1 - x0);
//         float yd = (y - y0) / (y1 - y0);
//         float zd = (z - z0) / (z1 - z0);

//         float c00 = c000*(1 - xd) + c100*xd;    
//         float c01 = c001*(1 - xd) + c101*xd;    
//         float c10 = c010*(1 - xd) + c110*xd;    
//         float c11 = c011*(1 - xd) + c111*xd;    

//         float c0 = c00*(1 - yd) + c10*yd;
//         float c1 = c01*(1 - yd) + c11*yd;

//         receiver_output[r] = c0*(1 - zd) + c1*zd;
//     }

//     receiver_output_file = receiver_output_folder + eikonal_method + "_data_nRec" + std::to_string(total_nodes) + "_shot_" + std::to_string(shot_id+1) + ".bin";
// }
