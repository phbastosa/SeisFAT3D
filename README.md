___

## Seismic First Arrival Toolkit

Modeling, inversion and migration using massive computational parallelism in object-oriented programming.

### Publications:

**Alves, P. H. B.**; Moreira, R. M.; Cetale, M.; Lopez, J. First arrival tomography application in low illumination time-lapse reservoir monitoring. In: Sociedade Brasileira de Geofísica. *18th
International Congress of the Brazilian Geophysical Society*, 2023. [link](https://sbgf.org.br/mysbgf/eventos/expanded_abstracts/18th_CISBGf/a8c88a0055f636e4a163a5e3d16adab7CISBGf_2023_tomography.pdf)

**Alves, P. H. B.**; Capuzzo, F.; Cetale, M.; Santos, L. 3D refraction travel times accuracy study in high contrasted media. In: Sociedade Brasileira de Geofísica. *IX Simpósio Brasileiro de Geofı́sica*, 2022. [link](https://sbgf.org.br/mysbgf/eventos/expanded_abstracts/IX_SimBGf/session/M%C3%A9todos%20Geof%C3%ADsicos%20e%20Geof%C3%ADsica%20Computacional/3D%20refraction%20travel%20times%20accuracy%20study%20in%20high%20contrasted%20media.pdf)

___

### Requirements:

- A linux distribution or the Windows Subsystem for linux
- A Nvidia Graphic Processing Unit (GPU)
- Cuda 12 or higher with compatible drivers
- C++ and Python 3 programming languages    
____

### Usage:

```console
seismic_imaging_3D$ cd run/
seismic_imaging_3D/run/$ ./program.sh -h

Usage:

    $ ./program.sh -compile        # Create executables 
    $ ./program.sh -modeling       # Perform eikonal solver          
    $ ./program.sh -inversion      # Perform adjoint-state tomography
    $ ./program.sh -migration      # Perform kirchhoff depth migration

Tests:

    $ ./program.sh -test_modeling       # Perform a small modeling experiment          
    $ ./program.sh -test_inversion      # Perform a small inversion experiment
    $ ./program.sh -test_migration      # Perform a small migration experiment         
```

The following tests was performed using a laptop with the configuration below:

* CPU: AMD Ryzen 5 2500U with Radeon Vega Mobile Gfx (8) @ 2.000GHz
* GPU: NVIDIA GeForce GTX 1050 Mobile 4GB 
* RAM: 16 GB       
___

### Modeling Benchmark:

#### Objective: Verify accuracy, performance and memory usage of the eikonal modeling methods 

First of all, you have to compile the code to generate the executables. Make sure you're inside the run folder.

```console
$ ./program -compile
```

After that, you just need to perform the test.

```console
$ ./program -test_modeling
```
The results will appear in .png images, as follows

![modelTest](https://github.com/phbastosa/first_break_imaging_3D/assets/44127778/d68bdfee-36de-4502-a343-d14106599539)

![accuracyTest](https://github.com/phbastosa/first_break_imaging_3D/assets/44127778/88ead2cd-8b9f-4d53-9260-c5d5ece89017)

![benchmarkTest](https://github.com/phbastosa/first_break_imaging_3D/assets/44127778/6ebe7fd0-51cb-481a-84e4-88ddeda92811)

___

### Inversion Benchmark:

As the same way in modeling test you have to compile the program and then

```console
$ ./program -test_inversion
```

#### The reference model: 

![trueModel](https://github.com/phbastosa/first_break_imaging_3D/assets/44127778/54699f42-acb2-41a6-a4c7-46e6f1d2559a) 

#### The initial model:

![initModel](https://github.com/phbastosa/first_break_imaging_3D/assets/44127778/f864022d-7a8c-4515-a22c-db8ab26ce58e)

#### The convergence curve for both, the classical and the adjoint state first arrival tomography

![curve](https://github.com/phbastosa/first_break_imaging_3D/assets/44127778/034a1006-f518-43ed-8fe8-1ca212adca59)

#### Some statistics of the procedures

|  Type                    |  Elapsed time  | RAM usage  | GPU memory usage | 
| ------------------------ | -------------- | ---------- | ---------------- |
| Modeling                 |      4.3 s     |   129 MB   |       4 MB       | 
| Least Squares tomograhy  |     17.3 s     |   142 MB   |       4 MB       | 
| Adjoint State tomography |     51.9 s     |   145 MB   |      10 MB       |  
-----------------------------------------------------------------------------

#### The final model generated by the classical tomography:

![final_model_ls](https://github.com/phbastosa/first_break_imaging_3D/assets/44127778/4b35929b-cd63-4b09-a56c-8fb9609d003e) 

#### The final model generated by the adjoint state tomography

![final_model_adj](https://github.com/phbastosa/first_break_imaging_3D/assets/44127778/26fd8d32-dbc9-4c17-95b0-3dea5a53847e)

___

### Migration Benchmark:

There is no implementation.

___

