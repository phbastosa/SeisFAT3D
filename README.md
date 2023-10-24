___

## Seismic First Arrival Toolkit

Modeling, inversion and migration using massive computational parallelism in object-oriented programming.
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
seismic_imaging_3D/run/$ ./program.sh -help

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
___

### Modeling Benchmark:

#### Objective: Verify accuracy, performance and memory usage of the eikonal modeling methods 

First of all, yiou have to compile the code to generate the executables. Make sure you're in run folder.

```console
$ ./program -compile
```

After that, you just need to perform the test.

```console
$ ./program -test_modeling
```
The results will appear in .png images, as follows

![modelTest](https://github.com/phbastosa/first_break_imaging_3D/assets/44127778/1f3a5a3b-a9d3-431c-9f8d-ddc3e3752dd1)

![accuracyTest](https://github.com/phbastosa/first_break_imaging_3D/assets/44127778/bcfbdf24-6236-4867-87f8-b53a840ba8dd)

![benchmarkTest](https://github.com/phbastosa/first_break_imaging_3D/assets/44127778/f2b0dfc8-3f7a-497e-acac-aae0b854fd38)
___

### Inversion Benchmark:

Basic inversion scheme with all steps to perform it.

___

### Migration Benchmark:

Basic migration scheme with all steps to perform it.

___

