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
```
___

### Modeling Benchmark:

Basic modeling scheme with all steps to perform it.

___

### Inversion Benchmark:

Basic inversion scheme with all steps to perform it.

___

### Migration Benchmark:

Basic migration scheme with all steps to perform it.

___

