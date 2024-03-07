___

## Seismic First Arrival Toolkit - SeisFAT3D

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
SeisFAT3D$ cd run/
SeisFAT3D/run/$ ./program.sh -h

Usage:

    $ ./program.sh -compile             # Create executables 
    $ ./program.sh -modeling            # Perform eikonal solver          
    $ ./program.sh -inversion           # Perform adjoint-state tomography
    $ ./program.sh -migration           # Perform kirchhoff depth migration

Tests:

    $ ./program.sh -test_modeling       # Perform a small modeling experiment          
    $ ./program.sh -test_inversion      # Perform a small inversion experiment
    $ ./program.sh -test_migration      # Perform a small migration experiment  
```
