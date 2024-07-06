___

## Seismic First Arrival Toolkit - SeisFAT3D

Modeling, inversion and migration using massive computational parallelism in object-oriented programming.

### Publications:

**Alves, P. H. B.**; Moreira, R. M.; Cetale, M.; Lopez, J. First arrival tomography application in low illumination time-lapse reservoir monitoring. In: Sociedade Brasileira de Geofísica. *18th International Congress of the Brazilian Geophysical Society*, 2023. [link](https://sbgf.org.br/mysbgf/eventos/expanded_abstracts/18th_CISBGf/a8c88a0055f636e4a163a5e3d16adab7CISBGf_2023_tomography.pdf)

**Alves, P. H. B.**; Capuzzo, F.; Cetale, M.; Santos, L. 3D refraction travel times accuracy study in high contrasted media. Brazilian Journal of Geophysics, Vol. 41, Nº 1, 2023. [link](https://sbgf.org.br/revista/index.php/rbgf/article/view/2295)
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

    $ ./program.sh -compile              
    $ ./program.sh -modeling                      
    $ ./program.sh -inversion           
    $ ./program.sh -migration           

Tests:

    $ ./program.sh -test_modeling                 
    $ ./program.sh -test_inversion      
    $ ./program.sh -test_migration        
    $ ./program.sh -test_wave_equation  

Visualization:

    $ ./program.sh -check_geometry               
```
