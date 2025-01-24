___

## Seismic First-Arrival Toolkit : SeisFAT3D

Modeling, inversion and migration using massive computational parallelism in object-oriented programming.

### Requirements:

- A linux distribution or the Windows Subsystem for linux
- A Nvidia Graphic Processing Unit (GPU)
- Cuda 12 or higher with compatible drivers
- C++ and Python 3 programming languages    
____
### Citation:

```console
@software{paulo_h_b_alves_2024_14170690,
  author       = {Paulo H. B. Alves},
  title        = {phbastosa/SeisFAT3D: Initial},
  month        = November,
  year         = 2024,
  publisher    = {Zenodo},
  version      = {1.0},
  doi          = {10.5281/zenodo.14170690},
  url          = {https://doi.org/10.5281/zenodo.14170690}
}
```
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
```

# Test results


### Modeling test
____

<p float="left">
  <img src="https://github.com/user-attachments/assets/2302a923-5d7b-4b40-a486-e0236bcee9c0" width="32.8%" />
  <img src="https://github.com/user-attachments/assets/1e5d2842-655e-4acb-98ca-79b0fdb78712" width="32.8%" />
  <img src="https://github.com/user-attachments/assets/b0ec784c-d5e6-4061-9057-f6eb514dbf2e" width="32.8%" />
</p>


![modeling_test_results](https://github.com/user-attachments/assets/3dbb3849-aed4-49a6-9617-6247e94a69f0)

![modeling_test_accuracy](https://github.com/user-attachments/assets/b97aa87f-685d-4487-9a8c-a49812a31e36)

____

### Inversion test
____

<p float="left">
  <img src="https://github.com/user-attachments/assets/18e734d1-e359-40bc-a457-77d652f1a765" width="32.8%" />
  <img src="https://github.com/user-attachments/assets/b7d04aaa-89ad-413a-bca4-ed84e778d53c" width="32.8%" />
  <img src="https://github.com/user-attachments/assets/1aac2342-c16d-4ac0-9c4f-aa10b91d5ccc" width="32.8%" />
</p>

<p float="left">
  <img src="https://github.com/user-attachments/assets/6d520107-94fd-44c2-8f3b-4d69b67e84d6" width="49.7%" />
  <img src="https://github.com/user-attachments/assets/2b9ca071-8e1a-4d3d-be76-135f6db621fd" width="49.7%" />
</p>

<p float="left">
  <img src="https://github.com/user-attachments/assets/3232bbd8-851a-4e8d-835a-798d2deec46d" width="49.7%" />
  <img src="https://github.com/user-attachments/assets/779db0c8-7ed4-4a00-97d6-dac0ecb7055b" width="49.7%" />
</p>

![inversion_test_convergence](https://github.com/user-attachments/assets/81b4d2ec-5d59-4160-8220-9fd0984f098a)

____

### Migration test
____

<p float="left">
  <img src="https://github.com/user-attachments/assets/92bdfc32-a539-4900-8b98-5d61edcdd9e9" width="32.8%" />
  <img src="https://github.com/user-attachments/assets/f8385137-37e2-4936-8bc0-cea9aeb751c7" width="32.8%" />
  <img src="https://github.com/user-attachments/assets/16eda223-2e72-4ee3-80a9-6c823a74e8bf" width="32.8%" />
</p>

![migration_test_result](https://github.com/user-attachments/assets/997f18cc-1243-4458-897f-fb485effbe36)
