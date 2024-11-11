___

## Seismic First Arrival Toolkit - SeisFAT3D

Modeling, inversion and migration using massive computational parallelism in object-oriented programming.

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
```

# Test results

### Modeling test

<p float="left">
  <img src="https://github.com/user-attachments/assets/2302a923-5d7b-4b40-a486-e0236bcee9c0" width="32.8%" />
  <img src="https://github.com/user-attachments/assets/1e5d2842-655e-4acb-98ca-79b0fdb78712" width="32.8%" />
  <img src="https://github.com/user-attachments/assets/b0ec784c-d5e6-4061-9057-f6eb514dbf2e" width="32.8%" />
</p>


![modeling_test_results](https://github.com/user-attachments/assets/8913e790-c2fe-4e95-aab1-592e54e9054d)

### Inversion test

<p float="left">
  <img src="https://github.com/user-attachments/assets/18e734d1-e359-40bc-a457-77d652f1a765" width="32.8%" />
  <img src="https://github.com/user-attachments/assets/b7d04aaa-89ad-413a-bca4-ed84e778d53c" width="32.8%" />
  <img src="https://github.com/user-attachments/assets/1aac2342-c16d-4ac0-9c4f-aa10b91d5ccc" width="32.8%" />
</p>

<p float="left">
  <img src="https://github.com/user-attachments/assets/6d520107-94fd-44c2-8f3b-4d69b67e84d6" width="49.7%" />
  <img src="https://github.com/user-attachments/assets/2b9ca071-8e1a-4d3d-be76-135f6db621fd" width="49.7%" />
</p>

![inversion_test_convergence](https://github.com/user-attachments/assets/949402ff-6645-47b6-ae00-f7bfe0f758be)

### Migration test


<p float="left">
  <img src="https://github.com/user-attachments/assets/1f1960ac-2611-427a-8972-71971819c4aa" width="32.8%" />
  <img src="https://github.com/user-attachments/assets/36a6b804-34eb-4fa7-ba61-9d06745b7d5c" width="32.8%" />
  <img src="https://github.com/user-attachments/assets/99da2e61-403c-44b5-9df7-ed39ea1ca565" width="32.8%" />
</p>

![migration_test_result](https://github.com/user-attachments/assets/1274cf57-5bd6-49a9-b224-0fdeee14c74f)