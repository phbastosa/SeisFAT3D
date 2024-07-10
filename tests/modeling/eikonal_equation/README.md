## Eikonal modeling benchmark:

### The following test was performed using a laptop with the configuration below:

* CPU: AMD Ryzen 5 2500U with Radeon Vega Mobile Gfx (8) @ 2.000GHz
* GPU: NVIDIA GeForce GTX 1050 Mobile 4GB 
* RAM: 16 GB       
___

First of all, you have to compile the code to generate the executables. Make sure you're inside the run folder.

```console
SeisFAT3D/run$ ./program -compile
```

After that, you just need to perform the test.

```console
SeisFAT3D/run$ ./program -test_modeling
```
The results will appear in .png images, as follows

### Objective: Verify accuracy, performance and memory usage of the eikonal modeling methods 

![modelTest](https://github.com/phbastosa/SeisFAT3D/assets/44127778/27e385e3-909e-47ad-b9b7-6f959968381b)

![accuracyTest](https://github.com/phbastosa/SeisFAT3D/assets/44127778/88ead2cd-8b9f-4d53-9260-c5d5ece89017)

![benchmarkTest](https://github.com/phbastosa/SeisFAT3D/assets/44127778/6ebe7fd0-51cb-481a-84e4-88ddeda92811)

