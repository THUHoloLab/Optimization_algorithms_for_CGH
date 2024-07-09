<p align="center">
<img src="https://github.com/Sui00004/Optimization-algorithms-for-computer-generated-holography/blob/main/1_CGH.jpg", height="300">
</p>

## CGH Optimization
> Non-convex optimization algorithms to synthesize computer-generated holograms

Computer-generated holography (CGH) involves computationally generating a hologram and optically reconstructing an object's wavefront, which provides an approach to digitally modulate a volumetric wavefront. 

Computer-generated holograms can be encoded on various types of holographic media, including diffractive optical elements, metasurfaces, and spatial light modulators. Even so, the algorithms for hologram synthesis can be universally applied. Existing optimization frameworks applied to CGH can be concluded with the following categories: alternative projections, first-order gradient descent, and second-order gradient descent.

<p align="center">
<img src="https://github.com/Sui00004/Optimization-algorithms-for-computer-generated-holography/blob/main/2_Frameworks.jpg", height="230">
</p>

Codes to achieve hologram synthesis with these optimization frameworks is provided here. The algorithms to realize functional hologram synthesis using different frameworks are written in MATLAB. If you use any of these codes, please cite it as follows:

> Sui X, He Z, Chu D, & Cao L. Non-convex Optimization for Inverse Problem Solving in Computer-generated Holography. Light: Science & Applications, 13(158), 1-23 (2024) [[link](https://www.nature.com/articles/s41377-024-01446-w)] [[bibtex](Inverse_CGH.bib)]

### Alternating projections
Alternative projections can be achieved by a pair of elementary projections repeatedly occurring in the optimization, which construct an iterative computation loop. Specially for CGH, alternating projections are applied to two enclosed sets associated with potential object solutions and potential hologram solutions. 

1. Gerchberg-Saxton (GS) algorithm: [[2D optimization (phase-only hologram)](https://github.com/THUHoloLab/Optimization_algorithms_for_CGH/blob/main/Optimization_algorithms_for_CGH/Alternative%20projection/Main1_GS_2D_FFT_POH.m)]
2. Iterative Fourier-transform algorithm (IFTA):
   * [[2D optimization (phase-only hologram)](https://github.com/THUHoloLab/Optimization_algorithms_for_CGH/blob/main/Optimization_algorithms_for_CGH/Alternative%20projection/Main2_IFTA_2D_FFT_POH.m)]
   * [[2D optimization with signal windows (phase-only hologram)](https://github.com/THUHoloLab/Optimization_algorithms_for_CGH/blob/main/Optimization_algorithms_for_CGH/Alternative%20projection/Main3_IFTA_2D_signalwindow_FFT_2D_POH.m)]
   * [[2D optimization with soft encoding (phase-only hologram)](https://github.com/THUHoloLab/Optimization_algorithms_for_CGH/blob/main/Optimization_algorithms_for_CGH/Alternative%20projection/Main4_IFTA_2D_soft_encoding_FFT_2D_POH.m)]
3. Iterative algorithm with angular spectrum theory:
   * [[2D optimization (complex hologram)](https://github.com/THUHoloLab/Optimization_algorithms_for_CGH/blob/main/Optimization_algorithms_for_CGH/Alternative%20projection/Main5_IFTA_2D_Angularspectrum_CH.m)]
   * [[2D optimization (phase-only hologram)](https://github.com/THUHoloLab/Optimization_algorithms_for_CGH/blob/main/Optimization_algorithms_for_CGH/Alternative%20projection/Main6_IFTA_2D_Angularspectrum_POH.m)].

### First-order gradient descent
The inverse problem of hologram synthesis in CGH can also be cast as the optimization of a parameterized objective function requiring minimization with respect to its parameters. Since the choice of the objective function is often stochastic and differentiable with respect to its parameters, stochastic gradient descent (SGD) is considered as an efficient and effective first-order gradient descent framework for optimization. 

1. Stochastic gradient descent (SGD) with single Fourier-transform propagation:
   * [[2D optimization (complex hologram)](https://github.com/THUHoloLab/Optimization_algorithms_for_CGH/blob/main/Optimization_algorithms_for_CGH/SGD/Main2_SGD_2D_FFT_CH.m)]
   * [[2D optimization (phase-only hologram)](https://github.com/THUHoloLab/Optimization_algorithms_for_CGH/blob/main/Optimization_algorithms_for_CGH/SGD/Main1_SGD_2D_FFT_POH.m)] 
3. Stochastic gradient descent (SGD) with angular spectrum theory:
   * [[2D optimization (complex hologram)](https://github.com/THUHoloLab/Optimization_algorithms_for_CGH/blob/main/Optimization_algorithms_for_CGH/SGD/Main4_SGD_2D_Angularspectrum_CH.m)]
   * [[2D optimization (phase-only hologram)](https://github.com/THUHoloLab/Optimization_algorithms_for_CGH/blob/main/Optimization_algorithms_for_CGH/SGD/Main3_SGD_2D_Angularspectrum_POH.m)]

### Second-order gradient descent
The second-order gradient descent is implemented by the quasi-Newton method here, which minimizes the loss function by constructing and storing a series of matrices that approximate the Hessian or inverse Hessian matrix of the loss function. 

1. The quasi-Newton method with single Fourier-transform propagation:
   * [[2D optimization (complex hologram)](https://github.com/THUHoloLab/Optimization_algorithms_for_CGH/blob/main/Optimization_algorithms_for_CGH/Quasi-Newton/Main2_quasiNewton_2D_FFT_CH.m)]
   * [[2D optimization (phase-only hologram)](https://github.com/THUHoloLab/Optimization_algorithms_for_CGH/blob/main/Optimization_algorithms_for_CGH/Quasi-Newton/Main1_quasiNewton_2D_FFT_POH.m)] 
2. Stochastic gradient descent (SGD) with angular spectrum theory:
   * [[2D optimization (complex hologram)](https://github.com/THUHoloLab/Optimization_algorithms_for_CGH/blob/main/Optimization_algorithms_for_CGH/Quasi-Newton/Main4_quasiNewton_2D_Angularspectrum_CH.m)]
   * [[2D optimization (phase-only hologram)](https://github.com/THUHoloLab/Optimization_algorithms_for_CGH/blob/main/Optimization_algorithms_for_CGH/Quasi-Newton/Main3_quasiNewton_2D_Angularspectrum_POH.m)]

In our demonstration, the 512×512 holograms are equally computed on a PC with an Intel Core i9-9900K 3.6 GHz CPU and 32.0 GB of RAM, and an NVIDIA GeForce RTX 2080Ti GPU. Two diffractive propagation methods including the FFT and ASM are compared, where the propagation distance of ASM is 50 mm. 

<p align="center">
<img src="https://github.com/Sui00004/Optimization-algorithms-for-computer-generated-holography/blob/main/2D_optimization.jpg", height="300">
</p>

### Optimization for 3D holograms
The flexibility of the optimization frameworks also brings with it a diversity of pipelines for the hologram synthesis of 3D objects. Many diffractive propagation models are applicable for volumetric optimization. Here, for a better illustration of different optimization pipelines, we unify them into the band-limited ASM.

1. The superposition method:
   * Alternating projections:
     - [[complex hologram](https://github.com/THUHoloLab/Optimization_algorithms_for_CGH/blob/main/Optimization_algorithms_for_CGH/Alternative%20projection/Main9_Superposition3D_Angularspectrum_CH.m)]
     - [[phase-only hologram](https://github.com/THUHoloLab/Optimization_algorithms_for_CGH/blob/main/Optimization_algorithms_for_CGH/Alternative%20projection/Main10_Superposition3D_Angularspectrum_POH.m)]
   * Stochastic gradient descent:
     - [[complex hologram](https://github.com/THUHoloLab/Optimization_algorithms_for_CGH/blob/main/Optimization_algorithms_for_CGH/SGD/Main8_SGD_Superposition3D_Angularspectrum_CH.m)]
     - [[phase-only hologram](https://github.com/THUHoloLab/Optimization_algorithms_for_CGH/blob/main/Optimization_algorithms_for_CGH/SGD/Main7_SGD_Superposition3D_Angularspectrum_POH.m)]
   * Quasi-Newton method:
     - [[complex hologram](https://github.com/THUHoloLab/Optimization_algorithms_for_CGH/blob/main/Optimization_algorithms_for_CGH/Quasi-Newton/Main6_quasiNewton_superposition3D_Angularspectrum_CH.m)]
     - [[phase-only hologram](https://github.com/THUHoloLab/Optimization_algorithms_for_CGH/blob/main/Optimization_algorithms_for_CGH/Quasi-Newton/Main5_quasiNewton_superposition3D_Angularspectrum_POH.m)]
2. The global method:
   * Alternating projections:
     - [[complex hologram](https://github.com/THUHoloLab/Optimization_algorithms_for_CGH/blob/main/Optimization_algorithms_for_CGH/Alternative%20projection/Main11_Global3D_Angularspectrum_CH.m)]
     - [[phase-only hologram](https://github.com/THUHoloLab/Optimization_algorithms_for_CGH/blob/main/Optimization_algorithms_for_CGH/Alternative%20projection/Main12_Global3D_Angularspectrum_POH.m)]
   * Stochastic gradient descent:
     - [[complex hologram](https://github.com/THUHoloLab/Optimization_algorithms_for_CGH/blob/main/Optimization_algorithms_for_CGH/SGD/Main6_SGD_Global3D_angularspectrum_CH.m)]
     - [[phase-only hologram](https://github.com/THUHoloLab/Optimization_algorithms_for_CGH/blob/main/Optimization_algorithms_for_CGH/SGD/Main5_SGD_Global3D_angularspectrum_POH.m)]
   * Quasi-Newton method:
     - [[complex hologram](https://github.com/THUHoloLab/Optimization_algorithms_for_CGH/blob/main/Optimization_algorithms_for_CGH/Quasi-Newton/Main8_quasiNewton_global3D_Angularspectrum_CH.m)]
     - [[phase-only hologram](https://github.com/THUHoloLab/Optimization_algorithms_for_CGH/blob/main/Optimization_algorithms_for_CGH/Quasi-Newton/Main7_quasiNewton_Global3D_Angularspectrum_POH.m)]  
3. The sequential method:
   * Alternating projections:
     - [[complex hologram](https://github.com/THUHoloLab/Optimization_algorithms_for_CGH/blob/main/Optimization_algorithms_for_CGH/Alternative%20projection/Main7_Sequential3D_Angularspectrum_CH.m)]
     - [[phase-only hologram](https://github.com/THUHoloLab/Optimization_algorithms_for_CGH/blob/main/Optimization_algorithms_for_CGH/Alternative%20projection/Main8_Sequential3D_Angularspectrum_POH.m)]



<p align="center">
<img src="https://github.com/Sui00004/Optimization-algorithms-for-computer-generated-holography/blob/main/3D_optimization.jpg", height="400">
</p>

### Questions and feedback

If you have questions, bug reports, or feature requests, please use the [Issues](https://github.com/Sui00004/Optimization-algorithms-for-computer-generated-holography/issues) section to report them.

