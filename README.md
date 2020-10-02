# Matlab Codes for Sparse Pre-Contrast T1 Mapping for High-Resolution Whole-Brain DCE-MRI
[Magnetic Resonance Engineering Laboratory (MREL)](https://mrel.usc.edu/)  
**University of Southern California**
## Code Structure
### Demo script
**demo.m**              Main script demonstrating sparse T1 mapping using a 256x256x12 anatomically realistic brain tumor Digital Reference Object (**DRO**).
### Functions
**spgr.m**              VFA images signal intensity computation.  
**genKspace.m**         VFA k-space signal computation.  
**addNoise.m**          Generate and apply synthesized noise.  
**applyU.m**            Generate and apply sparse sampling pattern.  
**P_SEN.m**             Reconstruction function.  
**P2sig.m**             Compute cost function value and gradient.  
**argmin.m**            Solver function.  
### Folders
**DRO**                 DRO and noise data.  
**Utils**               Utility functions, e.g. Fourier transform.  
**Results**             For results storage.  
**Images**              Sample images, including: ground truth M0 and T1 maps, reconstructed M0 and T1 maps, fractional difference maps.
### Relevant References
1. Bosca RJ, Jackson EF. Creating an anthropomorphic digital MR phantom - An extensible tool for comparing and evaluating quantitative imaging algorithms. Phys Med Biol. 2016;61(2):974-982. doi:10.1088/0031-9155/61/2/974
2. Bliesener Y, Lingala SG, Haldar JP, Nayak KS. Impact of (k,t) sampling on DCE MRI tracer kinetic parameter estimation in digital reference objects. Magn Reson Med. 2020;83(5):1625-1639. doi:10.1002/mrm.28024
3. Stanisz GJ, Odrobina EE, Pun J, et al. T1, T2 relaxation and magnetization transfer in tissue at 3T. Magn Reson Med. 2005;54(3):507-512. doi:10.1002/mrm.20605
4. Müller A, Jurcoane A, Kebir S, et al. Quantitative T1-mapping detects cloudy-enhancing tumor compartments predicting outcome of patients with glioblastoma. Cancer Med. 2017;6(1):89-99. doi:10.1002/cam4.966
5. Hattingen E, Müller A, Jurcoane A, et al. Value of quantitative magnetic resonance imaging T1-relaxometry in predicting contrast-enhancement in glioblastoma patients. Oncotarget. 2017;8(32):53542-53551. doi:10.18632/oncotarget.18612
6. Badve C, Yu A, Dastmalchian S, et al. MR fingerprinting of adult brain tumors: Initial experience. In: American Journal of Neuroradiology. Vol 38. American Society of Neuroradiology; 2017:492-499. doi:10.3174/ajnr.A5035

Got questions? Contact zhibozhu@usc.edu.
