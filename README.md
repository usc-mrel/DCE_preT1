# Matlab codes for Sparse Pre-Contrast T1 Mapping for High-Resolution Whole-Brain DCE-MRI
[Magnetic Resonance Engineering Laboratory (MREL)](https://mrel.usc.edu/)  
**University of Southern California**
## Code Structure
### Demo script
**demo.m**              Main script demonstrating sparse T1 mapping using a 256x256x12 anatomically realistic brain tumor Digital Reference Object (**DRO**).
### Functions
**spgr.m**              VFA images signal intensity computation.  
**genKspace.m**         VFA k-space signal computation.  
**applyNoise.m**        Generate and apply syntheszied noise.  
**applyU.m**            Generate and apply sparse sampling pattern.  
**P_SEN**               Reconstruction entry.  
**P2sig**               Cost function value and gradient computation.  
**argmin**              Solver function.  
### Folders
**DRO**                 DRO data generation.  
**Utils**               Utility functions, e.g. Fourier transform.  
**results**             For results storage.