This code accompanies the paper "Friction Modeling of Continuum Robots Through
Linear Complementarity Problem" (under review) by Jia Shen, Junhyoung Ha, and Yue Chen.  
This code implements the friction model for continuum robots by solving the Linear Complementarity Problem (LCP). It includes simulations that reproduce the case studies presented in the associated paper. 
The simulations are located in the "simulations" folder:
    - simu_ctr.m: A concentric tube robot (CTR) interacting with a plane.
    - simu_rod_plane.m: An elastic rod interacting with a plane.
    - simu_rod_pipe.m: An elastic rod inside a straight pipe.
    - simu_rod_ring_plane.m: An elastic rod interacting with both a plane and a ring.

Requirements:  
- Programming language: MATLAB, tested on Version R2023a.

Installation:  
1. Download or clone this repository.
2. Add all folders in the repository to the MATLAB path:
    In MATLAB, use the addpath function or go to Home > Set Path > Add with Subfolders and select the root folder of this project.

