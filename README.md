# SYNS_edges
Detect image edges and fit LiDAR surfaces

## Required functions

This code requires additional functions from MATLAB File Exchange. Please download these and add them to your MATLAB path before running the code.

* Circular Statistics Toolbox
https://au.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics

* SPHERE_GRID for Matlab
https://people.sc.fsu.edu/~jburkardt/m_src/sphere_grid/sphere_grid.html

* MATLAB to Point Cloud Library
https://au.mathworks.com/matlabcentral/fileexchange/40382-matlab-to-point-cloud-library

## SYNS dataset

The SYNS dataset can be downloaded at https://syns.soton.ac.uk/

Only the `rep1.hdr` Spheron image and `clean_cloud.pcd` LiDAR point cloud are required for each scene. The code assumes these files are located in the numbered folders in the `SYNS` folder (e.g., `SYNS/1/rep1.hdr`).

## To run

Run demo.m to do the entire edge processing pipeline on a single SYNS image.
