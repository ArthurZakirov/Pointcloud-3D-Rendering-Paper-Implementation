# **3D-Rendering**
<img src="img/meshlab_cat.gif" width="500"/>

 This repository contains my implementations of 3D rendering methods. Basic parts of the code such as data loading are taken from [3D AI Lab](https://www.3dunderstanding.org/index.html), however I implemented the algorithms myself.


### **Implicite Function and Rendering**
This repository contains two methods to create an implicit function that represents a given pointcloud:<br>
- [Hoppe'92](https://graphics.pixar.com/library/Reconstruction/paper.pdf)<br>
- [RadialBasisFunction](https://en.wikipedia.org/wiki/Radial_basis_function)<br>

In order to render the implicit function and discretize it, this repository uses the following method:
- [Marching cubes](https://www.researchgate.net/publication/202232897_Marching_Cubes_A_High_Resolution_3D_Surface_Construction_Algorithm)

### **Setup**
In order to run this project, you need to set up the following libraries with [cmake](https://cmake.org/).
- [ceres](http://ceres-solver.org/)
- [eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page)
- [flann](https://github.com/flann-lib/flann)
- [glog](https://github.com/google/glog)

### **Usage**
In order to run the rendering algorithm, make sure to put the raw point cloud in the ```Data``` directory and set the path inside the ```main.cpp``` to ```std::string filenameIn = "../../Data/surface_data/normalized.pcb";```.