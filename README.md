[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-22041afd0340ce965d47ae6ef1cefeee28c7c493a6346c4f15d667ab976d596c.svg)](https://classroom.github.com/a/DDm8twLC)
# Assignment 5 Description (11% of total grade) #

The task for this assignment is the implementation of Vector field visualization.

## Reading assignments ##

* Data Visualization book, Chapter 6.5, 6.6
* Bruno Jobard and Wilfrid Lefer:
  Creating Evenly-Spaced Streamlines of Arbitrary Density,
  Eurographics Workshop on Visualization in Scientific Computing 1997
  http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.29.9498

## Basic Tasks ##

* Vector field visualization in 2D using glyphs
* Streamline integration using Euler and Runge-Kutta methods
* Pathline integration using Euler and Runge-Kutta methods
* Color-coded scalar field visualization using colormaps

## Minimum Requirements ##

+ Glyph visualization: draw glyphs ( arrows or something that indicate direction ) (15 points)
    * Adjustable length: choose between arrows with constant length or arrows with the length depending on the speed (vector magnitude) (5 points)
    * Change (sub)sampling rate for arrows or glyphs (vector grid) (5 points)

+ Streamlines using a) Euler (15 points) and b) RK2 (10 points)
    * Every time you seed a point keep the streamline of that point (don't reset streamlines) (-3 points penalty if the streamlines are removed)
    * When you change a time-slice the streamline should update correctly (recalculate) (5 points)
    * Do bilinear interpolation of vectors (don't snap to nearest vector) (5 points)
    * Stopping condition depends on the accumulated length of the streamline or if the queried vector is almost zero or if you hit a boundary (-3 points penalty if those conditions are not implemented)
    * Do backward integration with forward integration (-3 points penalty if only one integration is implemented)

+ Pathlines using a) Euler (20 points) and b) RK2 (10 points)
    * Every time you seed a point keep the pathline of that point (don't reset pathlines) (-3 points penalty if the streamlines are removed)
    * Do trilinear interpolation of vectors (don't snap to nearest vector and timeslice) (5 points)
    * Stopping condition depends on the accumulated length of the streamline or if the queried vector is almost zero or if you hit a boundary (-3 points penalty if those conditions are not implemented)
    * Do backward integration with forward integration (-3 points penalty if only one integration is implemented)

+ Color-coded scalar field visualization (15 points)
    * Implement at least two colormaps (e.g. rainbow, cool-warm) applied to the background scalar field in the fragment shader (10 points)
    * User adjustable blend factor between grayscale and color-mapped scalar field (5 points)

* Adjustable dt value by user (5 points)



## Bonus ##
* Release multiple streamline seeds in horizontal or vertical rake (+5 points)
* Switch between different background scalar field images (+3 points)
* RK4 for both streamlines (+5 points) and pathlines (+5 points)

## Notes ##
* The description and the download location of the data can be found in the file README_data.txt
* The vector data is in the vector_array, where every 3 elements in the array will give you the xyz component of a vector. After one 2D slice is finished, and if there is more than one time step, the vector data continues for the next time step.
  So the size of the array = 3 * width * height * number_of_timesteps. More information can be found in the file README_data.txt found in the source code.
* There aren't prototypes for every function you might need. Create functions as you need them.


## Screenshots for Minimum Requirements Solution ##
Glyphs Visualization

![glyphs.png](https://bitbucket.org/repo/Mq6ygx/images/425286000-glyphs.png)

Streamlines

![streamlines.png](https://bitbucket.org/repo/Mq6ygx/images/154883112-streamlines.png)

Pathlines

![pathlines.png](https://bitbucket.org/repo/Mq6ygx/images/3862458026-pathlines.png)
