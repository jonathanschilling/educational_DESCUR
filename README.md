## Educational DESCUR

This is a stripped-down version of the implementation of DESCUR.
It is forked from the `master` branch of the corresponding [Stellarator-Tools](https://github.com/ORNL-Fusion/DESCUR) repository.

The `cmake` build system used here is borrowed
from [hiddenSymmetries/VMEC2000](https://github.com/hiddenSymmetries/VMEC2000)
and from [ORNL-Fusion/LIBSTELL](https://github.com/ORNL-Fusion/LIBSTELL).

## Building
This is a fairly standard CMake setup, if you are used to it.
Here is how it works:
 * Create a directory `build` in the main folder: `mkdir build`
 * Go into the `build` directory: `cd build`
 * Run CMake: `cmake ..`
 * Execute the actual build process: `make` (optional multi-threaded build: `make -j`)
 * The DESCUR executable `xdescur` is then located in `build/bin` with respect to the main folder.

## Example Execution
 * Create a directory `test` in the main folder: `mkdir test`
 * Change into the `test` dir: `cd test`
 * Run DESCUR: `../build/bin/xdescur`
 * Provide the following interactive information:
    ```
     Enter spectral convergence parameter
     =0 for polar,  =1 for equal arclength
     >1 for smaller spectral width
     =4 corresponds to STELLOPT choice): 4
     Use (default) VMEC-compatible compression (V)
     or Advanced Hirshman-Breslau compression (A): v
     Select source of curve data to be fit:
     0 :  Points from file
     2 :  Solove'ev Equilibrium
     3 :  Assorted 2-D shapes

    3
      Shape to fit: 1=ellipse; 2=D; 3=bean; 4=belt; 5=square; 6=D3D-asym; 7=heliac)
    2
    ```
 * `DESCUR` will then create a spectrally condensed representation of a D-shaped boundary:
    ```
     ORDERING SURFACE POINTS
     Average elongation = 1.4218E+00
     Raxis =   3.0000E+00 Zaxis =  -6.5321E-18
     Number of Theta Points Matched =  2500
     Number of Phi Planes =     1
     Max Poloidal Mode Number =    30
     Max Toroidal Mode Number =     1

                      Fitting toroidal plane #   1

     RAXIS =  3.000E+00 ZAXIS = -6.532E-18 R10 =  1.723E+00

     ITERATIONS    RMS ERROR    FORCE GRADIENT    <M>    MAX m   DELT
           1       4.837E-03       8.903E-04      1.17      30  1.00E+00
         100       3.216E-05       2.847E-05      1.16      30  9.55E-01
         200       1.510E-08       1.084E-08      1.16      30  9.55E-01
         206       1.034E-08       7.956E-09      1.16      30  9.55E-01

     ANGLE CONSTRAINTS WERE APPLIED FOR POLAR DAMPING EXPONENT (PEXP) =  4
     RM**2 + ZM**2 SPECTRUM COMPUTED WITH P =     4.00 AND Q =     1.00
     TIME:   7.98E-01 SEC.

      OUTPUTTING FOURIER COEFFICIENTS TO OUTCURVE FILE
      ```
 * One last question remains to be answered:
    ``` Do you wish to plot this data (Y/N)? n```
   (unless you have [`DESCUR_PLOT`](https://github.com/PrincetonUniversity/STELLOPT/tree/v251/DESCUR_PLOT) compiled and installed)
 * The spectrally condensed output is then in the file `outcurve`.
   Some more output can be found in the `plotout` file.

For more information, please consider also the corresponding page in the `STELLOPT` Wiki:
https://princetonuniversity.github.io/STELLOPT/DESCUR


