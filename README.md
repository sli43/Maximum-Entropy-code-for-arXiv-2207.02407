# Maximum-Entropy-code-for-arXiv-2207.02407

Copyright [2023] [Shaozhi Li]

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.


The Maximum Entropy code is used to calculate the Eliashberg function with the self-energy input. 
Contact lishaozhiphys@gmail.com
*******************************************************

The kernel functions for the real part and imagniary part of the self-energy are different.
The input samples are included in each diretory. One can use these testing data to run the code.

1. System requirements

   a. Linux system
   
   b. intel compiler with mkl library
   
   c. testing version: intel 2018.1.163
   
2. Installation guide

   use makefile to compile the code: make
   
3. Demo

   input file: 
   ***_band_input.dat
   *** denotes "alpha" and "beta" in the real part and imaginary part directories. The input file includes the measured self-energy for *** band
               
   output file: 
   Eliashberg_XXXX.dat. XXXX denotes the index number (1 or 2). This output includes the eliashberg function
   
   back_selfenergy_XXXX.dat. The selfenergy generated using the Eliashberg_XXXX.dat
   
   XXalpha_Prob.dat. The probability of each alpha value in the Maximum Entropy calculations.
  
   Using testing dataset:
  
   a. To use alpha band data, set the runid as 1 in main.f90 on line 6.
   
   b. To use beta band data, set the runid as 2 in main.f90 on line 6.
                 
   run time: 0.5 hour
                 
4. instructions for use

   go to each directory. Type "make" to compile the code. In the end, you will get an application named as "main". use "./main" to run the code
  
