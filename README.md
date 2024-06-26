# Impedance Examples

Scripts that generate data, solve inverse problems, and make plots
which appear in the paper "Reconstructing the shape and material
parameters of dissipative obstacles using an
impedance model" by Travis Askham and Carlos Borges 

## dependencies/reproducibility

- The examples were run using MATLAB R2023b
- Requires the Optimization Toolbox and Image Processing Toolbox
- Uses a copy of the inverse-obstacle-scattering-2d MATLAB package, the
specific version used is in the inverse-obstacle-scattering2d folder.
The latest version is available at
https://gitub.com/fastalgorithms/inverse-obstacle-scattering2d but the
scripts here may not be compatible with the latest version.
This package must be installed. See the instructions in the
inverse-obstacle-scattering2d folder. 
- Note: The files in the inverse-obstacle-scattering-2d folder were obtained
from the devel-impedance branch at commit
https://github.com/flatironinstitute/inverse-obstacle-scattering2d/commit/861a3de5fff11ba0b31ccd4c819b731ee6137065 
- Uses v1.0.0 of the chunkIE MATLAB package, which is included
in the chunkie folder.
The latest version is available at
https://gitub.com/fastalgorithms/chunkie but the
scripts here may not be compatible with the latest version.

## generating forward data

There are two folders imp-data and trans-data which contain scripts
for generating the data and solving inverse scattering problems with
impedance and transmission data, respectively.

The forward data can be generated in the imp-data/data-gen and
trans-data/data-gen folders. In those folders, there are functions
called "run_impedance_generator" and "run_transmission_generator".
These contain predefined test ID numbers that will generate
specific data sets. This can be done by calling

```
run_transmission_generator(test_id)
```
in MATLAB in the imp-data/data-gen folder. This will generate a
file in imp-data/data-out that can be used as data for an inverse
solve. Generating transmission data is similar.

The following test_id numbers correspond to tests that are
in the paper.

Data for first revision of paper:
- Section 4.1 (impedance): test_ids [16,17,18,19]
- Section 4.2.1 (transmission): test_ids [59,60,61,62,63]
- Section 4.2.2 (transmission): same as 4.2.1
- Section 4.2.3 (transmission): test_ids [66,67,68,69,70]
- Section 4.2.4 (transmission): same as 4.2.1 plus test_ids [64,65]
- Section 4.3 (transmission with corners using chunkIE): test_ids [71,72,73,74,75]
this was the only test that uses chunkIE to generate the data

## solving inverse problems

The inverse problems can be run in the imp-data/tests and trans-data/tests
folders. These folders contain a function called inversetest_runner
which will run the inverse solver with various default settings.

Tests for first revision of paper:
- Section 4.1 (impedance): for each data set, do
inversetest_runner(test_id) [ABV model],
inversetest_runner(test_id,true,[],[],[],true) [CH model with constraints],
and inversetest_runner(test_id,[],[],true) [Fourier model].
Also run inversetest_runner(16,[],[],[],[],[],[],[],[],true) to
run the constant model initializer example.
- Section 4.2.1 (transmission): for each data set, do
inversetest_runner(test_id)
- Section 4.2.2 (transmission): for each data set, do
inversetest_runner(test_id,[],[],[],[],[],[],[],1e-1)
- Section 4.2.3 (transmission): for each data set, do
inversetest_runner(test_id)
- Section 4.2.4 (transmission): for each data set, do
inversetest_neumann_runner(test_id) [does the neumann model] and
inversettest_runner(test_id,[],[],true,[],[],[],[],[],[],0)
[does the constant impedance model]
- Section 4.3 (transmission with a corner using chunkIE): for each
data set, do inversetest_runner(test_id)
and inversetest_runner(test_id,[],[],[],[],[],[],5e-2) which runs the
inverse impedance model with the same choice of C_H=0.9 as the other
runs and a different choice of C_H=0.95 (this is set by eps_curv=1-C_H)
For these problems (with corners), more curvature regularization was
(sometimes) helpful.

## making plots

All plot generating scripts are the in paper-figs folder. Most
require that you do the data generation and inverse solver runs
above (slow).

Plots for first revision of paper:
- Introduction: see scattering_plots_imp_and_trans.m
- Section 4.1: see reconstructions_diff_imp.m and
reconstructions_lamabv_withconstinitial.m. For the latter,
you'll have to edit the script according to the instructions
at the top.
- Section 4.2.1: see reconstructions_vary_delta.m and
residual_image_delta_and_freq.m. In each, set image_to_make=421
and run the script to get the reconstructions and error plots
as in the paper.
Also set image_to_make=4212 and run reconstructions_vary_delta.m
to get the recovered parameters as in the paper.
- Section 4.2.2: see reconstructions_vary_delta.m and
residual_image_delta_and_freq.m. In each, set image_to_make=422
and run the script. 
- Section 4.2.3: see reconstructions_vary_delta.m and
residual_image_delta_and_freq.m. In each, set image_to_make=423
and run the script.
- Section 4.2.4: see reconstructions_vary_delta.m and
residual_image_delta_and_freq.m. In each, set image_to_make=424
and run the script and set image_to_make=4242 and run the script.
- Section 4.2.4: see reconstructions_vary_delta.m and
residual_image_delta_and_freq.m. In each, set image_to_make=43
and run the script.

## running MATLAB on a server

nohup matlab -nodisplay -nosplash -r "try; inversetest_runner(test_id); catch; disp('failed'); end; exit" > out.txt &
