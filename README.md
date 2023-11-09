# Impedance Examples

Scripts that generate data, solve inverse problems, and make plots
which appear in the paper "Inverse obstacle scattering with impedance as a
model for transmission data"

## dependencies/reproducibility

Depends on the inverse-obstacle-scattering-2d MATLAB package, which
is available at https://gitub.com/fastalgorithms/inverse-obstacle-scattering2d.
The paper examples were run using the files in the devel-impedance branch
at commit https://github.com/flatironinstitute/inverse-obstacle-scattering2d/commit/861a3de5fff11ba0b31ccd4c819b731ee6137065

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

Data for version 1 of paper:
- Section 4.1 (impedance): test_ids [16,17,18,19]
- Section 4.2.1 (transmission): test_ids [59,60,61,62,63]
- Section 4.2.2 (transmission): same as 4.2.1
- Section 4.2.3 (transmission): test_ids [66,67,68,69,70]
- Section 4.2.4 (transmission): same as 4.2.1 plus test_ids [64,65]

## solving inverse problems

The inverse problems can be run in the imp-data/tests and trans-data/tests
folders. These folders contain a function called inversetest_runner
which will run the inverse solver with various default settings.

Tests for version 1 of paper:
- Section 4.1 (impedance): for each data set, do
inversetest_runner(test_id) [ABV model],
inversetest_runner(test_id,true,[],[],[],true) [CH model with constraints],
and inversetest_runner(test_id,[],[],true) [Fourier model]
- Section 4.2.1 (transmission): for each data set, do
inversetest_runner(test_id)
- Section 4.2.2 (transmission): for each data set, do
inversetest_runner(test_id,[],[],[],[],[],[],[],1e-1)
- Section 4.2.3 (transmission): for each data set, do
inversetest_runner(test_id)
- Section 4.2.4 (transmission): for each data set, do
inversetest_neumann_runner(test_id)

## making plots

All plot generating scripts are the in paper-figs folder. Most
require that you do the data generation and inverse solver runs
above (slow).

Plots for version 1 of paper:
- Introduction: see scattering_plots_imp_and_trans.m
- Section 4.1: see reconstructions_diff_imp.m
- Section 4.2.1: see reconstructions_vary_delta.m and
residual_image_delta_and_freq.m. In each, set image_to_make=421
and run the script. 
- Section 4.2.2: see reconstructions_vary_delta.m and
residual_image_delta_and_freq.m. In each, set image_to_make=422
and run the script. 
- Section 4.2.3: see reconstructions_vary_delta.m and
residual_image_delta_and_freq.m. In each, set image_to_make=423
and run the script. 
- Section 4.2.4: see reconstructions_vary_delta.m and
residual_image_delta_and_freq.m. In each, set image_to_make=424
and run the script. 
