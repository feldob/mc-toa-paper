Experiments on Test Set Selection (TSP = Test Selection Problem) for paper *Generic and Industrial Scale Many-Criteria Test Selection*.
Due to confidentiality, the system data required to reproduce the results is not available. All experimental results can be found [here](https://zenodo.org/record/7740864#.ZBMIG4DMJH4).

on linux:
=========
install Julia 1.5.3 and Gurobi 9.1 (with License) locally (see `install.sh`).

execute bash script through `./experiment_script.sh`

In case of problems with out of memory errors, [activate swap via e.g. swap file](https://askubuntu.com/questions/33697/how-do-i-add-swap-after-system-installation) (seems easiest):

Take a coffee... or two, the experiments take a couple of days to run.

learning the framework:
=======================
an example has been added as the `example.jl` file. The code is documentend to clarify the use and how to adjust it.
