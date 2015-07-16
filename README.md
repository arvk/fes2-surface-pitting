# FeS<sub>2</sub>-surface-kMC #

This is a simple Fortran program for simulating dynamics of point defects (Fe and S vacancies) at different temperatures on the (100) surface of pyrite, FeS<sub>2</sub> using kinetic Monte Carlo. This code acts as supplementary information to the journal article **Dynamics of point defect formation, clustering and pit initiation on the pyrite surface**, *Electrochimica Acta*, 127, 416-426, 2014. DOI: [10.1016/j.electacta.2014.02.048](http://dx.doi.org/10.1016/j.electacta.2014.02.048). The article also provides sources for the barriers used in the code.


Installing and Running
-------

* Compile with the MPI wrapper using `mpif90 -o kmc.exe main.f90`
* Create a folder for each temperature. For the default code in the repo, you should do `mkdir 298K 299K 393K 394K 443K 444K 483K 484K 513K 514K 543K 544K 573K 574K 603K 604K`
* Copy the `param.in` file into each folder 
* Run using `mpiexec -n <no of cores> kmc.exe`


Output
-------

* Defect distributions on the FeS<sub>2</sub> surface are printed out in `fp-*.txt` files after a fixed number of kMC steps. The latest version of the file is stored as `forplot.txt`
* Use the gnuplot script in the `vis` folder to plot out a simulated STM-image of the surface


Repo
-----

This repo is available at Bitbucket (https://bitbucket.org/arvk/fes2-surface-kmc) and github.mit.edu (https://github.mit.edu/aravindk/fes2-surface-kmc)