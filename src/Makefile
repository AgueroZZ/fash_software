# Compile all shared objects for TMB models
all: Gaussian_ind.so Gaussian_dep.so Gaussian_ind_fixed.so Gaussian_dep_fixed.so Poisson_ind.so Poisson_ind_fixed.so

# Compile Gaussian_ind.cpp into Gaussian_ind.so
Gaussian_ind.so: Gaussian_ind.cpp
	Rscript --vanilla -e "TMB::compile('Gaussian_ind.cpp', '-O0 -g')"

# Compile Gaussian_ind_fixed.cpp into Gaussian_ind_fixed.so
Gaussian_ind_fixed.so: Gaussian_ind_fixed.cpp
	Rscript --vanilla -e "TMB::compile('Gaussian_ind_fixed.cpp', '-O0 -g')"

# Compile Gaussian_dep.cpp into Gaussian_dep.so
Gaussian_dep.so: Gaussian_dep.cpp
	Rscript --vanilla -e "TMB::compile('Gaussian_dep.cpp', '-O0 -g')"

# Compile Gaussian_dep_fixed.cpp into Gaussian_dep_fixed.so
Gaussian_dep_fixed.so: Gaussian_dep_fixed.cpp
	Rscript --vanilla -e "TMB::compile('Gaussian_dep_fixed.cpp', '-O0 -g')"

# Compile Poisson_ind.cpp into Poisson_ind.so
Poisson_ind.so: Poisson_ind.cpp
	Rscript --vanilla -e "TMB::compile('Poisson_ind.cpp', '-O0 -g')"

# Compile Poisson_ind.cpp into Poisson_ind.so
Poisson_ind_fixed.so: Poisson_ind_fixed.cpp
	Rscript --vanilla -e "TMB::compile('Poisson_ind_fixed.cpp', '-O0 -g')"

# Clean up generated files
clean:
	rm -f *.so *.o
