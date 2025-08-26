# LAMMPS input script for Carbon (diamond) substrate
# ----------------------------------------------------

# STEP 1: Initialization
# ----------------------
units           metal
boundary        p p p
atom_style      atomic

# STEP 2: System Definition
# -------------------------
lattice         diamond 3.567
region          box block 0 10 0 10 0 10
create_box      1 box
create_atoms    1 box

# STEP 3: Force Field and Minimization
# ------------------------------------
pair_style      tersoff
pair_coeff      * * SiC.tersoff C

mass            1 12.011

minimize        1.0e-5 1.0e-7 1000 10000

# STEP 4: Equilibration
# ---------------------------------------
velocity        all create 300.0 4928459 mom yes rot yes
fix             1 all nvt temp 300.0 300.0 0.1
timestep        0.001
thermo          500
thermo_style    custom step temp pe etotal press
dump            1 all atom 1000 dump.carbon.atom
run             10000

# ----------------------------------------------------
# End of Script