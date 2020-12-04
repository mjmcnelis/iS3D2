
cd ../..

# copy chosen particles and momentum tables
cp tests/modified_distribution/chosen_particles.dat PDG
cp tests/modified_distribution/central/tables/pT_table.dat tables/momentum
cp tests/modified_distribution/central/tables/phi_table.dat tables/momentum
cp tests/modified_distribution/central/tables/y_table.dat tables/momentum
cp tests/modified_distribution/central/tables/eta_table.dat tables/spacetime_rapidity


# compile iS3D
sh cleanMakeCPU.sh

cd scripts/modified_distribution

sh run_modified_distribution_central.sh small grad none
sh run_modified_distribution_central.sh small grad shear
sh run_modified_distribution_central.sh small grad bulk
sh run_modified_distribution_central.sh small grad shear_bulk

sh run_modified_distribution_central.sh small ce shear
sh run_modified_distribution_central.sh small ce bulk
sh run_modified_distribution_central.sh small ce shear_bulk

sh run_modified_distribution_central.sh small ptm shear
sh run_modified_distribution_central.sh small ptm bulk
sh run_modified_distribution_central.sh small ptm shear_bulk

sh run_modified_distribution_central.sh small ptb shear
sh run_modified_distribution_central.sh small ptb bulk
sh run_modified_distribution_central.sh small ptb shear_bulk

sh run_modified_distribution_central.sh large grad none
sh run_modified_distribution_central.sh large grad shear
sh run_modified_distribution_central.sh large grad bulk
sh run_modified_distribution_central.sh large grad shear_bulk

sh run_modified_distribution_central.sh large ce shear
sh run_modified_distribution_central.sh large ce bulk
sh run_modified_distribution_central.sh large ce shear_bulk

sh run_modified_distribution_central.sh large ptm shear
sh run_modified_distribution_central.sh large ptm bulk
sh run_modified_distribution_central.sh large ptm shear_bulk

sh run_modified_distribution_central.sh large ptb shear
sh run_modified_distribution_central.sh large ptb bulk
sh run_modified_distribution_central.sh large ptb shear_bulk




