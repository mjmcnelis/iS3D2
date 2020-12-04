echo
echo 'Running smooth particle spectra for modified distribution paper......'
echo
echo 'Bulk viscosity =' $1
echo 'df mode =' $2
echo 'viscous correction =' $3
echo

cd ../..


# copy surface and iS3D parameters
cp tests/modified_distribution/noncentral/$1_bulk/input/surface.dat input
cp tests/modified_distribution/noncentral/$1_bulk/parameters/$2/$3/iS3D_parameters.dat .


# clear results and run iS3D
sh clear_results.sh
./iS3D.e

# copy results to tests
cp results/continuous/vn* tests/modified_distribution/noncentral/$1_bulk/data/$2/$3

# $1 = (small, large)
# $2 = (grad, ce, ptm, ptb, fa, famod)
# $3 = (none, shear, bulk, shear_bulk)	for $2 = (grad, ce, ptm, ptb)
#	   (none, shear) 					for $2 = (fa, famod)
