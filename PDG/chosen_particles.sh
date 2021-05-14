
if [ $1 = "pikp" ]; then
	echo "Setting chosen particles to (pi+,K+,p)"
    cp chosen_particles_pikp.dat chosen_particles.dat

elif [ $1 = "urqmd" ]; then
	echo "Setting chosen particles to UrQMD (v3.3+)"
    cp chosen_particles_urqmd_v3.3+.dat chosen_particles.dat

elif [ $1 = "smash" ]; then
	echo "Setting chosen particles to SMASH"
    cp chosen_particles_smash.dat chosen_particles.dat

elif [ $1 = "box" ]; then
	echo "Setting chosen particles to SMASH Box"
    cp chosen_particles_box.dat chosen_particles.dat

fi


# $1 = pikp, urqmd, smash, smash_box