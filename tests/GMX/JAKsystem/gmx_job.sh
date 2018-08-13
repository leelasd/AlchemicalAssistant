GMX_BIN='/usr/local/gromacs/bin/'
ifile=$1
$GMX_BIN/gmx grompp -f em.mdp -c ${ifile}.pdb -p topol.top -o em.tpr -maxwarn 1
$GMX_BIN/gmx mdrun -deffnm em
echo '1 2 3 4 5 6 7 8 9 0' | $GMX_BIN/gmx energy -f em.edr > LOG_GMX
