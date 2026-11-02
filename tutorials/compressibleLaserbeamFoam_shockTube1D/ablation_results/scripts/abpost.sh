source /usr/lib/openfoam/openfoam2506/etc/bashrc
OD="/mnt/c/Users/g05228sr/OneDrive - The University of Manchester/Third Year/LIBs LaserbeamFoam Project/Simulation Files/shockTube_work/ablation_results"
mkdir -p "$OD/2d"
A=~/ablation; rm -rf $A; mkdir -p $A/runs; cd $A
M=~/shockTube_mm/runs
for c in A Af B; do
  ln -s $M/${c}_full runs/${c}_V0; ln -s $M/${c}_full_mm runs/${c}_V1; ln -s $M/${c}_full_nc+mm runs/${c}_V2
  ln -s $M/${c}_full_o10+nc+mm runs/${c}_V3; ln -s $M/${c}_contErr_o10+nc+mm runs/${c}_V4
done
python3 -I ~/shockTube_mm/energy_check.py $(ls -d runs/*) > /dev/null
cp energy_table.txt "$OD/"
sed -n '/^Table 1:/,/^$/p;/^Table 1b/,/^err E_int/p' energy_table.txt
echo "== negative c_i:"; for r in runs/*; do n=$(cat $r/[0-9]*/c.air $r/[0-9]*/c.metal1vapour 2>/dev/null | grep -c -E '^-[0-9]'); echo "  $(basename $r): $n"; done
# 2D: cell centres for plotting
cd ~/plumeEdge2D/runs
for r in P_full P_full_mm P_full_nc+mm P_full_o10+nc+mm P_contErr_o10+nc+mm; do
  (cd $r; postProcess -funcs "(writeCellCentres writeCellVolumes)" > log.geom 2>&1; echo "$r geom $?: $(grep ClockTime log.solver | tail -1)"; grep -m1 PROBLEM log.solver)
done
