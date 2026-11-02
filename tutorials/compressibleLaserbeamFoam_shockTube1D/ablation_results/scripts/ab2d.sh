source /usr/lib/openfoam/openfoam2506/etc/bashrc
cd ~/plumeEdge2D
printf '%s\n' P:full P:full:mm P:full:nc+mm P:full:o10+nc+mm P:contErr:o10+nc+mm \
  | xargs -P 4 -I{} sh -c 'NP=2 END=2e-8 WRITE=1e-9 ./Allrun {} > log.allrun.$(echo {} | tr ":+" "_p") 2>&1'
cat log.allrun.P_full* log.allrun.P_contErr_o10pncpmm
