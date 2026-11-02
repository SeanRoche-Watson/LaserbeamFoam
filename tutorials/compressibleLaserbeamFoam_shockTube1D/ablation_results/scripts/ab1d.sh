source /usr/lib/openfoam/openfoam2506/etc/bashrc
cd ~/shockTube_mm
# V0-V3 for A/Af/B not already run (V4 = *_contErr_o10+nc+mm and B V0/V1 = B_full, B_full_mm exist)
printf '%s\n' Af:full:o10+nc+mm Af:full:nc+mm Af:full:mm Af:full \
              A:full A:full:mm A:full:nc+mm A:full:o10+nc+mm B:full:nc+mm B:full:o10+nc+mm \
  | xargs -P 4 -I{} sh -c './Allrun {} > log.allrun.$(echo {} | tr ":+" "_p") 2>&1'
for f in log.allrun.A*_full* log.allrun.B_full_nc* log.allrun.B_full_o10pnc*; do echo "== $f"; grep -E "solver exit|PROBLEM" $f; done
