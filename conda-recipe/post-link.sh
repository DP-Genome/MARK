#!/bin/bash
# Printed by conda after the package is installed.
cat >> "${PREFIX}/.messages.txt" <<'MSG'

  MARK is installed. One more step is required: install freebayes, MARK's variant caller, in its own
  environment named mark_fb (it cannot share this environment). MARK finds it by itself:

      conda create -n mark_fb -c conda-forge -c bioconda freebayes=1.3.10 -y

  Then: conda activate <this environment> and run MARKLaunch.py, MARK.sh or MARK-I.sh.

MSG
