#!/bin/bash
cd ..
set -o errexit
rm -rf PeakSegDisk-release
cp -r PeakSegDisk PeakSegDisk-release
grep -v Remotes PeakSegDisk/DESCRIPTION > PeakSegDisk-release/DESCRIPTION
rm PeakSegDisk-release/tests/testthat/*
cp PeakSegDisk/tests/testthat/test-CRAN*.R PeakSegDisk-release/tests/testthat
R CMD build PeakSegDisk-release | tee build.out
PKG_TGZ=$(grep building build.out|sed "s/.*\(PeakSegDisk.*.tar.gz\).*/\1/")
echo $PKG_TGZ
R CMD INSTALL $PKG_TGZ
R CMD check --as-cran $PKG_TGZ
