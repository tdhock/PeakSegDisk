#!/bin/bash
cd ..
set -o errexit
rm -rf PeakSegDisk-release
cp -r PeakSegDisk PeakSegDisk-release
grep -v Remotes PeakSegDisk/DESCRIPTION > PeakSegDisk-release/DESCRIPTION
rm PeakSegDisk-release/tests/testthat/*
cp PeakSegDisk/tests/testthat/test-CRAN*.R PeakSegDisk-release/tests/testthat
PKG_TGZ=$(R CMD build PeakSegDisk-release|grep building|sed 's/.*‘//'|sed 's/’.*//')
R CMD INSTALL $PKG_TGZ
R CMD check --as-cran $PKG_TGZ
