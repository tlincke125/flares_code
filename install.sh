#!/usr/bin/env sh

######################################################################
# @author      : thli2739 (theo.j.lincke@gmail.com)
# @file        : install
# @created     : Monday Jan 10, 2022 15:02:57 MST
#
# @description : Builds then installs flares 
######################################################################


python3 -m build
pip3 install ./dist/flares-segmentation-tlincke125-0.0.1.tar.gz
