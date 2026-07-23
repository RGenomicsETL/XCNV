#!/bin/sh
set -eu

printf '%s\n' \
  'XCNV is now an R package.' \
  'Install dependencies in your project environment, then run: R CMD INSTALL .' \
  'Annotation and model resources are external; use XCNV::xcnv_resources(path).' \
  'This script does not download resources, install packages, compile bedtools, or change permissions.'
