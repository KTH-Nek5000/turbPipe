# Nek5000 v.17 case for simulation and post-processing of fully-developed turbulent pipe flow

## General:
  - Nek5000 version: v.17

## Documentation:
  - See the [technical report](http://urn.kb.se/resolve?urn=urn:nbn:se:kth:diva-265021)

## To compile:
  1. cd compile/
  2. set the settings in makenek and compile_script
  3. ./compile_script --all

## To run:
  1. cd run/
  2. set the number of processors in runBash_local.sh
  3. bash runBash_local.sh

## To generate mesh for the pipe:
  - Use the gmsh-script [PipeMesh](https://github.com/KTH-Nek5000/PipeMesh) 
