#!/bin/bash

# Directories

HOMEDIR=$PWD
BINDIR=$HOMEDIR/bin
RESDIR=$HOMEDIR/res
MATDIR=/home/lello/Tesi/CalibrationMatlab2/Cpp

# Commands

  make remove
  make

  $BINDIR/*.exe

  cp $RESDIR/*.txt $MATDIR/.
