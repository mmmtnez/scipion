#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# General script for Xmipp-based pre-processing of micrographs: 
#  - downsampling
#  - power spectral density (PSD) and CTF estimation on the micrograph
#
# For each micrograph given, this script will perform 
# the requested operations below.
# For each micrograph a subdirectory will be created
#
# Author: Carlos Oscar Sorzano, July 2011
#

# {begin_of_header}

#{eval} expandCommentRun()

#-----------------------------------------------------------------------------
# {section} Input Data
#-----------------------------------------------------------------------------
# {run}(import_movies) Import Movie Run
ImportRun = ''
# {expert} Window size
""" 
Window size (shifts are assumed to be constant within this window).
"""
WinSize = 150

# Use GPU version
""" Set to true if you want the GPU implementation
"""
DoGPU =True

# {eval} expandExpert()

#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
from protocol_align_movies import *
if __name__ == '__main__':
    protocolMain(ProtAlignMovies)
