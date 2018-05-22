### paPAM
paPAM is an acoustic particle motion analysis program for fish researchers which is written in Matlab.
- For more information on how and why to analyse particle motion, you can read [this paper](https://www.researchgate.net/publication/292677434_Particle_motion_The_missing_link_in_underwater_acoustic_ecology) as a good starting point.
- For information on how to use paPAM or to get a better idea of what it can do, please check the [user manual](https://raw.githubusercontent.com/RTbecard/paPAM/master/User%20Manual/User%20Manual.pdf).

If you happen to find any bugs in paPAM, please report them [here](https://github.com/RTbecard/paPAM/issues).  Pull requests and general suggestions for improving paPAM are also welcome!

### Compiled binaries

If you want to use paPAM but do not have a version of Matlab installed, you can download the following packages to run paPAM as a standalone compiled program.  In fact, we reccomend using the compiled program rather than the Matlab scripts, as this avoids any issues which are related to running Matlab scripts which were created on a different version.  Please check the user manual for instructions on how to download and install Matlab Compiler Runtime, a free program which is necessary to run the compiled versions of paPAM.

- [Linux (0.87)](https://github.com/RTbecard/paPAM/raw/master/Compiled%20Binaries/MCR_Linux_0.87.zip)
- [Mac (0.87)](https://github.com/RTbecard/paPAM/raw/master/Compiled%20Binaries/MCR_Mac_0.87.zip)
- [Windows (0.901)](https://github.com/RTbecard/paPAM/raw/master/Compiled%20Binaries/MCR_PC_0.901.zip)

### About us

paPAM was created through a collaborative effort between individuals from the following institutions.

<img src="User Manual/Uni_logo.png" width="500" />

### Version History
- 0.901: Fixed bug where "skip first seconds" is set to 0, causing batch analysis to fail.
