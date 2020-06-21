# Extraction of objective material barriers to the transport of momentum and vorticity in two-dimensional flows

## License

This software is made public for research use only. It may be modified and redistributed under the terms of the GNU General Public License.

## Algorithm

The current scripts implement theoretical results developed by the Haller Group at ETH Zurich. More specifically, the code presented here demonstrates the computation and extraction of objective material barriers to the transport of momentum and vorticity in a two-dimensional decaying turbulence simulation, as they are presented in [1].

## References
[1] G. Haller, S. Katsanoulis, M. Holzner, B. Frohnapfel & D. Gatti, Objective material barriers to the transport of momentum and vorticity. submitted (2020)

## Installation notes

1) After you unzipped the files to mydir, put the Current Directory in Matlab to mydir

2) Download the two dimensional turbulence dataset (velocity and vorticity fields in the folder "Data") and the advected trajectories ("particles1100.mat") for the Lagrangian calculations from: https://polybox.ethz.ch/index.php/s/L0b0piQEqhOqNjg and save them to mydir. The combined size of the entire dataset is approximately 10GB.

3) For a comparison of the classic FTLE/PRA and the aFTLE/aPRA run sequentially every section of the scripts inside the folder "aFTLE-aPRA computation".

4) For the extraction of material barriers to the transport of momentum and vorticity run sequentially every section of the scripts inside the folder "Extraction of vortex boundaries".

Tested on Matlab R2019a.

NOTE: This code may be improved and is subject to several changes. Therefore, we suggest to visit this page and check if the version you downloaded is up to date.  

Maintained by Stergios Katsanoulis,

katsanos at ethz dot ch

June 22, 2020.
