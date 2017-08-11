pypenelope-emerg-eng-ang
========================

C++ Code to process the 'pe-trajectories.dat' from pypenelope to extract the energy corresponding to the angle of the emerging particles


pyPENELOPE is an open-source software to facilitate the use of the Monte Carlo code PENELOPE and its main program PENEPMA in the field of microanalysis. It consists in a graphical user interface (GUI) to setup materials, geometry, simulation parameters and position of the detectors as well as to display the simulation’s results. An application programming interface (API) is also available to create a large series of simulations using object-oriented programming and to interpret efficiently the results from these simulations.

Download, Homepage and further information:

* http://pypenelope.sourceforge.net/
* https://launchpad.net/pypenelope


This small C++ program 'extract_emerging.cpp' reads the 'pe-trajectories' given from pypenelope (shower) and extracts the direction and energies of the emerging particles.

It generates a file containing the cartesian x/y/z information of the last positions when the particle exits the surface (e.g. backscattered) and its energy.

With this information the exit-angle and the exit energy can be calculated.

Also a file with the exit angles in spherical (phi, theta, r) coordinates is given where the energies in a specific solid angle are summarized (into bins).


Important:
==========

* variable *tilt_angle* must be minus the angle, the surface is tilted. If the particles hit the angle under an angle of e.g. 30 degrees, the variable *tilt_angle* must be set to -30
* the variable *control_particle* is the amount of incident particles. If secondary electrons(SE) are enabled, the variable *must* be set manually to the amount of incident particles by hand. If not the counted amount of incident particles is wrong and the normalization is getting ruined!


description of the outputfiles
==============================

The outputfiles are described below:



*backscatang_cartesian.dat:*
-------------------------
Contains the cartesian coordinates(x,y,z) of the last two positions of a backscattered particle and its energy. With this information the exiting vector out of the sample can be calculated. With this information the energy backscattered (at specifit angles) is gathered.



*info.dat:*
---------
Contains general information about the processing done (e.g. important settings and parameters)


*energy_raw_data.dat:*
-------------------
For every backscattered particle the exiting (backscattered) energy is written in this file. If the binning of the energies is not satisfying the energies can be rebinned with a different binsize.


*energy_probability.dat:*
-----------------------
The energy from 0 to the maximum backscattered energy(which is usually the incident energy) is divided into bins of the amount of (backscattered energy)/(amount of bins). The energies are collected into the bins. This leads to a distribution of the backscatter energy versus the probability. This all is then normalized.


*binned_probability_spherical.dat:*
---------------------------------
The surface of a sphere is divided into bins. Every bin collects the hits going through the surface of the bin(with a specific solid angle). So this file contains the backscattered probability versus the angle theta and phi. The angles are just 'standard spherical coordinates' where theta goes from [0°;180°] and phi from [0°;360°]


*binned_energy_spherical.dat:*
----------------------------
Same as the 'binned_probability_spherical.dat' but here the enery is summarized into bins.


*probability_summarized.dat.dat:*
-------------------------------
Data from 'binned_probability_spherical.dat' is summarized versus the angle phi.  -> Angle theta versus (backscatter) probability


*energy_summarized.dat:*
----------------------
Data from 'binned_energy_spherical.dat' is summarized versus the angle phi. -> Angle theta versus (backscatter) energy


*energy_probability_cut.dat:*
----------------------
generates the backscatter energy vs. probability data (same as *energy_probability.dat*) for a restriced range of the angle theta. Angles between [-pi:pi] and theta:[theta_low:theta_high] are collected here. Variables *theta_low* and *theta_high* must be set here! Angle theta_low must be lower than theta_high


*energy_probability_cut_negative.dat:*
----------------------
Same as *energy_probability.dat* but for angles between [-pi:pi] and theta:[theta_high:pi]. *energy_probability_cut.dat* and *energy_probability_cut_negative.dat* added gives *energy_probability.dat*
