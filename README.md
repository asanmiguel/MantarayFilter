# MantarayFilter

This code analyzes experimental and simulation results for a microfluidic Manta Ray inspired filtration system. 

**ParticleCounts_25and15.m:** 
Code analyzes particle collection images of experiments using red fluorescent 25 µm and green fluorescent 15 µm particles. Particle counts are obtained by binarizing and counting objects of a certain size within the red and green channels of individual images.  Efficiency and concentration ratio are calculated for each inlet flow rate and saved to a file for later plotting/visualization. 

**ParticleRange_Counts.m:** 
Code analyzes particle collection images of experiments using green fluorescent particles with a size range of 10-29 µm. For each inlet flow rate and sample location (inlet, out 1, out 2), particles are detected, and their size is recorded for later analysis in ParticleRange_Stats.m.

**ParticleRange_Stats.m:** 
The code imports particle range counts (ParticleRange_Counts.m) and calculates efficiencies for particles within a certain size range based on predetermined bin size.  The code uses function binRangeSize.m to bin particle counts prior to efficiency calculation to allow analyzing various bin sizes simultaneously. 

**binRangeSize.m:** 
This function imports particle range counts and bins particle counts based on size. 

**BentLobe_InflectionPts.m; OblongLobe_InflectionPts.m:** 
These codes import velocity field simulation data obtained from Ansys Fluent for analysis.  Using the known coordinates of the geometry, the velocity field can be analyzed anywhere in the device.  The code finds the lobe with the highest outflow velocity. Next, it measures the height of the local max velocity and inflection point in the velocity profile at the respective lobe.  This code loops over every inlet flow rate for the Oblong or Bent Lobe device to predict filtration success.
