# The GENerator of Reduced Organic Aerosol Mechanisms (GENOA v2.0)

GENOA v2.0 is an algorithm that generates semi-explicit chemical mechanisms from explicit mechanisms, with a specific focus on SOA (Secondary Organic Aerosol) formation.

Compared to GENOA v1.0, GENOA v2.0 adopts a parallel reduction framework to identify the most optimal reductions from competitive candidates, and can reduce chemical mechanisms from multiple aerosol precursors.

- Last update: 2023/07/25

 
Requirements:
--------------

1.	[python 3.5 or later](https://www.python.org/)


2.	[numpy 1.11.0 or later](https://numpy.org/)


3.	All the requirements to run [ssh-aerosol](https://sshaerosol.wordpress.com/), including the construction tool [SCONS](http://www.scons.org/wiki/SconsTutorial1)

4.	[Open Babel 3.0.1 or later](http://openbabel.org/) (Optional: used for computing aerosol properties from SMILES structures)

5.	[UManSysProp](https://github.com/waveform-computing/umansysprop) (Optional: used for computing aerosol properties from SMILES structures. Requires an update to be used with Python 3)

6.	[matplotlib 1.5.1 or later](https://matplotlib.org/) (Optional: used for postprocessing)

7.	[basemap 1.2.1 or later](https://matplotlib.org/basemap/) (Optional: used for postprocessing)
 




Learning more about the GENOA algorithms:
--------------

- GENOA v1.0


Refer to the GENOA_v1.0_Manual and [Wang et al.,2022](https://doi.org/10.5194/gmd-15-8957-2022) to learn more about GENOA v1.0.

[![GENOAv1.0 code archive](https://zenodo.org/badge/481260565.svg)](https://zenodo.org/badge/latestdoi/481260565)


- GENOA v2.0

The GENOA_v2.0_Manual will be made available soon.

