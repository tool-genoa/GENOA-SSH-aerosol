# This is a python tool that converts and reduces MCM chemical schemes to SSH chemical schemes. 

- Author: Zhizhao Wang

- Last update: 16/12/2021

Reduction Options:
------------------

1.	lump species by their functional groups, MWs, saturation vapor pressure.

2.	remove all the further reactions from non-volatile species (Psat < threshold).

3. 	remove species/reactions if the gas/aerosol concentrations of them and their derived species are under the threshold.

4.	remove KDEC reaction.

5. 	remove species/reactions if the minimum Psat of all derived species is still high.

6.	removing MCM reactions if the number of the derived species is under the defined threshold.


 
Required third party python libraries:
-------------------------------------

1.	OpenBabel


2.	netCDF4


3.	UManSysProp




Source code:
------------

1.	__init__.py



2.	Module.py

	store classes (Reaction, Species, Kinetic)


3.	MolProperty.py

	functions related to molecule properties: MWs, functional groups, Psat, Hvap


4.	DataStream.py

	functions that read MCM reaction and species files and write SSH required files 


5.	ReductionStrategy.py

	functions that merge or reduce species/reactions


6.	SSHSetUp.py
	functions that apply the new chemical scheme to SSH snd run simulations in SSH


7.	SSHResultProcess.py

	functions that are used for processing simulation results in SSH


8.	InitConcExtract.py

	functions that extract data from NetCDF file


9.	KineticRate.py

	functions that transfer MCM kinetic rate to SPACK format


10.	files/
	store files that required in the program


Example of input MCM files: FromMCM/
-------------------------------------

1.	reacion list

2.	species list



Example of output SSH files: ToSSH/
------------------------------------

1.	[name].species

	species list for SPACK

2.	[name].reactions

	reactions list for SPACK

3.	[name].aer

	aerosol species-list for SSH

4.	[name].RO2

	RO2 species list for SSH chem()

5.	[name].cst

	list of species with constant conc. during the simulations for SSH 

6.	[name].photolysis

	list of photolysis reactions for SSH 	

