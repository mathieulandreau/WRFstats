# Article 1 : Observation and LES modeling of an offshore atmospheric undular bore generated during sea-breeze initiation near a peninsula

This folder contains input files from the simulation label 14\_20200515 t30, which is the one used for this article. It also contains Python postprocessing notebook to draw the figures

## input

Contains
- **namelist.wps** : namelist of wps
- **namelist.input** : namelist of WRF 
- **namelist_rst_real.input**, **namelist_rst.input** : WRF namelist of the restart run. The simulation has been run in two parts (first until May 18 and then until May 18)
- **LAI_modif.py** : This script is called between ./real.exe and ./wrf.exe to add the higher resolution LAI and VEGRA field (see article 1)
- **slurm.job**, **slurm_LAI.job**, **slurm_real.job** : sbatch programs to run on the cluster
- **tslist** : list of the timeseries locations
- **wrfio_d0X.txt** : list the desired and undesired output variables

How to use :
- link namelist.wps in the WPS working directory
- run WPS
- link all the other files in the WRF working directory
- run real
- run LAI modif
- run WRF

## postproc

Contains
- **1430_18_article1_GW_carac2.ipynb**
- **1430_18_article1_GW_stats4-5_v2.ipynb**
- **1430_18_article1.ipynb** : Figures 
- **1430_19_GW_var.ipynb**

How to use
- These scripts need the py\_wrf\_arps library to run, the version of py\_wrf\_arps is written in the first cell

