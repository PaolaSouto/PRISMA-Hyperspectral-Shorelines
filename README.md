# PRISMA-Hyperspectral-Shorelines
These codes allow to extract subpixel shorelines with two different methods from PRISMA hyperspectral images as described in the publication: "Subpixel shoreline mapping tool from PRISMA hyperspectral remotely-sensed images.”

https://doi.org/10.3390/rs15082117

## Data

The codes work with PRISMA hyperspectral images distributed by the "Agenzia Spaziale Italiana" (ASI) after request. For the image license request access to the link https://www.asi.it/scienze-della-terra/prisma/ and click on "portale PRISMA" as indicated with the red arrow on the image below.

<img src="https://github.com/PaolaSouto/PRISMA-Hyperspectral-Shorelines/blob/main/imgs/AccessPrisma.jpg">

A short project of, at least, 150 words must be sent, and ASI will release the license within a few days. 

## 1. Installation

### Create Anaconda environment

Several steps are required before running the codes.

1. Install *Anaconda*, which can be download from https://www.anaconda.com/download.

Using the shell (Mac and Linux):

2. Go to the folder where this repository was downloaded by tipping: 'cd repository_directory'
3. Create an environment containing the required Python packages: `conda create --name PRISMA_sds --file prisma_sds.txt
4. Activate the environment: `conda activate PRISMA_sds`
5. Launch jupyter notebook: `jupyter notebook`

### Before running the jupyter notebooks

Check the pre-requirements.

MANDATORY :warning:

* The folder structure (SCENES folder)
* The PRISMA images 
* Beach baseline
* Profiles
* Image grid

To set up the folders and files for obtaining the SDS, run the **Preprocessing.ipynb** jupyter notebook.

OPTIONAL

* Pre-defined files
* Bounding Box file
* Nothing else

### Run the jupyter notebooks

Shorelines can be obtained using two methods, each with its own jupyter notebook: the profiles and the k-means approach.


* #### Profiles approach: Profiles_method.ipynb




* #### K-means approach: kmeans_method.ipynb




### Outputs

The results for each algorithm were one .txt and one .geojson with the Satellite-Deriv⍺ed Shoreline of the analyzed image. The .txt is given in image coordinates, whereas the geojson is in real-world coordinates.


The Satellite-Derived-Shorelines are stored in the directory ../SCENES/COUNTRY/BEACH_NUM/Results/SDS/ with each method having its own folders.



