# Spectral flow cytometry-based intracellular cytokine staining

This is a detailed description of the data analysis pipeline for spectral flow cytometry-based intracellular cytokine staining.

## 1) Getting started

### 1.1) Connect this repository to a local directory

1) Scroll to the top of this page and click the blue "Code" botton
2) Copy the link under "Clone with SSH"
3) Open your terminal and navigate to your local folder where you want to clone the git repository to (using the cd and ls commands)
4) Type `git clone [copied link from step 2]` &rarr; Press enter

You have now cloned the repository to your computer. 

Now it's time to open the project in R studio:

### 1.2) Open project in RStudio

5) Open R studio
6) Press "File" &rarr; "New Project..." &rarr; "Existing Directory" &rarr; "Browse" &rarr; choose the local folder from step 3 
7) Press "Create Project"

You have now opened the project in R studio and are ready to proceed with data analysis.

But first, please read the good practices below to ensure quality and safety:

## 2) Good practices

1) Do **not** share data on GitLab, only code
2) Do **not** share personal information (e.g. AUID in path to your data files - instead, use relative paths)
3) To best fulfill point 1 and 2, create a `0_data` folder on your computer which is **not** uploaded to GitLab (add the folder to a `.gitignore` file; see section 2.1)
4) Follow the core principles below:

![alt text](assets/project_structure.png){ width=70% }

**Figure 1: Core principles of project structure for coherent data analysis pipelines**. Created with BioRender.com. Adapted from https://towardsdatascience.com/how-to-keep-your-research-projects-organized-part-1-folder-structure-10bd56034d3a

### 2.1) Adding files or folders to `.gitignore` file
1) In the terminal, navigate to your project folder with the cloned repository from step 3 above (using the cd and ls commands)
2) Check if you have a .gitignore file by writing `ls -a` 
3) If no .gitignore file is listed, create a .gitignore file by writing `touch .gitignore`, else skip this step
4) Check the contents of the .gitignore file by writing `cat .gitignore`
5) If the file/folder you want not to be ignored by Git (i.e. not uploaded) is not present in the .gitignore file, write `echo 'FILENAME/FOLDERNAME' >> .gitignore`
6) Check contents of the .gitignore file (see step 4) to confirm the filename or foldername was added successfully

Now anything listed in the .gitignore file will be ignored by GitLab and therefore not uploaded. 

## 3) Data analysis

Below is a description of the analysis pipeline. The pipeline is inspired by:

- *den Braanker H, Bongenaar M and Lubberts E (2021) How to Prepare Spectral Flow Cytometry Datasets for High Dimensional Data Analysis: A Practical Workflow. Front. Immunol. 12:768113. doi: 10.3389/fimmu.2021.768113*

- *Van Gassen S, Callebaut B, Van Helden M, Lambrecht B, Demeester P, Dhaene T, Saeys Y (2015). “FlowSOM: Using self-organizing maps for visualization and interpretation of cytometry data.” Cytometry Part A, 87(7), 636-645. https://onlinelibrary.wiley.com/doi/full/10.1002/cyto.a.22625.*


### Install packages
**File to use**: `1_scripts/install_packages.R`

If you have never run spectral flow cytometry data analysis on your device before, you need to install some of the most essential R packages through **BiocManager**, **devtools** and **remote**. To do so, use the **install_packages.R** script.
For installation of conventional R packages from a repository, use the `install.packages()` function. 

### Step A: Load data
**File to use**: `1_scripts/1_load_transform.Rmd`

The first step of data analysis is importing the unmixed FCS files into R. The **flowCore** package provides basic structures for working with flow cytometry data in R and imports FCS files as a **flowSet**. 
A flowSet is a class for storing flow cytometry raw data from quantitative cell-based assays (i.e., data from multiple .fcs files). A flowSet consist of one or more **flowFrames**. 
A flowFrame is a class for storing observed quantitative properties for a population of cells from a FACS run (i.e., data from one .fcs file). There are three parts of the data in a flowFrame:

1.	A numeric matrix of the raw measurements with rows=events and columns=parameters
2.	Annotation for the parameters (e.g., measurement channels, stains, dynamic range)
3.	Additional annotation provided through keywords in the .fcs file

![alt text](assets/flowSet.png)

**Figure 4: Import FCS files as a flowSet in RStudio**. Created with BioRender.com. Adapted from https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-106. 

### Step B: Data transformation
**File to use**: `1_scripts/1_load_transform.Rmd`

As the physical process of exciting, emitting and detecting fluorescence signals on a flow cytometer results in increased variance of fluorescence signal when mean fluorescence intensity increases (signal variance is inhomogenous), transformation of data is needed to stabilise the variance. Data is transformed using *logicle transformation* from the **flowWorkspace** package.


### Step C: Pregating
**File to use**: `1_scripts/2_pregate.R`

The second step of data analysis is pregating to remove debris, doublets, and dead cells. A fourth and optional gating can be done on your cell type of interest (e.g. CD45+ positive leukocytes):

![alt text](assets/gatingscheme_pregating.png)

**Figure 5: Gating strategy**. 

Gating is done using the **CytoExploreR** package with interactive pop-up windows to draw you gates. CytoExploreR need a **gatingSet** object as input, which can be created directly from the transformed flowSet. Drawn gates are saved in a **gatingTemplate** (a csv file containing information on parent population, gate name, gate coordinates, etc.) which can be reused on other samples. 

After gating, the end gate population is saved as a flowSet and used for downstream analysis. 

### Step D: Quality control 
**File to use**: `1_scripts/3_qualitycontrol.R`

Automatic quality control is performed using the **PeacoQC** package for removal of low-quality events, e.g.

- Temporary shift in signal (clogs/slow uptake)
- Permanent shift in signal (speed/flow rate change)
- Monotonic decrease in signal (contamination)
- Monotonic increase in signal (contamination)

An overview of the PeacoQC algorithm is shown below:

![alt text](assets/peacoQC.jpg)

**Figure 6: Chart of the PeacoQC algorithm**. From: https://pubmed.ncbi.nlm.nih.gov/34549881/. 

The PeacoQC algorithm created QC fcs files, which can be read back into RStudio as a flowSet for downstream analysis.

### Step E: Gate cell types
**File to use**: `1_scripts/4_gateCellTypes.Rmd`

![alt text](assets/gatingscheme_celltype.png){ width=80% }

**Figure 7: Gating strategy of manual gating**

![alt text](assets/gatingtree_celltypes.png){ width=80% }

**Figure 8: Gating tree of manual gating**. 

### Step F: Manual gating strategy
**File to use**: `1_scripts/5_gateMarkers.R`


### Step G: Clustering and dimensionality reduction
**File to use**: `1_scripts/6_dimred_clustering.Rmd`

