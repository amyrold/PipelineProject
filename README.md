# PipelineProject
# Packages needed:
1. Bowtie2
2. SPAdes
3. Biopython (for entrez)
4. python import of "os" (done via script)
5. python import of subprocess (done via script)

# Comments about pipeline
The whole script is self-contained. PipelineProject.py can be run from any directory and will be able to create necessary folders, download raw or test data, and perform required analyses. There is test data provided in the github repo with appropriate relative pathing, whether it is downloaded or not, the pipeline will run as expected. When you run the file, it will prompt "raw" or "test" to determine which data to download/use. From there, everything else will be performed automatically. 

# Prep:
- clone github repo
- if data/test data are being downloaded fresh, the PipelineProject.py file can be run from any diretory
- if the user is using the pre-downloaded test data, the test data must be inside a folder (5_data_test) within the same directory as PipelineProject.py (otherwise the pipeline will download the test data itself).

# Code usage:
- Use "python3 PipelineProject.py" to run the .py file
- you will then be prompted for "test" or "raw" (choose test)
- everything else will be carried out by the script. 
- the resulting PipelineProject.log file will be produced inside the same directory as the .py file
