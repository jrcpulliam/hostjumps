***** README.txt
***** JRCP (Last modified: 29 January 2009)

This file accompanies the supplementary files to Pulliam and Dushoff (2009), which are available for download at: 
http://lalashan.mcmaster.ca/hostjumps/

The original paper is available at: http://www.journals.uchicago.edu/doi/full/10.1086/596510

*** CONTENTS ***

This file contains the following information:

    * Supplementary file descriptions
        - Search documentation
        - Data files
        - Code files
    * How to use the supplementary files
        - Search documentation
        - Data files
        - Setting up your R session (a very basic introduction)
        - Code files
    * Copyright and licensing information
    * Contact information

*** SUPPLEMENTARY FILE DESCRIPTIONS ***

Three types of files (not including this readme file) are available for download on the supplementary materials website: search documentation, data files, and code files. The files can be downloaded individually or as a complete package in the form of a zipped archive. A description of each file follows.

  ** Search documentation **

PLEASE NOTE: Rich text files (files with the extension .rtf) can be read in a text editor, but the files 
available here include additional formatting (such as text highlighting) that will only appear if the 
document is opened by programs that support this type of formatting. These files are also available as 
portable document files (files with the extension .pdf), which preserve all formatting but are 
non-modifiable.

  Families and Genera infecting vertebrates.rtf -- This file lists all viral family and genera that are known to infect vertebrates. This document was used to determine which viral groups to conduct searches for in order to determine whether a virus species would be included in the database. For all viral families with members known to infect vertebrates, searches were conducted for all species in any genera with members known to infect mammals, unless all members of a genus require co-infection with another viral species for replication (marked "dependent") or all members of a genus were known a priori to be maintained in humans (marked "human"). 

  Pulliam and Dushoff inclusion searches.rtf -- This file lists all searches conducted to determine which viral species would be included in the database and to determine human infection status for those included in the database. If inclusion criteria were not clear-cut, a note is included to describe the decision that was made and provide justification for that decision. In addition, if searches were not performed for a species, a note is made as to why the species was excluded from searches (eg, the species is dependent on co-infection for replication or the species is maintained in humans).

  ** Data files **

  database.txt -- This file is a tab-delimited data file containing the following information on each viral species known to infect domestic artiodactyls and not excluded based on the inclusion criteria (see Methods section in the text of the paper for details of the criteria used to determine inclusion in the database): Family (the viral family to which the species belongs), Genus (the viral genus to which the species belongs), Species (an abbreviation that uniquely identifies the viral species), Status (the taxonomic status, either "approved" or "tentative" or the viral species, as listed in the 8th Report of the International Committee on the Taxonomy of Viruses), and Human (human infection status, determined according to criteria described in Supplementary Table S1).

  database.with.flu.txt -- This file is the database.txt file with one additional viral species included (Influenza A). This file is used for supplementary analyses to demonstrate the robustness of the results to the inclusion of Influenza A in the database. See the Methods and Results sections of the  associated paper for details.

  vir.fam.txt -- This file is a tab-delimited data file that is used to map the viral traits of interest (segmentation, genomic material, and site of replication) onto viral species within a family. The file contains the following information for each viral family: Family (name of the viral family), Segments (genome segmentation, coded as "multiple" or "single"), Nucleic Acid (genomic material, coded as "RNA" or "DNA"), Site (site of replication, coded as "nucleus" or "cytoplasm"), Seg (genome segmentation, coded as binary), GM (genomic material, coded as binary), and SR (site of replication, coded as binary). See the Methods section of the associated paper for details of coding.

  vir.gen.txt -- This file is a tab-delimited data file that is used to match viral genera to family-level traits in order to perform genus-level permutation tests. This file contains the following information for each viral genus: Family (name of the viral family) and Genus (name of the viral genus).

  ** Code files **

PLEASE NOTE: All code is written in R. R is a statistical programming language and software package that is distributed under a GNU General Public License. R documentation and software is available for free download through the R Project for Statistical Computing website at http://www.r-project.org. The software is available as an executable file for a wide variety of operating systems and computer architectures, and a compilable binary is also available should you need it.

  sgcp.functions.R -- This file defines a number of functions that are used in the two scripts for analysis (sgcp.model.R and sgcp.hypothesis.R).

  sgcp.data.prep.R -- This file is used to import the database into R from the data files and prepare the data for analysis.

  sgcp.model.R -- This file is used for all analyses based on model-based prediction. The file also produces Figures 1 and 2, which appear in the associated paper.

  sgcp.hypothesis.R -- This file is used for all analyses based on hypothesis testing via permutation tests.


*** HOW TO USE THE SUPPLEMENTARY FILES ***

  ** Search documentation **

  The search documentation is provided to allow others to reproduce or update our database and to facilitate follow-up studies. These files are not needed for analysis.

  ** Data files **

  These files contain the data needed for analysis. To run an analysis, you will need both vir.fam.txt and vir.gen.txt and either database.txt or database.with.flu.txt (depending on whether you are running the analysis with or without including Influenza A) to be located in the same directory as the code files used for analysis. You will also need this directory to be the working directory for your R session (see below for instructions).

  ** Setting up your R session (a very basic introduction) **

  Once you have placed the appropriate files into a directory and opened an R session, you need to make sure that R knows where to access the files. To do this, first ask R the location of its current working directory, by typing the following command at the prompt:

    getwd()

If the directory returned by R is different from the directory where your files are stored, you can change the working directory using the setwd() command. To learn about the setwd() command, type the following command at the prompt:

    ?setwd

This will open the help file associated with the setwd() command. You will want to provide the setwd() command with the path of the directory where your files are stored, though exactly how you specify the path will depend on your operating system and the version of R you are using. For example, on a Mac, you might use the following command to access a directory called "Pulliam and Dushoff 2009" in your home directory:

    setwd("~/Pulliam and Dushoff 2009")

On a PC, on the other hand you might use a command that looks something like this to access a directory called "Pulliam and Dushoff 2009" located directly on your C drive:

    setwd("C://Pulliam and Dushoff 2009")

You can use the output returned by the getwd() command to help you determine how to format the path name on your computer.

  ** Code files **

  These files contain all the scripts needed for analysis and for producing graphics. The files sgcp.functions.R and sgcp.data.prep.R are run from within the analysis files (sgcp.model.R and sgcp.hypothesis.R), so they need to be located in the specified work directory, even though you do not need to run them directly.

  Once you have specified the working directory correctly, you can run the provided scripts. To perform model-based prediction analyses and reproduce Figures 1 and 2 from the associated paper, use the following command:

    source("sgcp.model.R")

By default, the figures will be saved as PDF files in the specified working directory and the model output and control parameters used will be saved as an R data file ("sgcp.models.Rdata") in the working directory.

  To perform hypothesis testing analyses, use the following command:

    source("sgcp.hypothesis.R")

By default, 100,000 permutations of the data will be used for each of the 9 permutation tests performed, so be aware that it may take a substantial amount of time for the script to finish running. Also note that the default seed for randomization is the same as was used to produce the results reported in the associated paper, so if you run the script as is, you should reproduce the exact results reported despite the random nature of the permutation test. By default, the p-values from the permutation tests and control parameters will be saved as an R data file (sgcp.pval.Rdata) in the working directory.

  Remember that by default all analyses are performed using the version of the database that excludes Influenza A. To include Influenza A in your analyses, you will need to change the value of the control parameter INCLUDE.FLU from FALSE to TRUE in the appropriate script. Control parameters are generally defined near the beginning of of the scripts and are denoted by the use of capital letters. To change other aspects of the analysis or graphics, you may have to make modifications in the body of the code. For example to change which subset of the database is used for analysis (eg, if you would like to code virus species with serological evidence of infecting humans as able to infect humans), you should modify the appropriate line in the data preparation file (sgcp.data.prep.R). See comments within the files for some suggestions on how to make modifications for further exploration of the data set.


*** COPYRIGHT AND LICENSING INFORMATION ***

¬Â© 2009, Juliet R.C. Pulliam and Jonathan Dushoff. Some Rights Reserved.

All supplementary materials are Copyright Juliet R.C. Pulliam and Jonathan Dushoff 2008 and are made available under a Creative Commons Attribution Noncommercial Share Alike United States 3.0 License. You are free to copy, distribute, and modify this work, provided that you do not use the work for commercial purposes, you credit the authors, and you include a copy of this notice.  If this work is used for academic purposes, the associated paper should be cited:

Pulliam, JRC and Dushoff, J (2009) Ability to replicate in the cytoplasm predicts zoonotic transmission of 
livestock viruses. Journal of Infectious Diseases 199(4): 565-568. DOI: 10.1086/596510

If you modify this work, you may distribute the resulting derivative work only under the same or a similar license to this one and with proper attribution of the original work, as stated above. To see further details of this license, including a link to the legal code, visit http://creativecommons.org/licenses/by-nc-sa/3.0/.


*** CONTACT INFORMATION ***

Please contact Dr. Juliet Pulliam with questions regarding the files provided and their authorized usage:

Juliet Pulliam, PhD
RAPIDD Working Group on Modeling Zoonoses
Fogarty International Center
National Institutes of Health
Bethesda, MD 20892
Phone: (301) 402-5203
Email: pulliamjuliet@mail.nih.gov



