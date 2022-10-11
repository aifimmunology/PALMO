# PALMO (Platform for Analyzing Longitudinal Multi-omics data)

# <a name="introduction"></a> Introduction
PALMO `(Platform for Analyzing Longitudinal Multi-omics data)` is a platform for anayzing longitudinal data from bulk as well as single cell. It allows to identify inter-, intra-donor variations in genes over longitudinal time points. The analysis can be done on bulk expression dataset without known celltype information or single cell with celltype/user-defined groups. It allows to infer stable and variable features in given donor and each celltype (or user defined group). The outlier analysis can be performed to identify techinical/biological perturbed samples in donor/participant. Further, differential analysis can be performed to deciher time-wise changes in gene expression in a celltype.

<br> ![img](https://github.com/aifimmunology/PALMO/blob/data/data/vignette/PALMO-workflow.png) <br>
General workflow and analysis schema of **PALMO**. It can work with longitudinal data obtained from bulk such as clinical, bulk RNAseq, proteomic or single cell dataset from scRNAseq, and scATACseq.

# <a name="library"></a> Install package and load library

To install library, simply run
    
    From CRAN (https://cran.r-project.org/web/packages/PALMO/index.html)
    install.packages("PALMO")
    library("PALMO")
    
    OR
    
    library("devtools")
    install_github("aifimmunology/PALMO")
    library("PALMO")

# <a name="manual"></a> Reference manual

PALMO reference manual can be obtained from [https://cran.r-project.org/web/packages/PALMO/PALMO.pdf](https://cran.r-project.org/web/packages/PALMO/PALMO.pdf)

# <a name="example-main"></a> Tutorials

There are couple of tutorials presented to help users to run PALMO on bulk and single cell data. The examples includes:

* [Tutorial-1: Plasma proteome longitudinal data](https://github.com/aifimmunology/PALMO/blob/main/Vignette-PALMO.pdf#page=3)
* [Tutorial-2: scRNA longitudinal data](https://github.com/aifimmunology/PALMO/blob/main/Vignette-PALMO.pdf#page=11)
* [Tutorial-3: scATAC longitudinal data](https://github.com/aifimmunology/PALMO/blob/main/Vignette-PALMO.pdf#page=23)
* [Tutorial-4: Multi-modal data integration](https://github.com/aifimmunology/PALMO/blob/main/Vignette-PALMO.pdf#page=31)
* [Tutorial-5: Longitudinal COVID dataset (CNP0001102)](https://github.com/aifimmunology/PALMO/blob/main/Vignette-PALMO.pdf#page=33)
* [Tutorial-6: Differential Gene analysis in longitudinal data](https://github.com/aifimmunology/PALMO/blob/main/Vignette-PALMO.pdf#page=41)
* [Tutorial-7: Mouse brain dataset (GSE129788)](https://github.com/aifimmunology/PALMO/blob/main/Vignette-PALMO.pdf#page=43)
* [Tutorial-8: TCRB profiling dataset(GSE156980)](https://github.com/aifimmunology/PALMO/blob/main/Vignette-PALMO.pdf#page=52)
* [PALMO Data Structure]()(https://github.com/aifimmunology/PALMO/blob/main/Vignette-PALMO.pdf#page=56)
* [Other examples](https://github.com/aifimmunology/PALMO/blob/main/Vignette-PALMO.pdf#page=56) 

# <a name="authors"></a> Authors

[Suhas Vasaikar](https://github.com/suhasaii), [Aarthi talla](https://github.com/aarthitallaAI) and [Xiaojun Li](https://github.com/Xiaojun-Li) designed the PALM algorithm. [Suhas Vasaikar](https://github.com/suhasaii) implemented the PALMO package.

# <a name="license"></a> License
PALMO is licensed under the [MIT-License](https://github.com/git/git-scm.com/blob/main/MIT-LICENSE.txt) and Allen Institute Software License [AISL](https://github.com/AllenInstitute/ghinfo/blob/master/LICENSE). Please read [Allen Institute Terms of use](https://alleninstitute.org/legal/terms-use/).

