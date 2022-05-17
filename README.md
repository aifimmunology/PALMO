# PALMO (Platform for Analyzing Longitudinal Multi-omics data)

# <a name="introduction"></a> Introduction
PALMO `(Platform for Analyzing Longitudinal Multi-omics data)` is a platform for anayzing longitudinal data from bulk as well as single cell. It allows to identify inter-, intra-donor variations in genes over longitudinal time points. The analysis can be done on bulk expression dataset without known celltype information or single cell with celltype/user-defined groups. It allows to infer stable and variable features in given donor and each celltype (or user defined group). The outlier analysis can be performed to identify techinical/biological perturbed samples in donor/participant. Further, differential analysis can be performed to deciher time-wise changes in gene expression in a celltype.

<br> ![img](https://github.com/aifimmunology/PALMO/blob/data/data/vignette/PALMO-workflow.png) <br>
General workflow and analysis schema of **PALMO**. It can work with longitudinal data obtained from bulk such as clinical, bulk RNAseq, proteomic or single cell dataset from scRNAseq, and scATACseq.

# <a name="library"></a> Install package and load library

To install library, simply run
   
    library("devtools")
    install_github("aifimmunology/PALMO")
    library("PALMO")

# <a name="example-main"></a> Tutorials

There are couple of tutorials presented to help users to run PALMO on bulk and single cell data. The examples includes:

* [Tutorial-1: Plasma proteome longitudinal data](https://github.com/aifimmunology/PALMO/blob/main/ReferenceManual-PALMO-v0.1.0.pdf#page=3)
* [Tutorial-2: scRNA longitudinal data](https://github.com/aifimmunology/PALMO/blob/main/ReferenceManual-PALMO-v0.1.0.pdf#page=11)
* [Tutorial-3: scATAC longitudinal data](https://github.com/aifimmunology/PALMO/blob/main/ReferenceManual-PALMO-v0.1.0.pdf#page=21)
* [Tutorial-4: Multi-modal data integration](https://github.com/aifimmunology/PALMO/blob/main/ReferenceManual-PALMO-v0.1.0.pdf#page=29)
* [Tutorial-5: CNP0001102 data (longitudinal COVID dataset)](https://github.com/aifimmunology/PALMO/blob/main/ReferenceManual-PALMO-v0.1.0.pdf#page=31)
* [Tutorial-6: GSE129788 data (Mouse brain)](https://github.com/aifimmunology/PALMO/blob/main/ReferenceManual-PALMO-v0.1.0.pdf#page=39)
* [Tutorial-7: Differential Gene analysis in longitudinal data](https://github.com/aifimmunology/PALMO/blob/main/ReferenceManual-PALMO-v0.1.0.pdf#page=47)


# <a name="authors"></a> Authors

[Suhas Vasaikar](https://github.com/suhasaii), [Aarthi talla](https://github.com/aarthitallaAI) and [Xiaojun Li](https://github.com/Xiaojun-Li) designed the PALM algorithm. [Suhas Vasaikar](https://github.com/suhasaii) implemented the PALM package.

# <a name="license"></a> License
PALM is licensed under the [MIT-License](https://github.com/git/git-scm.com/blob/main/MIT-LICENSE.txt).
