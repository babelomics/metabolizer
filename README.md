# The Metabolizer

Here you can find the R script of the Metabolizer tool and some toy examples of its usage. 
In addition to this R code, you can also perform the analysis of the metabolic module activity using the Metabolizer webserver at: http://metabolizer.babelomics.org/

Here you can find the pathway modules of approx. 3400 organisms that are ready for analysis by Metabolizer. 
https://drive.google.com/file/d/1-k1r5MPqz1yknOnOu1p73t_nzaZrfdSp/view?usp=sharing

For each organism you will find: 

**->** org_module_data_April2019.RData (modules in graph format; this file is mandatory for the analysis and used by *metabolizer* function through *moduleinfo* argument. Check https://github.com/babelomics/metabolizer/blob/develop/metabolizer.git/script/analysis/example.R)

**->** org_geneIDs.RData (for gene id conversion; this file is mandatory for the analysis)

> org_KEGG_module_and_path_table_April2019.txt and org_moduleinfo.RData (detailed table of the modules for this organism; this file is mandatory for the analysis)
https://drive.google.com/file/d/1-k1r5MPqz1yknOnOu1p73t_nzaZrfdSp/view?usp=sharing

Please, use the following convention to cite Metabolizer tool:
Cubuk, C., Hidalgo, M., Amadoz, A., Rian, K., Salavert, F., Pujana, M., Mateo, F., Herranz, C., Carbonell-Caballero, J. and Dopazo, J. (2019). Differential metabolic activity and discovery of therapeutic targets using summarized metabolic pathway models. npj Systems Biology and Applications, 5(1), doi: 10.1038/s41540-019-0087-2.
https://www.nature.com/articles/s41540-019-0087-2

You might also be interested in the following study, do not miss it :)
Cubuk, C. et al. Gene expression integration into pathway modules reveals a pan-cancer metabolic landscape. Cancer Res. 2705, 2017 (2018).
http://cancerres.aacrjournals.org/content/78/21/6059

