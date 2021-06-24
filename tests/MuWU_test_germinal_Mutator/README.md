Tests include all input files.  
Options are already pre-selected for *Mutator* TE and this particular test dataset.  

To run enter:  
`snakemake --use-conda --cores xx`  
`snakemake --use-conda --cores xx --conda-prefix conda_envs` (when using the singularity container)  

This test dataset simulates a run using the GRID method which also identifies the subset of heritable germinal insertions.  

Note that this particular test will build a bowtie2 index for the complete maize v4 reference genome.  
A complete run takes ~30 min on an AMD EPYC 7662 with 48 cores.  

Exemplary output from `results/insertions_table_final/germinal_identified_insertions_annotated.csv`:  

`stock` information refers to the intersection in the grid matrix and thus the F2 family (in the *BonnMu* project) in which the germinal insertion/mutation can be uniquely found.

```
"GeneID","Chr","Start","End","Sample","InsertionStart","InsertionEnd","StartReads","EndReads","Gene_length","stock"
"Zm00001d027324","1",2959502,2972168,"Row-01",2959609,2959617,2,4,12667,"D-0081"
"Zm00001d027324","1",2959502,2972168,"Col-04",2959609,2959617,3,15,12667,"D-0081"
"Zm00001d028211","1",26525210,26529825,"Row-03",26529757,26529765,16,2,4616,"D-0127"
"Zm00001d028211","1",26525210,26529825,"Col-02",26529757,26529765,2,4,4616,"D-0127"
"Zm00001d028532","1",38393485,38403506,"Row-02",38393512,38393520,2,3,10022,"D-0102"
"Zm00001d028532","1",38393485,38403506,"Col-01",38393512,38393520,2,2,10022,"D-0102"
"Zm00001d029980","1",97331766,97337093,"Row-01",97331843,97331851,2,2,5328,"D-0081"
"Zm00001d029980","1",97331766,97337093,"Col-04",97331843,97331851,6,2,5328,"D-0081"
"Zm00001d030164","1",109844576,109871086,"Row-03",109844709,109844717,2,2,26511,"D-0127"
"Zm00001d030164","1",109844576,109871086,"Col-02",109844709,109844717,3,2,26511,"D-0127"
"Zm00001d030533","1",142260864,142278394,"Row-01",142273323,142273331,5,8,17531,"D-0079"
"Zm00001d030533","1",142260864,142278394,"Col-02",142273323,142273331,6,8,17531,"D-0079"
"Zm00001d031013","1",173396175,173402763,"Row-01",173396182,173396190,4,4,6589,"D-0080"
"Zm00001d031013","1",173396175,173402763,"Col-03",173396182,173396190,3,2,6589,"D-0080"
...
```
