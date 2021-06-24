Tests include all input files.  
Options are already pre-selected for karma TE and this particular test dataset - see [RelocaTE2](https://github.com/JinfengChen/RelocaTE2).  

To run enter:  
`snakemake --use-conda --cores xx`  
`snakemake --use-conda --cores xx --conda-prefix conda_envs` (when using the singularity container)  

karma transposons in this dataset result in 12 bp TSDs.  

The RelocaTE2 test data is untargeted and in fact many different TEs are found and analyzed in the analysis by RelocaTE2.  

In this test we focus on karma transposons since these create the longest TSDs of all TEs in the test data and we wanted to test MuWU with TSD targets > 9 bp (*Mutator*).  

The test runs the MuWU pipeline on the test dataset and trims the karma consensus sequence from all reads if it is found.  
This way, a subset of reads is created that starts/ends with the TSD region - creating a library suitable for the MuWU insertion identification algorithm.  

Finally, we compare the MuWU found karma insertions with those reported by RelocaTE2.

In terms of recall we identify 21 of 23 RelocaTE2 insertions.  
As the input data is not targeted towards karma TEs we have overall little read support and thus reduced the `overlap_support` parameter which results in many additional potential insertion sites.  

In conclusion, MuWU can easily be adapted to other TEs - however a targeted approach is required for specfic TSD coverage which allows increasing the `overlap_support` parameter which in turn leads to a final set of higher confidence and with less false positives.

See our test with [ITIS test data](https://github.com/tgstoecker/MuWU/tree/master/tests/MuWU_ITIS_test_with_mping_TE) which is derived from the genome resequencing project of Japonica A123(SRR631734), which have TE mping activated. The library is thus enriched for mping TE making MuWU more suitable.  
We show that our top candidates are the same as the high confidence set as reported by [ITIS](https://github.com/Chuan-Jiang/ITIS).  

Output generated by the test:

```
RelocaTE2 finds 23 karma insertions in this dataset:

     V1   V2                             V3      V4      V5 V6 V7 V8
1  Chr3 rice transposable_element_attribute  109857  109868  +  .  . 
2  Chr3 rice transposable_element_attribute  155951  155962  -  .  .
3  Chr3 rice transposable_element_attribute  190294  190305  -  .  .
4  Chr3 rice transposable_element_attribute  355582  355593  -  .  .
5  Chr3 rice transposable_element_attribute  400634  400645  +  .  .
6  Chr3 rice transposable_element_attribute  481214  481225  +  .  .
7  Chr3 rice transposable_element_attribute  513469  513480  -  .  .
8  Chr3 rice transposable_element_attribute  667120  667131  -  .  .
9  Chr3 rice transposable_element_attribute  676498  676509  -  .  .
10 Chr3 rice transposable_element_attribute  775878  775889  +  .  .
11 Chr3 rice transposable_element_attribute 1180074 1180085  -  .  .
12 Chr3 rice transposable_element_attribute 1180304 1180315  -  .  .
13 Chr3 rice transposable_element_attribute 1217980 1217991  -  .  .
14 Chr3 rice transposable_element_attribute 1222316 1222327  -  .  .
15 Chr3 rice transposable_element_attribute 1472716 1472727  +  .  .
16 Chr3 rice transposable_element_attribute 1574112 1574123  -  .  .
17 Chr3 rice transposable_element_attribute 1581467 1581478  +  .  .
18 Chr3 rice transposable_element_attribute 1680003 1680014  +  .  .
19 Chr3 rice transposable_element_attribute 1696136 1696147  -  .  .
20 Chr3 rice transposable_element_attribute 1804303 1804314  +  .  .
21 Chr3 rice transposable_element_attribute 1835580 1835591  -  .  .
22 Chr3 rice transposable_element_attribute 1880399 1880410  -  .  .
23 Chr3 rice transposable_element_attribute 1892934 1892945  +  .  .
                                                                                                    V9
1          ID=Chr3.109868.spanners;avg_flankers=3;spanners=0;type=homozygous;TE=karma;TSD=AGTATACATGGT
2    ID=Chr3.155962.spanners;avg_flankers=4;spanners=0;type=homozygous;TE=karma/RIRE3;TSD=CTTTTTCTGCTA
3          ID=Chr3.190305.spanners;avg_flankers=5;spanners=0;type=homozygous;TE=karma;TSD=TGACTTTTTAAA
4          ID=Chr3.355593.spanners;avg_flankers=5;spanners=0;type=homozygous;TE=karma;TSD=TCGATCGCCTCC
5          ID=Chr3.400645.spanners;avg_flankers=3;spanners=0;type=homozygous;TE=karma;TSD=TTTCTCTAATTC
6        ID=Chr3.481225.spanners;avg_flankers=3.5;spanners=0;type=homozygous;TE=karma;TSD=TTAAACTGAGGG
7        ID=Chr3.513480.spanners;avg_flankers=4.5;spanners=0;type=homozygous;TE=karma;TSD=TCACGAAACAAT
8          ID=Chr3.667131.spanners;avg_flankers=4;spanners=0;type=homozygous;TE=karma;TSD=TGATATTTGCAA
9        ID=Chr3.676509.spanners;avg_flankers=3.5;spanners=0;type=homozygous;TE=karma;TSD=GGGACTTAAAAG
10         ID=Chr3.775889.spanners;avg_flankers=5;spanners=0;type=homozygous;TE=karma;TSD=CCCAGACCTGTG
11        ID=Chr3.1180085.spanners;avg_flankers=2;spanners=0;type=homozygous;TE=karma;TSD=AAGGCCGTGGCG
12        ID=Chr3.1180315.spanners;avg_flankers=5;spanners=0;type=homozygous;TE=karma;TSD=GTAATCAGCTAG
13        ID=Chr3.1217991.spanners;avg_flankers=3;spanners=0;type=homozygous;TE=karma;TSD=CCCCAGAATATT
14      ID=Chr3.1222327.spanners;avg_flankers=2.5;spanners=0;type=homozygous;TE=karma;TSD=TCAGTCATTGGT
15        ID=Chr3.1472727.spanners;avg_flankers=3;spanners=0;type=homozygous;TE=karma;TSD=CATCGAGAGGGT
16 ID=Chr3.1574123.spanners;avg_flankers=5;spanners=0;type=homozygous;TE=karma/Retro1;TSD=TTATTTTCATTA
17        ID=Chr3.1581478.spanners;avg_flankers=3;spanners=0;type=homozygous;TE=karma;TSD=CTCGGTCGGTAT
18        ID=Chr3.1680014.spanners;avg_flankers=4;spanners=0;type=homozygous;TE=karma;TSD=TGGTCCCGGTGG
19      ID=Chr3.1696147.spanners;avg_flankers=2.5;spanners=0;type=homozygous;TE=karma;TSD=GGGCAGAAATGG
20      ID=Chr3.1804314.spanners;avg_flankers=4.5;spanners=0;type=homozygous;TE=karma;TSD=GCATGGGGTTCT
21        ID=Chr3.1835591.spanners;avg_flankers=5;spanners=0;type=homozygous;TE=karma;TSD=TAGGTCTGCCAT
22        ID=Chr3.1880410.spanners;avg_flankers=3;spanners=0;type=homozygous;TE=karma;TSD=GGGCCAGCCGCA
23        ID=Chr3.1892945.spanners;avg_flankers=2;spanners=0;type=homozygous;TE=karma;TSD=GATGTGGTGCAG

Overlapping with MuWU results:

   Chr.x Start.x   End.x Chr.y Start.y   End.y Sample StartReads EndReads
1   Chr3  109857  109868  Chr3  109857  109868  karma          2        1
2   Chr3  109857  109868  Chr3  109858  109869  karma          1        1
3   Chr3  155951  155962  Chr3  155951  155962  karma          2        4
4   Chr3  190294  190305  Chr3  190294  190305  karma          4        4
5   Chr3  355582  355593  Chr3  355582  355593  karma          4        3
6   Chr3  400634  400645  Chr3  400634  400645  karma          1        1
7   Chr3  481214  481225  Chr3  481214  481225  karma          2        1
8   Chr3  513469  513480  Chr3  513469  513480  karma          2        2
9   Chr3  667120  667131  Chr3  667120  667131  karma          2        1
10  Chr3  676498  676509  Chr3  676498  676509  karma          4        3
11  Chr3  775878  775889  Chr3  775878  775889  karma          1        3
12  Chr3 1180074 1180085  Chr3 1180074 1180085  karma          4        2
13  Chr3 1180304 1180315  Chr3 1180304 1180315  karma          2        1
14  Chr3 1217980 1217991  Chr3 1217980 1217991  karma          2        1
15  Chr3 1574112 1574123  Chr3 1574112 1574123  karma          4        1
16  Chr3 1581467 1581478  Chr3 1581467 1581478  karma          3        1
17  Chr3 1680003 1680014  Chr3 1680003 1680014  karma          4        3
18  Chr3 1696136 1696147  Chr3 1696136 1696147  karma          2        2
19  Chr3 1804303 1804314  Chr3 1804303 1804314  karma          3        4
20  Chr3 1835580 1835591  Chr3 1835580 1835591  karma          6        2
21  Chr3 1880399 1880410  Chr3 1880399 1880410  karma          4        2
22  Chr3 1892934 1892945  Chr3 1892934 1892945  karma          1        3

Consider, that MuWU would normally be used after a targeted PCR approach.
We nevertheless identify 21 of the 23 karma insertions as identified by RelocaTE2
```