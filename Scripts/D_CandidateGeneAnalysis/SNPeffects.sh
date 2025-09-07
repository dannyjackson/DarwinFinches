# Code used to investigate specific mutations in genes, looking especially for coding or CpG mutations

cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/genelist/gene_analyses/indv_genes/cra
FST="/xdisk/mcnew/finches/dannyjackson/finches/analyses/fst/crapre_crapost/1/crapre_crapost.1.fst"

# ANGPT1
## Cluster of high FST snps between exon 1 and 2, no high FST snps in an exon
## NC_044572.1:133469688-133623609	
### chr 2
awk -F'\t' 'BEGIN {OFS="\t"} NR==1 || ($2 == 2 && $3 >= 133469688 && $3 <= 133623609)' $FST > ANGPT1.fst.tsv


sort -k5 ANGPT1.fst.tsv | tail -n 20
region  chr     midPos  Nsites  fst
(1097950,1097950)(133552154,133552154)(133552154,133552155)     2       133552154       2       0.301759 # not a CpG or exon...
(1097962,1097962)(133554495,133554495)(133554495,133554496)     2       133554495       2       0.303404 # not a CpG or exon...
(1098605,1098605)(133611559,133611559)(133611559,133611560)     2       133611559       2       0.309636 # not a CpG or exon...
(1097699,1097699)(133525845,133525845)(133525845,133525846)     2       133525845       2       0.327387 # not a CpG or exon...
(1097657,1097657)(133520546,133520546)(133520546,133520547)     2       133520546       2       0.329191 # not a CpG or exon...
(1098628,1098628)(133616393,133616393)(133616393,133616394)     2       133616393       2       0.335755 # not a CpG or exon...
(1097868,1097868)(133540708,133540708)(133540708,133540709)     2       133540708       2       0.336797 # not a CpG or exon...
(1097701,1097701)(133526063,133526063)(133526063,133526064)     2       133526063       2       0.342401 # not a CpG or exon but very close to exon 5
(1098695,1098695)(133622792,133622792)(133622792,133622793)     2       133622792       2       0.385923 # not a CpG or exon...
(1098252,1098252)(133585675,133585675)(133585675,133585676)     2       133585675       2       0.386016 # not a CpG or exon...
(1098690,1098690)(133622521,133622521)(133622521,133622522)     2       133622521       2       0.391499 # not a CpG or exon...
(1098548,1098548)(133607177,133607177)(133607177,133607178)     2       133607177       2       0.422368 # not a CpG or exon...
(1098693,1098693)(133622572,133622572)(133622572,133622573)     2       133622572       2       0.430909 # not a CpG or exon...
(1098287,1098287)(133590204,133590204)(133590204,133590205)     2       133590204       2       0.461506 # not a CpG or exon...
(1098495,1098495)(133603632,133603632)(133603632,133603633)     2       133603632       2       0.484333 # not a CpG or exon...
(1098288,1098288)(133591034,133591034)(133591034,133591035)     2       133591034       2       0.501876 # not a CpG or exon...
(1098339,1098339)(133594381,133594381)(133594381,133594382)     2       133594381       2       0.514035 # not a CpG or exon...
(1098417,1098417)(133600392,133600392)(133600392,133600393)     2       133600392       2       0.523367 # not a CpG or exon...
(1098315,1098315)(133593164,133593164)(133593164,133593165)     2       133593164       2       0.577674 # not a CpG or exon... #



~/programs/angsd_elgato/angsd/angsd \
  -GL 1 \
  -doGeno 2 \
  -doPost 1 \
  -doCounts 1 \
  -doMajorMinor 3 \
  -doMaf 1 \
  -doDepth 1 \
  -setMinDepthInd 4 \
  -minQ 30 \
  -minMapQ 30 \
  -r NC_044572.1:133552154 \
  -sites /xdisk/mcnew/finches/dannyjackson/finches/referencelists/allsnps.sites_headless.mafs \
  -bam /xdisk/mcnew/finches/dannyjackson/finches/referencelists/crapostbams.txt \
  -anc /xdisk/mcnew/finches/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa \
  -out single_site

zcat single_site.mafs.gz single_site.geno.gz 

#### Exon 9
133469688-133471933
awk -F'\t' 'BEGIN {OFS="\t"} NR==1 || ($2 == 2 && $3 >= 133469688 && $3 <= 133471933)' ANGPT1.fst.tsv | sort -k5 -g # max FST = 0.054519
#### Exon 8
133480883-133481013
awk -F'\t' 'BEGIN {OFS="\t"} NR==1 || ($2 == 2 && $3 >= 133480883 && $3 <= 133481013)' ANGPT1.fst.tsv | sort -k5 -g # none
#### Exon 7
133506126-133506292
awk -F'\t' 'BEGIN {OFS="\t"} NR==1 || ($2 == 2 && $3 >= 133506126 && $3 <= 133506292)' ANGPT1.fst.tsv | sort -k5 -g # none
#### Exon 6
133520578-133520679
awk -F'\t' 'BEGIN {OFS="\t"} NR==1 || ($2 == 2 && $3 >= 133520578 && $3 <= 133520679)' ANGPT1.fst.tsv | sort -k5 -g # none
#### Exon 5
133526072-133526199
awk -F'\t' 'BEGIN {OFS="\t"} NR==1 || ($2 == 2 && $3 >= 133526072 && $3 <= 133526199)' ANGPT1.fst.tsv | sort -k5 -g # max FST = 0.015850
#### Exon 4
133543331-133543563
awk -F'\t' 'BEGIN {OFS="\t"} NR==1 || ($2 == 2 && $3 >= 133543331 && $3 <= 133543563)' ANGPT1.fst.tsv | sort -k5 -g # max FST = 0.052670
#### Exon 3
133552262-133552383
awk -F'\t' 'BEGIN {OFS="\t"} NR==1 || ($2 == 2 && $3 >= 133552262 && $3 <= 133552383)' ANGPT1.fst.tsv | sort -k5 -g # max FST = 0.011172
#### Exon 2
133553271-133553426
awk -F'\t' 'BEGIN {OFS="\t"} NR==1 || ($2 == 2 && $3 >= 133553271 && $3 <= 133553426)' ANGPT1.fst.tsv | sort -k5 -g # max FST = 0.018203
#### Exon 1
133622956-133623609
awk -F'\t' 'BEGIN {OFS="\t"} NR==1 || ($2 == 2 && $3 >= 133622956 && $3 <= 133623609)' ANGPT1.fst.tsv | sort -k5 -g # max FST = 0.017017

# HPSE
## NC_044574.1:59804605-59815556	
awk -F'\t' 'BEGIN {OFS="\t"} NR==1 || ($2 == 4 && $3 >= 59804605 && $3 <= 59815556)' $FST > HPSE.fst.tsv
NC_044574.1 59815472 C T
NC_044574.1 59815473 A T
NC_044574.1 59815483 A C
sort -k5 HPSE.fst.tsv | tail
# second highest snp, 59815456, is in the 3' noncoding but transcribed region of exon 12 (3'UTR)
59815450-59815460
# get info from the vcf
~/programs/angsd/angsd \
  -GL 1 \
  -doGeno 2 \
  -doPost 1 \
  -doCounts 1 \
  -doMajorMinor 3 \
  -doMaf 1 \
  -doDepth 1 \
  -setMinDepthInd 4 \
  -minQ 30 \
  -minMapQ 30 \
  -r NC_044574.1:59815450-59815460 \
  -sites /xdisk/mcnew/finches/dannyjackson/finches/referencelists/allsnps.sites_headless.mafs \
  -bam /xdisk/mcnew/finches/dannyjackson/finches/referencelists/craprebams.txt \
  -anc /xdisk/mcnew/finches/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa \
  -out single_site_59815456

0 = homozygous major
1 = heterozygous
2 = homozygous minor
3 = missing
-1 = missing

# an increase in major sites, loss of minor allele in response to selection at this site.
zcat single_site_59815456.mafs.gz single_site_59815456.geno.gz
# POST
chromo  position        major   minor   anc     knownEM nInd
NC_044574.1     59815456        G       C       G       0.127546        8
NC_044574.1     59815456        1       1       0       0       -1      0       0       0       0
# PRE
chromo  position        major   minor   anc     knownEM nInd
NC_044574.1     59815456        G       C       G       0.497128        4
NC_044574.1     59815456        2       0       1       -1      1

#### Exon 1
59804605-59804943
awk -F'\t' 'BEGIN {OFS="\t"} NR==1 || ($2 == 4 && $3 >= 59804605 && $3 <= 59804943)' HPSE.fst.tsv | sort -k5 -g # max FST = 0.099870
#### Exon 2
59806429-59806574
awk -F'\t' 'BEGIN {OFS="\t"} NR==1 || ($2 == 4 && $3 >= 59806429 && $3 <= 59806574)' HPSE.fst.tsv | sort -k5 -g # max FST = 0.022266
#### Exon 3
59807574-59807699
awk -F'\t' 'BEGIN {OFS="\t"} NR==1 || ($2 == 4 && $3 >= 59807574 && $3 <= 59807699)' HPSE.fst.tsv | sort -k5 -g # max FST = 0.057582
#### Exon 4
59808056-59808229 
awk -F'\t' 'BEGIN {OFS="\t"} NR==1 || ($2 == 4 && $3 >= 59808056 && $3 <= 59808229)' HPSE.fst.tsv | sort -k5 -g # none
#### Exon 5
59808673-59808844
awk -F'\t' 'BEGIN {OFS="\t"} NR==1 || ($2 == 4 && $3 >= 59808673 && $3 <= 59808844)' HPSE.fst.tsv | sort -k5 -g # max FST = 0.012084
#### Exon 6
59809360-59809407
awk -F'\t' 'BEGIN {OFS="\t"} NR==1 || ($2 == 4 && $3 >= 59809360 && $3 <= 59809407)' HPSE.fst.tsv | sort -k5 -g # max FST = 0.014061
#### Exon 7
59810789-59810882
awk -F'\t' 'BEGIN {OFS="\t"} NR==1 || ($2 == 4 && $3 >= 59810789 && $3 <= 59810882)' HPSE.fst.tsv | sort -k5 -g # none
#### Exon 8
59811510-59811616
awk -F'\t' 'BEGIN {OFS="\t"} NR==1 || ($2 == 4 && $3 >= 59811510 && $3 <= 59811616)' HPSE.fst.tsv | sort -k5 -g # none
#### Exon 9
59812128-59812242
awk -F'\t' 'BEGIN {OFS="\t"} NR==1 || ($2 == 4 && $3 >= 59812128 && $3 <= 59812242)' HPSE.fst.tsv | sort -k5 -g # none
#### Exon 10
59813879-59813997
awk -F'\t' 'BEGIN {OFS="\t"} NR==1 || ($2 == 4 && $3 >= 59813879 && $3 <= 59813997)' HPSE.fst.tsv | sort -k5 -g # none
#### Exon 11
59814131-59814277
awk -F'\t' 'BEGIN {OFS="\t"} NR==1 || ($2 == 4 && $3 >= 59814131 && $3 <= 59814277)' HPSE.fst.tsv | sort -k5 -g # max FST = 0.025356
#### Exon 12
59815031-59815556
awk -F'\t' 'BEGIN {OFS="\t"} NR==1 || ($2 == 4 && $3 >= 59815031 && $3 <= 59815556)' HPSE.fst.tsv | sort -k5 -g # max FST = 0.327240
# (499664,499664)(59815456,59815456)(59815456,59815457)   4       59815456        2       0.327240

# ITGA2B
## There is a block of high FST in the region after ITGA2B
## Gene boundaries
## NC_044597.1:1174211-1182568	
### chr 27
awk -F'\t' 'BEGIN {OFS="\t"} NR==1 || ($2 == 27 && $3 >= 1174211 && $3 <= 1182568)' $FST > ITGA2B.fst.tsv
sort -k5 ITGA2B.fst.tsv | tail
## Window boundaries
## NC_044597.1:1100000-1200000	
### chr 27
awk -F'\t' 'BEGIN {OFS="\t"} NR==1 || ($2 == 27 && $3 >= 1100000 && $3 <= 1200000)' $FST > ITGA2B.fst.window.tsv
sort -k5 ITGA2B.fst.window.tsv | tail

(7791,7791)(1185472,1185472)(1185472,1185473)   27      1185472 2       0.237483
(7770,7770)(1182713,1182713)(1182713,1182714)   27      1182713 2       0.258380
(7850,7850)(1195848,1195848)(1195848,1195849)   27      1195848 2       0.265959
(7824,7824)(1191787,1191787)(1191787,1191788)   27      1191787 2       0.282932
(7793,7793)(1185687,1185687)(1185687,1185688)   27      1185687 2       0.292889
(7871,7871)(1199065,1199065)(1199065,1199066)   27      1199065 2       0.300759
(7785,7785)(1183852,1183852)(1183852,1183853)   27      1183852 2       0.352795
(7842,7842)(1194683,1194683)(1194683,1194684)   27      1194683 2       0.359753
(7860,7860)(1197354,1197354)(1197354,1197355)   27      1197354 2       0.372513
(7805,7805)(1188499,1188499)(1188499,1188500)   27      1188499 2       0.375165
(7790,7790)(1185310,1185310)(1185310,1185311)   27      1185310 2       0.375760
(7825,7825)(1191916,1191916)(1191916,1191917)   27      1191916 2       0.380167
(7813,7813)(1189733,1189733)(1189733,1189734)   27      1189733 2       0.382977
(7792,7792)(1185599,1185599)(1185599,1185600)   27      1185599 2       0.387932
(7818,7818)(1190750,1190750)(1190750,1190751)   27      1190750 2       0.390260
(7782,7782)(1183583,1183583)(1183583,1183584)   27      1183583 2       0.401805
(7821,7821)(1191212,1191212)(1191212,1191213)   27      1191212 2       0.408694
(7779,7779)(1183519,1183519)(1183519,1183520)   27      1183519 2       0.432111
(7788,7788)(1184438,1184438)(1184438,1184439)   27      1184438 2       0.467511
(7833,7833)(1193890,1193890)(1193890,1193891)   27      1193890 2       0.501025
(7851,7851)(1195930,1195930)(1195930,1195931)   27      1195930 2       0.530085
(7780,7780)(1183526,1183526)(1183526,1183527)   27      1183526 2       0.548984
(7775,7775)(1183287,1183287)(1183287,1183288)   27      1183287 2       0.580132
(7829,7829)(1192746,1192746)(1192746,1192747)   27      1192746 2       0.592209
(7849,7849)(1195693,1195693)(1195693,1195694)   27      1195693 2       0.626504
(7834,7834)(1193921,1193921)(1193921,1193922)   27      1193921 2       0.633616
(7832,7832)(1193786,1193786)(1193786,1193787)   27      1193786 2       0.641068
(7787,7787)(1184274,1184274)(1184274,1184275)   27      1184274 2       0.661058
(7774,7774)(1183262,1183262)(1183262,1183263)   27      1183262 2       0.665178



# PODXL
## NC_044586.1:67365108-67411859
### chr 1A
grep '1A' $FST | awk -F'\t' 'BEGIN {OFS="\t"} NR==1 || ($3 >= 67365108 && $3 <= 67411859)' > PODXL.fst.tsv
sort -k5 PODXL.fst.tsv | tail -n 30
# missense mutation: 67381351

(577462,577462)(67374315,67374315)(67374315,67374316)   1A      67374315        2       0.401530
(577439,577439)(67372984,67372984)(67372984,67372985)   1A      67372984        2       0.404963
(577415,577415)(67370799,67370799)(67370799,67370800)   1A      67370799        2       0.407155
(577476,577476)(67374999,67374999)(67374999,67375000)   1A      67374999        2       0.436989
(577421,577421)(67371273,67371273)(67371273,67371274)   1A      67371273        2       0.439638 # within exon 7
(577562,577562)(67380375,67380375)(67380375,67380376)   1A      67380375        2       0.441156
(577592,577592)(67381955,67381955)(67381955,67381956)   1A      67381955        2       0.441159 # within exon 2
(577554,577554)(67379695,67379695)(67379695,67379696)   1A      67379695        2       0.446048
(577585,577585)(67381351,67381351)(67381351,67381352)   1A      67381351        2       0.448900 # within exon 2, near 67381354
(577436,577436)(67372732,67372732)(67372732,67372733)   1A      67372732        2       0.457514
(577557,577557)(67379978,67379978)(67379978,67379979)   1A      67379978        2       0.460745
(577558,577558)(67379982,67379982)(67379982,67379983)   1A      67379982        2       0.477208
(577555,577555)(67379743,67379743)(67379743,67379744)   1A      67379743        2       0.508188
(577544,577544)(67379458,67379458)(67379458,67379459)   1A      67379458        2       0.511110
(577474,577474)(67374839,67374839)(67374839,67374840)   1A      67374839        2       0.512435
(577586,577586)(67381354,67381354)(67381354,67381355)   1A      67381354        2       0.515044 # within exon 2
(577473,577473)(67374803,67374803)(67374803,67374804)   1A      67374803        2       0.526073
(577583,577583)(67381237,67381237)(67381237,67381238)   1A      67381237        2       0.528400
(577584,577584)(67381246,67381246)(67381246,67381247)   1A      67381246        2       0.532308
(577468,577468)(67374527,67374527)(67374527,67374528)   1A      67374527        2       0.536468
(577573,577573)(67380967,67380967)(67380967,67380968)   1A      67380967        2       0.566537
(577469,577469)(67374587,67374587)(67374587,67374588)   1A      67374587        2       0.566669
(577576,577576)(67381056,67381056)(67381056,67381057)   1A      67381056        2       0.595505
(577463,577463)(67374317,67374317)(67374317,67374318)   1A      67374317        2       0.607642
(577572,577572)(67380966,67380966)(67380966,67380967)   1A      67380966        2       0.621728
(577574,577574)(67380968,67380968)(67380968,67380969)   1A      67380968        2       0.624351
(577453,577453)(67373799,67373799)(67373799,67373800)   1A      67373799        2       0.635237
(577461,577461)(67374192,67374192)(67374192,67374193)   1A      67374192        2       0.635423
(577587,577587)(67381476,67381476)(67381476,67381477)   1A      67381476        2       0.667502 # within exon 2
(577458,577458)(67374004,67374004)(67374004,67374005)   1A      67374004        2       0.667848



~/programs/angsd/angsd \
  -GL 1 \
  -doGeno 2 \
  -doPost 1 \
  -doCounts 1 \
  -doMajorMinor 3 \
  -doMaf 1 \
  -doDepth 1 \
  -setMinDepthInd 4 \
  -minQ 30 \
  -minMapQ 30 \
  -r NC_044586.1:67381351 \
  -sites /xdisk/mcnew/finches/dannyjackson/finches/referencelists/allsnps.sites_headless.mafs \
  -bam /xdisk/mcnew/finches/dannyjackson/finches/referencelists/crapostbams.txt \
  -anc /xdisk/mcnew/finches/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa \
  -out single_site_post

zcat single_site_post.mafs.gz

~/programs/angsd/angsd \
  -GL 1 \
  -doGeno 2 \
  -doPost 1 \
  -doCounts 1 \
  -doMajorMinor 3 \
  -doMaf 1 \
  -doDepth 1 \
  -setMinDepthInd 4 \
  -minQ 30 \
  -minMapQ 30 \
  -r NC_044586.1:67381351 \
  -sites /xdisk/mcnew/finches/dannyjackson/finches/referencelists/allsnps.sites_headless.mafs \
  -bam /xdisk/mcnew/finches/dannyjackson/finches/referencelists/craprebams.txt \
  -anc /xdisk/mcnew/finches/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa \
  -out single_site_pre

zcat single_site_post.mafs.gz

zcat single_site_post.geno.gz 
zcat single_site_pre.geno.gz 



# look at TEAD4 and YAP1
# YAP1
## NC_044571.1:75789510-75877796
### chr 1
awk -F'\t' 'BEGIN {OFS="\t"} NR==1 || ($2 == 1 && $3 >= 75789510 && $3 <= 75877796)' $FST > YAP1.fst.tsv
sort -k5 YAP1.fst.tsv | tail

(667208,667208)(75865819,75865819)(75865819,75865820)   1       75865819        2       0.285528 # not coding, not cpg
(667191,667191)(75863909,75863909)(75863909,75863910)   1       75863909        2       0.287913 # not coding, not cpg
(667243,667243)(75872061,75872061)(75872061,75872062)   1       75872061        2       0.295235 # not coding, not cpg
(667216,667216)(75867659,75867659)(75867659,75867660)   1       75867659        2       0.298887 # intronic between exon 2 and 3; CpG!
(667207,667207)(75865564,75865564)(75865564,75865565)   1       75865564        2       0.299544 # intronic between exon 2 and 3; CpG!
(667242,667242)(75871868,75871868)(75871868,75871869)   1       75871868        2       0.300368 # not coding, not cpg
(667203,667203)(75865086,75865086)(75865086,75865087)   1       75865086        2       0.306153 # not coding, not cpg
(667210,667210)(75866111,75866111)(75866111,75866112)   1       75866111        2       0.367343 # not coding, not cpg
(667192,667192)(75864035,75864035)(75864035,75864036)   1       75864035        2       0.393496 # not coding, not cpg

chromo  position        major   minor   anc     knownEM nInd
NC_044571.1     75865564        A       G       A       0.054887        9
NC_044571.1     75865564        0       0       0       0       0       0       0       0       1
NC_044571.1     75865564        1       0       1       2       0

chromo  position        major   minor   anc     knownEM nInd
NC_044571.1     75867659        A       G       A       0.056462        9
NC_044571.1     75867659        0       0       0       0       0       0       0       0       1
NC_044571.1     75867659        1       0       1       2       0

~/programs/angsd/angsd \
  -GL 1 \
  -doGeno 2 \
  -doPost 1 \
  -doCounts 1 \
  -doMajorMinor 3 \
  -doMaf 1 \
  -doDepth 1 \
  -setMinDepthInd 4 \
  -minQ 30 \
  -minMapQ 30 \
  -r NC_044571.1:75865819 \
  -sites /xdisk/mcnew/finches/dannyjackson/finches/referencelists/allsnps.sites_headless.mafs \
  -bam /xdisk/mcnew/finches/dannyjackson/finches/referencelists/crapostbams.txt \
  -anc /xdisk/mcnew/finches/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa \
  -out single_site_post

zcat single_site_post.mafs.gz

~/programs/angsd/angsd \
  -GL 1 \
  -doGeno 2 \
  -doPost 1 \
  -doCounts 1 \
  -doMajorMinor 3 \
  -doMaf 1 \
  -doDepth 1 \
  -setMinDepthInd 4 \
  -minQ 30 \
  -minMapQ 30 \
  -r NC_044571.1:75865819 \
  -sites /xdisk/mcnew/finches/dannyjackson/finches/referencelists/allsnps.sites_headless.mafs \
  -bam /xdisk/mcnew/finches/dannyjackson/finches/referencelists/craprebams.txt \
  -anc /xdisk/mcnew/finches/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa \
  -out single_site_pre

zcat single_site_post.mafs.gz

zcat single_site_post.geno.gz 
zcat single_site_pre.geno.gz 


# TEAD4
## NC_044571.1:87095605-87147397
### chr 1
cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/genelist/gene_analyses/indv_genes/par
FST="/xdisk/mcnew/finches/dannyjackson/finches/analyses/fst/parpre_parpost/1/parpre_parpost.1.fst"

awk -F'\t' 'BEGIN {OFS="\t"} NR==1 || ($2 == 1 && $3 >= 87095605 && $3 <= 87147397)' $FST > TEAD4.fst.tsv
sort -k5 TEAD4.fst.tsv | tail

(759621,759621)(87127610,87127610)(87127610,87127611)   1       87127610        2       0.203289
(759482,759482)(87108547,87108547)(87108547,87108548)   1       87108547        2       0.203629

## Window boundaries
## NC_044571.1:87100000-87200000	
awk -F'\t' 'BEGIN {OFS="\t"} NR==1 || ($2 == 1 && $3 >= 87050000 && $3 <= 87200000)' $FST > TEAD4.window.fst.tsv
sort -k5 TEAD4.window.fst.tsv | tail

## hmm, all of these are within TSPAN9
87170612 - 87173371
NC_044571.1:87153063-87322660
(759924,759924)(87171219,87171219)(87171219,87171220)   1       87171219        2       0.612386
(759923,759923)(87171218,87171218)(87171218,87171219)   1       87171218        2       0.612393
(759913,759913)(87170797,87170797)(87170797,87170798)   1       87170797        2       0.612440
(759953,759953)(87173371,87173371)(87173371,87173372)   1       87173371        2       0.612504
(759914,759914)(87170816,87170816)(87170816,87170817)   1       87170816        2       0.612641
(759916,759916)(87170869,87170869)(87170869,87170870)   1       87170869        2       0.612678
(759915,759915)(87170845,87170845)(87170845,87170846)   1       87170845        2       0.613120
(759940,759940)(87171968,87171968)(87171968,87171969)   1       87171968        2       0.618860
(759909,759909)(87170612,87170612)(87170612,87170613)   1       87170612        2       0.653203



# BMPER, NUDT4, and TENM3 in all three species
## CRA
cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/genelist/gene_analyses/indv_genes/cra
FST="/xdisk/mcnew/finches/dannyjackson/finches/analyses/fst/crapre_crapost/1/crapre_crapost.1.fst"
### BMPER
#### NC_044572.1:46478775-46625955
##### CHR 2
awk -F'\t' 'BEGIN {OFS="\t"} NR==1 || ($2 == 2 && $3 >= 46478775 && $3 <= 46625955)' $FST > BMPER.fst.tsv
sort -k5 BMPER.fst.tsv | tail

(400910,400910)(46614822,46614822)(46614822,46614823)   2       46614822        2       0.235733 # not coding, not CpG
(400456,400456)(46545222,46545222)(46545222,46545223)   2       46545222        2       0.239838 # not coding, not CpG
(400163,400163)(46511099,46511099)(46511099,46511100)   2       46511099        2       0.246331 # CpG site!!! Noncoding, between exon 3 and 4
(400150,400150)(46509534,46509534)(46509534,46509535)   2       46509534        2       0.286284 # not coding, not CpG
(400151,400151)(46509653,46509653)(46509653,46509654)   2       46509653        2       0.296369 # not coding, not CpG
(400154,400154)(46509974,46509974)(46509974,46509975)   2       46509974        2       0.301060 # not coding, not CpG
(400365,400365)(46531886,46531886)(46531886,46531887)   2       46531886        2       0.309886 # not coding, not CpG; long string of 5'-TTTTT-3'
(400939,400939)(46621349,46621349)(46621349,46621350)   2       46621349        2       0.327547 # not coding, not CpG; long string of 5'-TTTTT-3'
(400158,400158)(46510683,46510683)(46510683,46510684)   2       46510683        2       0.396960 # not coding, not CpG


chromo  position        major   minor   anc     knownEM nInd
NC_044572.1     46511099        C       T       T       0.109415        7
NC_044572.1     46511099        0       0       1       -1      -1      0       0       0       0
NC_044572.1     46511099        1       1       1       1       1
# 5 were heterozygous in pre, now 1 heterozygous and 6 homozygous for major; selection for CpG site.

~/programs/angsd/angsd \
  -GL 1 \
  -doGeno 2 \
  -doPost 1 \
  -doCounts 1 \
  -doMajorMinor 3 \
  -doMaf 1 \
  -doDepth 1 \
  -setMinDepthInd 4 \
  -minQ 30 \
  -minMapQ 30 \
  -r NC_044572.1:46511099 \
  -sites /xdisk/mcnew/finches/dannyjackson/finches/referencelists/allsnps.sites_headless.mafs \
  -bam /xdisk/mcnew/finches/dannyjackson/finches/referencelists/crapostbams.txt \
  -anc /xdisk/mcnew/finches/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa \
  -out single_site_post

zcat single_site_post.mafs.gz

~/programs/angsd/angsd \
  -GL 1 \
  -doGeno 2 \
  -doPost 1 \
  -doCounts 1 \
  -doMajorMinor 3 \
  -doMaf 1 \
  -doDepth 1 \
  -setMinDepthInd 4 \
  -minQ 30 \
  -minMapQ 30 \
  -r NC_044572.1:46511099 \
  -sites /xdisk/mcnew/finches/dannyjackson/finches/referencelists/allsnps.sites_headless.mafs \
  -bam /xdisk/mcnew/finches/dannyjackson/finches/referencelists/craprebams.txt \
  -anc /xdisk/mcnew/finches/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa \
  -out single_site_pre

zcat single_site_post.mafs.gz

zcat single_site_post.geno.gz 
zcat single_site_pre.geno.gz 


## FOR
cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/genelist/gene_analyses/indv_genes/for
FST="/xdisk/mcnew/finches/dannyjackson/finches/analyses/fst/forpre_forpost/1/forpre_forpost.1.fst"
### BMPER
#### NC_044572.1:46478775-46625955
##### CHR 2
awk -F'\t' 'BEGIN {OFS="\t"} NR==1 || ($2 == 2 && $3 >= 46478775 && $3 <= 46625955)' $FST > BMPER.fst.tsv
sort -k5 BMPER.fst.tsv | tail -n 30

(400749,400749)(46574549,46574549)(46574549,46574550)   2       46574549        2       0.301754 # not coding, not CpG
(400810,400810)(46580828,46580828)(46580828,46580829)   2       46580828        2       0.303085 # CpG site between exon 13 and 14 (just after 13)
(400849,400849)(46584598,46584598)(46584598,46584599)   2       46584598        2       0.327521 # CpG site between exon 9 and 10
(400863,400863)(46587306,46587306)(46587306,46587307)   2       46587306        2       0.328390 # not coding, not CpG
(400813,400813)(46580935,46580935)(46580935,46580936)   2       46580935        2       0.331082 # CpG site between exon 13 and 14 (just after 13)
(400621,400621)(46553585,46553585)(46553585,46553586)   2       46553585        2       0.343005 # not coding, not CpG
(400819,400819)(46581259,46581259)(46581259,46581260)   2       46581259        2       0.345283 # not coding, not CpG
(400840,400840)(46583776,46583776)(46583776,46583777)   2       46583776        2       0.350549 # not coding, not CpG
(400850,400850)(46585031,46585031)(46585031,46585032)   2       46585031        2       0.355726 # not coding, not CpG
(400861,400861)(46587062,46587062)(46587062,46587063)   2       46587062        2       0.359803 # not coding, not CpG
(400924,400924)(46602867,46602867)(46602867,46602868)   2       46602867        2       0.376037 # CpG site between exon 14 and 15
(400783,400783)(46576940,46576940)(46576940,46576941)   2       46576940        2       0.402065 # CpG site between exon 12 and 13
(400770,400770)(46575714,46575714)(46575714,46575715)   2       46575714        2       0.406229 # not coding, not CpG
(400778,400778)(46576584,46576584)(46576584,46576585)   2       46576584        2       0.411312 # not coding, not CpG
(400782,400782)(46576742,46576742)(46576742,46576743)   2       46576742        2       0.411911 # not coding, not CpG
(400780,400780)(46576711,46576711)(46576711,46576712)   2       46576711        2       0.413786 # not coding, not CpG
(400775,400775)(46576067,46576067)(46576067,46576068)   2       46576067        2       0.425120 # not coding, not CpG
(400765,400765)(46575540,46575540)(46575540,46575541)   2       46575540        2       0.468533 # not coding, not CpG
(400772,400772)(46575836,46575836)(46575836,46575837)   2       46575836        2       0.473998 # not coding, not CpG


~/programs/angsd/angsd \
  -GL 1 \
  -doGeno 2 \
  -doPost 1 \
  -doCounts 1 \
  -doMajorMinor 3 \
  -doMaf 1 \
  -doDepth 1 \
  -setMinDepthInd 4 \
  -minQ 30 \
  -minMapQ 30 \
  -r NC_044572.1:46580828 \
  -sites /xdisk/mcnew/finches/dannyjackson/finches/referencelists/allsnps.sites_headless.mafs \
  -bam /xdisk/mcnew/finches/dannyjackson/finches/referencelists/forpostbams.txt \
  -anc /xdisk/mcnew/finches/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa \
  -out single_site_post

zcat single_site_post.mafs.gz

~/programs/angsd/angsd \
  -GL 1 \
  -doGeno 2 \
  -doPost 1 \
  -doCounts 1 \
  -doMajorMinor 3 \
  -doMaf 1 \
  -doDepth 1 \
  -setMinDepthInd 4 \
  -minQ 30 \
  -minMapQ 30 \
  -r NC_044572.1:46580828 \
  -sites /xdisk/mcnew/finches/dannyjackson/finches/referencelists/allsnps.sites_headless.mafs \
  -bam /xdisk/mcnew/finches/dannyjackson/finches/referencelists/forprebams.txt \
  -anc /xdisk/mcnew/finches/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa \
  -out single_site_pre

zcat single_site_post.mafs.gz

zcat single_site_post.geno.gz 
zcat single_site_pre.geno.gz 


## PAR
cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/genelist/gene_analyses/indv_genes/par
FST="/xdisk/mcnew/finches/dannyjackson/finches/analyses/fst/parpre_parpost/1/parpre_parpost.1.fst"
### BMPER
#### NC_044572.1:46478775-46625955
##### CHR 2
awk -F'\t' 'BEGIN {OFS="\t"} NR==1 || ($2 == 2 && $3 >= 46478775 && $3 <= 46625955)' $FST > BMPER.fst.tsv
sort -k5 BMPER.fst.tsv | tail

(400836,400836)(46578649,46578649)(46578649,46578650)   2       46578649        2       0.245141 # CpG site between exon 12 and 13
(400592,400592)(46547030,46547030)(46547030,46547031)   2       46547030        2       0.245156 # not coding, not CpG
(400858,400858)(46581243,46581243)(46581243,46581244)   2       46581243        2       0.254531 # not coding, not CpG
(400310,400310)(46513495,46513495)(46513495,46513496)   2       46513495        2       0.257733 # CpG site between exon 3 and 4, very close to 4 (46514490, 995 sites)
(400794,400794)(46574902,46574902)(46574902,46574903)   2       46574902        2       0.258011 # CpG, adjacent to noncoding nonCpG site two rows down, between exon 9 and 10
(400838,400838)(46578843,46578843)(46578843,46578844)   2       46578843        2       0.275800 # not coding, not CpG
(400793,400793)(46574901,46574901)(46574901,46574902)   2       46574901        2       0.303619 # not coding, not CpG
(400719,400719)(46564789,46564789)(46564789,46564790)   2       46564789        2       0.309934 # not coding, not CpG
(400852,400852)(46580840,46580840)(46580840,46580841)   2       46580840        2       0.330260 # not coding, not CpG


46574901
46574902

~/programs/angsd/angsd \
  -GL 1 \
  -doGeno 2 \
  -doPost 1 \
  -doCounts 1 \
  -doMajorMinor 3 \
  -doMaf 1 \
  -doDepth 1 \
  -setMinDepthInd 4 \
  -minQ 30 \
  -minMapQ 30 \
  -r NC_044572.1:46513495 \
  -sites /xdisk/mcnew/finches/dannyjackson/finches/referencelists/allsnps.sites_headless.mafs \
  -bam /xdisk/mcnew/finches/dannyjackson/finches/referencelists/parpostbams.txt \
  -anc /xdisk/mcnew/finches/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa \
  -out single_site_post


~/programs/angsd/angsd \
  -GL 1 \
  -doGeno 2 \
  -doPost 1 \
  -doCounts 1 \
  -doMajorMinor 3 \
  -doMaf 1 \
  -doDepth 1 \
  -setMinDepthInd 4 \
  -minQ 30 \
  -minMapQ 30 \
  -r NC_044572.1:46513495 \
  -sites /xdisk/mcnew/finches/dannyjackson/finches/referencelists/allsnps.sites_headless.mafs \
  -bam /xdisk/mcnew/finches/dannyjackson/finches/referencelists/parprebams.txt \
  -anc /xdisk/mcnew/finches/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa \
  -out single_site_pre

zcat single_site_post.mafs.gz

zcat single_site_post.geno.gz 
zcat single_site_pre.geno.gz 

# 
~/programs/angsd/angsd \
  -GL 1 \
  -doGeno 2 \
  -doPost 1 \
  -doCounts 1 \
  -doMajorMinor 3 \
  -doMaf 1 \
  -doDepth 1 \
  -setMinDepthInd 4 \
  -minQ 30 \
  -minMapQ 30 \
  -r NC_044572.1:46574901 \
  -sites /xdisk/mcnew/finches/dannyjackson/finches/referencelists/allsnps.sites_headless.mafs \
  -bam /xdisk/mcnew/finches/dannyjackson/finches/referencelists/parpostbams.txt \
  -anc /xdisk/mcnew/finches/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa \
  -out single_site_post_901


~/programs/angsd/angsd \
  -GL 1 \
  -doGeno 2 \
  -doPost 1 \
  -doCounts 1 \
  -doMajorMinor 3 \
  -doMaf 1 \
  -doDepth 1 \
  -setMinDepthInd 4 \
  -minQ 30 \
  -minMapQ 30 \
  -r NC_044572.1:46574901 \
  -sites /xdisk/mcnew/finches/dannyjackson/finches/referencelists/allsnps.sites_headless.mafs \
  -bam /xdisk/mcnew/finches/dannyjackson/finches/referencelists/parprebams.txt \
  -anc /xdisk/mcnew/finches/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa \
  -out single_site_pre_901

  chromo	position	major	minor	anc	knownEM	nInd							
NC_044572.1	46574901	C	A	C	0.856534	7							
NC_044572.1	46574901	1	2	2	2	1	2	2					
NC_044572.1	46574901	1	1	2	1	0	1	1	0	1	1	0	1
chromo	position	major	minor	anc	knownEM	nInd							
NC_044572.1	46574902	G	A	A	0.143376	7							
NC_044572.1	46574902	1	0	0	0	1	0	0					
NC_044572.1	46574902	0	1	0	1	2	1	1	2	1	1	2	1
													
	Post	C|A ; G|A	AA | GG	AA; GG	AA; GG	C|A ; G|A	AA; GG	AA;GG					
	Pre:	C|A ; G	C|A ; G|A	AA; GG	C|A ; G|A	CC; AA	C|A ; G|A	C|A ; G|A	CC; AA	C|A ; G|A	C|A ; G|A	CC; AA	C|A ; G|A

    # Kind of wild but it looks like this site is not a CpG because the mutants on 901 and 902 are on the same haplotype.

## CRA
cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/genelist/gene_analyses/indv_genes/cra
FST="/xdisk/mcnew/finches/dannyjackson/finches/analyses/fst/crapre_crapost/1/crapre_crapost.1.fst"
### NUDT4
#### NC_044586.1:28599179-28626256	
##### CHR 1A
awk -F'\t' 'BEGIN {OFS="\t"} NR==1 || ($2 == "1A" && $3 >= 28599179 && $3 <= 28626256)' $FST > NUDT4.fst.tsv
sort -k5 NUDT4.fst.tsv | tail

(234505,234505)(28622459,28622459)(28622459,28622460)   1A      28622459        2       0.096126 # not coding, not CpG
(234434,234434)(28611894,28611894)(28611894,28611895)   1A      28611894        2       0.112274 # not coding, not CpG
(234470,234470)(28618132,28618132)(28618132,28618133)   1A      28618132        2       0.112483 # not coding, not CpG
(234399,234399)(28603505,28603505)(28603505,28603506)   1A      28603505        2       0.155575 # not coding, not CpG
(234403,234403)(28604333,28604333)(28604333,28604334)   1A      28604333        2       0.156459 # not coding, not CpG
(234395,234395)(28603313,28603313)(28603313,28603314)   1A      28603313        2       0.244754 # CpG, between exon 4 and 5
(234370,234370)(28599273,28599273)(28599273,28599274)   1A      28599273        2       0.282672 # not coding, not CpG


~/programs/angsd/angsd \
  -GL 1 \
  -doGeno 2 \
  -doPost 1 \
  -doCounts 1 \
  -doMajorMinor 3 \
  -doMaf 1 \
  -doDepth 1 \
  -setMinDepthInd 4 \
  -minQ 30 \
  -minMapQ 30 \
  -r NC_044586.1:28622471 \
  -sites /xdisk/mcnew/finches/dannyjackson/finches/referencelists/allsnps.sites_headless.mafs \
  -bam /xdisk/mcnew/finches/dannyjackson/finches/referencelists/crapostbams.txt \
  -anc /xdisk/mcnew/finches/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa \
  -out single_site_post


~/programs/angsd/angsd \
  -GL 1 \
  -doGeno 2 \
  -doPost 1 \
  -doCounts 1 \
  -doMajorMinor 3 \
  -doMaf 1 \
  -doDepth 1 \
  -setMinDepthInd 4 \
  -minQ 30 \
  -minMapQ 30 \
  -r NC_044586.1:28622471 \
  -sites /xdisk/mcnew/finches/dannyjackson/finches/referencelists/allsnps.sites_headless.mafs \
  -bam /xdisk/mcnew/finches/dannyjackson/finches/referencelists/craprebams.txt \
  -anc /xdisk/mcnew/finches/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa \
  -out single_site_pre

zcat single_site_post.mafs.gz

zcat single_site_post.geno.gz 
zcat single_site_pre.geno.gz 

## FOR 
cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/genelist/gene_analyses/indv_genes/for
FST="/xdisk/mcnew/finches/dannyjackson/finches/analyses/fst/forpre_forpost/1/forpre_forpost.1.fst"

### NUDT4
#### NC_044586.1:28599179-28626256	
##### CHR 1A
awk -F'\t' 'BEGIN {OFS="\t"} NR==1 || ($2 == "1A" && $3 >= 28599179 && $3 <= 28626256)' $FST > NUDT4.fst.tsv
sort -k5 NUDT4.fst.tsv | tail -n 20

(234581,234581)(28621591,28621591)(28621591,28621592)   1A      28621591        2       0.205953 # CpG, between exon 1 and 2
(234583,234583)(28621731,28621731)(28621731,28621732)   1A      28621731        2       0.281011 # not coding, not CpG
(234595,234595)(28622638,28622638)(28622638,28622639)   1A      28622638        2       0.281034 # not coding, not CpG


~/programs/angsd/angsd \
  -GL 1 \
  -doGeno 2 \
  -doPost 1 \
  -doCounts 1 \
  -doMajorMinor 3 \
  -doMaf 1 \
  -doDepth 1 \
  -setMinDepthInd 4 \
  -minQ 30 \
  -minMapQ 30 \
  -r NC_044586.1:28621591 \
  -sites /xdisk/mcnew/finches/dannyjackson/finches/referencelists/allsnps.sites_headless.mafs \
  -bam /xdisk/mcnew/finches/dannyjackson/finches/referencelists/forpostbams.txt \
  -anc /xdisk/mcnew/finches/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa \
  -out single_site_post


~/programs/angsd/angsd \
  -GL 1 \
  -doGeno 2 \
  -doPost 1 \
  -doCounts 1 \
  -doMajorMinor 3 \
  -doMaf 1 \
  -doDepth 1 \
  -setMinDepthInd 4 \
  -minQ 30 \
  -minMapQ 30 \
  -r NC_044586.1:28621591 \
  -sites /xdisk/mcnew/finches/dannyjackson/finches/referencelists/allsnps.sites_headless.mafs \
  -bam /xdisk/mcnew/finches/dannyjackson/finches/referencelists/forprebams.txt \
  -anc /xdisk/mcnew/finches/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa \
  -out single_site_pre

zcat single_site_post.mafs.gz

zcat single_site_post.geno.gz 
zcat single_site_pre.geno.gz 

## PAR
cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/genelist/gene_analyses/indv_genes/par
FST="/xdisk/mcnew/finches/dannyjackson/finches/analyses/fst/parpre_parpost/1/parpre_parpost.1.fst"

### NUDT4
#### NC_044586.1:28599179-28626256	
##### CHR 1A
awk -F'\t' 'BEGIN {OFS="\t"} NR==1 || ($2 == "1A" && $3 >= 28599179 && $3 <= 28626256)' $FST > NUDT4.fst.tsv
sort -k5 NUDT4.fst.tsv | tail -n 20

(234502,234502)(28601561,28601561)(28601561,28601562)   1A      28601561        2       0.220293
(234497,234497)(28600324,28600324)(28600324,28600325)   1A      28600324        2       0.221644
(234557,234557)(28612172,28612172)(28612172,28612173)   1A      28612172        2       0.222042
(234516,234516)(28603347,28603347)(28603347,28603348)   1A      28603347        2       0.225766
(234521,234521)(28604044,28604044)(28604044,28604045)   1A      28604044        2       0.227664
(234531,234531)(28607310,28607310)(28607310,28607311)   1A      28607310        2       0.228620
(234569,234569)(28614924,28614924)(28614924,28614925)   1A      28614924        2       0.234236 # not coding, not CpG
(234526,234526)(28605185,28605185)(28605185,28605186)   1A      28605185        2       0.234581 # not coding, not CpG
(234496,234496)(28600154,28600154)(28600154,28600155)   1A      28600154        2       0.241479 # 5`UTR, CpG
(234506,234506)(28601802,28601802)(28601802,28601803)   1A      28601802        2       0.246576 # 5`UTR, not CpG
(234561,234561)(28614083,28614083)(28614083,28614084)   1A      28614083        2       0.246583 # not coding, not CpG
(234568,234568)(28614781,28614781)(28614781,28614782)   1A      28614781        2       0.249894 # CpG, between exon 1 and 2
(234528,234528)(28606861,28606861)(28606861,28606862)   1A      28606861        2       0.256904 # not coding, not CpG
(234491,234491)(28599511,28599511)(28599511,28599512)   1A      28599511        2       0.260273 # not coding, not CpG
(234534,234534)(28607563,28607563)(28607563,28607564)   1A      28607563        2       0.266425 # not coding, not CpG
(234538,234538)(28607834,28607834)(28607834,28607835)   1A      28607834        2       0.266934 # not coding, not CpG
(234525,234525)(28605052,28605052)(28605052,28605053)   1A      28605052        2       0.280145 # not coding, not CpG
(234505,234505)(28601647,28601647)(28601647,28601648)   1A      28601647        2       0.319059 # CpG within 5'UTR
(234560,234560)(28613783,28613783)(28613783,28613784)   1A      28613783        2       0.402059 # not coding, not CpG


~/programs/angsd/angsd \
  -GL 1 \
  -doGeno 2 \
  -doPost 1 \
  -doCounts 1 \
  -doMajorMinor 3 \
  -doMaf 1 \
  -doDepth 1 \
  -setMinDepthInd 4 \
  -minQ 30 \
  -minMapQ 30 \
  -r NC_044586.1:28600154 \
  -sites /xdisk/mcnew/finches/dannyjackson/finches/referencelists/allsnps.sites_headless.mafs \
  -bam /xdisk/mcnew/finches/dannyjackson/finches/referencelists/parpostbams.txt \
  -anc /xdisk/mcnew/finches/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa \
  -out single_site_post


~/programs/angsd/angsd \
  -GL 1 \
  -doGeno 2 \
  -doPost 1 \
  -doCounts 1 \
  -doMajorMinor 3 \
  -doMaf 1 \
  -doDepth 1 \
  -setMinDepthInd 4 \
  -minQ 30 \
  -minMapQ 30 \
  -r NC_044586.1:28600154 \
  -sites /xdisk/mcnew/finches/dannyjackson/finches/referencelists/allsnps.sites_headless.mafs \
  -bam /xdisk/mcnew/finches/dannyjackson/finches/referencelists/parprebams.txt \
  -anc /xdisk/mcnew/finches/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa \
  -out single_site_pre

zcat single_site_post.mafs.gz

zcat single_site_post.geno.gz 
zcat single_site_pre.geno.gz 



## CRA
cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/genelist/gene_analyses/indv_genes/cra
FST="/xdisk/mcnew/finches/dannyjackson/finches/analyses/fst/crapre_crapost/1/crapre_crapost.1.fst"
### TENM3
#### NC_044574.1:32971132-34255509
##### CHR 4
awk -F'\t' 'BEGIN {OFS="\t"} NR==1 || ($2 == 4 && $3 >= 32971132 && $3 <= 34255509)' $FST > TENM3.fst.tsv
sort -k5 TENM3.fst.tsv | tail -n 20

(295723,295723)(33979217,33979217)(33979217,33979218)   4       33979217        2       0.324344 # noncoding, not CpG
(294717,294717)(33852051,33852051)(33852051,33852052)   4       33852051        2       0.329928 # noncoding, not CpG
(294675,294675)(33848188,33848188)(33848188,33848189)   4       33848188        2       0.333416 # noncoding, not CpG
(294721,294721)(33852521,33852521)(33852521,33852522)   4       33852521        2       0.339110 # noncoding, not CpG
(288949,288949)(33030382,33030382)(33030382,33030383)   4       33030382        2       0.350997 # noncoding, not CpG
(295707,295707)(33978373,33978373)(33978373,33978374)   4       33978373        2       0.351489 # noncoding, not CpG
(294702,294702)(33851203,33851203)(33851203,33851204)   4       33851203        2       0.354574 # noncoding, not CpG
(288733,288733)(32989837,32989837)(32989837,32989838)   4       32989837        2       0.365599 # noncoding, not CpG
(294800,294800)(33864292,33864292)(33864292,33864293)   4       33864292        2       0.367998 # noncoding, not CpG
(288903,288903)(33021517,33021517)(33021517,33021518)   4       33021517        2       0.371277 # noncoding, not CpG
(288861,288861)(33014876,33014876)(33014876,33014877)   4       33014876        2       0.382262 # noncoding, not CpG
(295744,295744)(33979932,33979932)(33979932,33979933)   4       33979932        2       0.394456 # noncoding, not CpG
(292638,292638)(33525518,33525518)(33525518,33525519)   4       33525518        2       0.407956 # noncoding, not CpG
(294691,294691)(33850679,33850679)(33850679,33850680)   4       33850679        2       0.442180 # noncoding, not CpG
(296462,296462)(34050117,34050117)(34050117,34050118)   4       34050117        2       0.472807 # noncoding, not CpG
(296461,296461)(34050116,34050116)(34050116,34050117)   4       34050116        2       0.484182 # noncoding, not CpG
(294686,294686)(33850312,33850312)(33850312,33850313)   4       33850312        2       0.581519 # noncoding, not CpG
(294690,294690)(33850622,33850622)(33850622,33850623)   4       33850622        2       0.605379 # noncoding, not CpG
(294682,294682)(33849547,33849547)(33849547,33849548)   4       33849547        2       0.656030 # CpG, between exon 4,5

~/programs/angsd/angsd \
  -GL 1 \
  -doGeno 2 \
  -doPost 1 \
  -doCounts 1 \
  -doMajorMinor 3 \
  -doMaf 1 \
  -doDepth 1 \
  -setMinDepthInd 4 \
  -minQ 30 \
  -minMapQ 30 \
  -r NC_044574.1:33979217 \
  -sites /xdisk/mcnew/finches/dannyjackson/finches/referencelists/allsnps.sites_headless.mafs \
  -bam /xdisk/mcnew/finches/dannyjackson/finches/referencelists/crapostbams.txt \
  -anc /xdisk/mcnew/finches/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa \
  -out single_site_post


~/programs/angsd/angsd \
  -GL 1 \
  -doGeno 2 \
  -doPost 1 \
  -doCounts 1 \
  -doMajorMinor 3 \
  -doMaf 1 \
  -doDepth 1 \
  -setMinDepthInd 4 \
  -minQ 30 \
  -minMapQ 30 \
  -r NC_044574.1:33979217 \
  -sites /xdisk/mcnew/finches/dannyjackson/finches/referencelists/allsnps.sites_headless.mafs \
  -bam /xdisk/mcnew/finches/dannyjackson/finches/referencelists/craprebams.txt \
  -anc /xdisk/mcnew/finches/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa \
  -out single_site_pre

zcat single_site_post.mafs.gz

zcat single_site_post.geno.gz 
zcat single_site_pre.geno.gz 



## FOR
cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/genelist/gene_analyses/indv_genes/for
FST="/xdisk/mcnew/finches/dannyjackson/finches/analyses/fst/forpre_forpost/1/forpre_forpost.1.fst"
### TENM3
#### NC_044574.1:32971132-34255509
##### CHR 4
awk -F'\t' 'BEGIN {OFS="\t"} NR==1 || ($2 == 4 && $3 >= 32971132 && $3 <= 34255509)' $FST > TENM3.fst.tsv
sort -k5 TENM3.fst.tsv | tail -n 20

(295145,295145)(33901150,33901150)(33901150,33901151)   4       33901150        2       0.341189 # maybe
(293464,293464)(33649395,33649395)(33649395,33649396)   4       33649395        2       0.347293 # not coding, not CpG
(289895,289895)(33149057,33149057)(33149057,33149058)   4       33149057        2       0.347509 # CpG, between exon 12 and 13
(296698,296698)(34066846,34066846)(34066846,34066847)   4       34066846        2       0.350949 # not coding, not CpG
(296390,296390)(34039585,34039585)(34039585,34039586)   4       34039585        2       0.353533 # not coding, not CpG
(293494,293494)(33652369,33652369)(33652369,33652370)   4       33652369        2       0.357574 # not coding, not CpG
(297508,297508)(34166619,34166619)(34166619,34166620)   4       34166619        2       0.363279 # not coding, not CpG
(293496,293496)(33652654,33652654)(33652654,33652655)   4       33652654        2       0.363988 # not coding, not CpG
(293424,293424)(33646370,33646370)(33646370,33646371)   4       33646370        2       0.369249 # not coding, not CpG
(295481,295481)(33941284,33941284)(33941284,33941285)   4       33941284        2       0.381590 # not coding, not CpG
(293536,293536)(33655484,33655484)(33655484,33655485)   4       33655484        2       0.385297 # not coding, not CpG
(296745,296745)(34072314,34072314)(34072314,34072315)   4       34072314        2       0.390856 # not coding, not CpG
(293439,293439)(33647325,33647325)(33647325,33647326)   4       33647325        2       0.391125 # not coding, not CpG
(296835,296835)(34083956,34083956)(34083956,34083957)   4       34083956        2       0.392178 # not coding, not CpG
(296713,296713)(34067655,34067655)(34067655,34067656)   4       34067655        2       0.400474 # not coding, not CpG
(296736,296736)(34069815,34069815)(34069815,34069816)   4       34069815        2       0.412365 # CpG site, between exon 3 and 4
(293493,293493)(33652266,33652266)(33652266,33652267)   4       33652266        2       0.416604 # CpG site, between exon 6 and 7
(293514,293514)(33653784,33653784)(33653784,33653785)   4       33653784        2       0.429384 # CpG site, between exon 6 and 7
(296678,296678)(34063771,34063771)(34063771,34063772)   4       34063771        2       0.487079 # not coding, not CpG


~/programs/angsd/angsd \
  -GL 1 \
  -doGeno 2 \
  -doPost 1 \
  -doCounts 1 \
  -doMajorMinor 3 \
  -doMaf 1 \
  -doDepth 1 \
  -setMinDepthInd 4 \
  -minQ 30 \
  -minMapQ 30 \
  -r NC_044574.1:33149057 \
  -sites /xdisk/mcnew/finches/dannyjackson/finches/referencelists/allsnps.sites_headless.mafs \
  -bam /xdisk/mcnew/finches/dannyjackson/finches/referencelists/forpostbams.txt \
  -anc /xdisk/mcnew/finches/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa \
  -out single_site_post


~/programs/angsd/angsd \
  -GL 1 \
  -doGeno 2 \
  -doPost 1 \
  -doCounts 1 \
  -doMajorMinor 3 \
  -doMaf 1 \
  -doDepth 1 \
  -setMinDepthInd 4 \
  -minQ 30 \
  -minMapQ 30 \
  -r NC_044574.1:33149057 \
  -sites /xdisk/mcnew/finches/dannyjackson/finches/referencelists/allsnps.sites_headless.mafs \
  -bam /xdisk/mcnew/finches/dannyjackson/finches/referencelists/forprebams.txt \
  -anc /xdisk/mcnew/finches/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa \
  -out single_site_pre

zcat single_site_post.mafs.gz

zcat single_site_post.geno.gz 
zcat single_site_pre.geno.gz 


## PAR
cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/genelist/gene_analyses/indv_genes/par
FST="/xdisk/mcnew/finches/dannyjackson/finches/analyses/fst/parpre_parpost/1/parpre_parpost.1.fst"
### TENM3
#### NC_044574.1:32971132-34255509
##### CHR 4
awk -F'\t' 'BEGIN {OFS="\t"} NR==1 || ($2 == 4 && $3 >= 32971132 && $3 <= 34255509)' $FST > TENM3.fst.tsv
sort -k5 TENM3.fst.tsv | tail -n 20

(291066,291066)(33290428,33290428)(33290428,33290429)   4       33290428        2       0.363123 # not coding, not CpG
(298112,298112)(34248653,34248653)(34248653,34248654)   4       34248653        2       0.363803 # CpG, between exon 1 and 2
(291078,291078)(33291644,33291644)(33291644,33291645)   4       33291644        2       0.368302 # not coding, not CpG
(294217,294217)(33758634,33758634)(33758634,33758635)   4       33758634        2       0.370254 # not coding, not CpG
(291065,291065)(33290382,33290382)(33290382,33290383)   4       33290382        2       0.371981 # not coding, not CpG
(291077,291077)(33291567,33291567)(33291567,33291568)   4       33291567        2       0.372455 # not coding, not CpG
(297929,297929)(34224005,34224005)(34224005,34224006)   4       34224005        2       0.374003 # not coding, not CpG
(291076,291076)(33291533,33291533)(33291533,33291534)   4       33291533        2       0.375889 # not coding, not CpG
(291080,291080)(33291805,33291805)(33291805,33291806)   4       33291805        2       0.378515 # not coding, not CpG
(291081,291081)(33291835,33291835)(33291835,33291836)   4       33291835        2       0.387184 # not coding, not CpG
(291082,291082)(33291842,33291842)(33291842,33291843)   4       33291842        2       0.389059 # not coding, not CpG
(291045,291045)(33288670,33288670)(33288670,33288671)   4       33288670        2       0.406305 # not coding, not CpG
(298012,298012)(34235390,34235390)(34235390,34235391)   4       34235390        2       0.408015 # CpG, between exon 1 and 2
(297882,297882)(34217203,34217203)(34217203,34217204)   4       34217203        2       0.411867 # not coding, not CpG
(297910,297910)(34220802,34220802)(34220802,34220803)   4       34220802        2       0.412324 # not coding, not CpG
(297954,297954)(34227944,34227944)(34227944,34227945)   4       34227944        2       0.432350 # not coding, not CpG
(293844,293844)(33700865,33700865)(33700865,33700866)   4       33700865        2       0.456431 # not coding, not CpG
(290379,290379)(33210901,33210901)(33210901,33210902)   4       33210901        2       0.459470 # not coding, not CpG
(297941,297941)(34225827,34225827)(34225827,34225828)   4        v        2       0.485842 # CpG, between exon 1 and 2


~/programs/angsd/angsd \
  -GL 1 \
  -doGeno 2 \
  -doPost 1 \
  -doCounts 1 \
  -doMajorMinor 3 \
  -doMaf 1 \
  -doDepth 1 \
  -setMinDepthInd 4 \
  -minQ 30 \
  -minMapQ 30 \
  -r NC_044574.1:34248653 \
  -sites /xdisk/mcnew/finches/dannyjackson/finches/referencelists/allsnps.sites_headless.mafs \
  -bam /xdisk/mcnew/finches/dannyjackson/finches/referencelists/parpostbams.txt \
  -anc /xdisk/mcnew/finches/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa \
  -out single_site_post


~/programs/angsd/angsd \
  -GL 1 \
  -doGeno 2 \
  -doPost 1 \
  -doCounts 1 \
  -doMajorMinor 3 \
  -doMaf 1 \
  -doDepth 1 \
  -setMinDepthInd 4 \
  -minQ 30 \
  -minMapQ 30 \
  -r NC_044574.1:34248653 \
  -sites /xdisk/mcnew/finches/dannyjackson/finches/referencelists/allsnps.sites_headless.mafs \
  -bam /xdisk/mcnew/finches/dannyjackson/finches/referencelists/parprebams.txt \
  -anc /xdisk/mcnew/finches/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa \
  -out single_site_pre

zcat single_site_post.mafs.gz

zcat single_site_post.geno.gz 
zcat single_site_pre.geno.gz 


# plot change in freq of CpG by species by gene
ParallelGenes_BMPER_NUDT4_TENM3.csv

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)

# Read the data
df <- read.csv("ParallelGenes_BMPER_NUDT4_TENM3.csv")

# Filter out any rows missing Pre-freq or Post-freq
df_clean <- df %>%
  filter(!is.na(`Pre.freq`), !is.na(`Post.freq`))

# Create a long-format dataframe for plotting pre and post frequencies
df_long <- df_clean %>%
  pivot_longer(cols = c(`Pre.freq`, `Post.freq`),
               names_to = "Timepoint",
               values_to = "Frequency") %>%
  mutate(Timepoint = factor(Timepoint, levels = c("Pre.freq", "Post.freq")))

# Add unique ID for each line (to group by site and species)
df_long <- df_long %>%
  mutate(SiteID = paste(Gene, Species, Site, sep = "_"))

# Define custom colors for species
species_colors <- c(
  "CRA" = "#F47F7C",
  "PAR" = "#A5C965",
  "FOR" = "#48AEAD"
)

# Plot
ggplot(df_long, aes(x = Timepoint, y = Frequency, group = SiteID, color = Species)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 2) +
  facet_grid(rows=vars(Gene), scales = "free_y") +
  scale_color_manual(values = species_colors) +
  theme_minimal(base_size = 14) +
  labs(title = "Pre vs Post Frequencies per Gene and Species",
       x = "Timepoint",
       y = "Allele Frequency",
       color = "Species") +
  theme(strip.text = element_text(face = "bold"))

  ggsave(filename = paste0("CpG.freqchange.pdf"), 
       width = 5, height = 10, units = "in")




### FORTIS GO SIGNIFICANT GENES
# CHST4, CELF1, TERF2IP, PRTFDC1, ADAT1, FARS2, MCM9
cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/genelist/gene_analyses/indv_genes/for

FST="/xdisk/mcnew/finches/dannyjackson/finches/analyses/fst/forpre_forpost/1/forpre_forpost.1.fst"

# CHST4
# NC_044581.1:60789-89600	; CHR 11
awk -F'\t' 'BEGIN {OFS="\t"} NR==1 || ($2 == 11 && $3 >= 60789 && $3 <= 89600)' $FST > CHST4.fst.tsv
sort -k5 CHST4.fst.tsv | tail -n 20

(420,420)(78911,78911)(78911,78912)     11      78911   2       0.126085
(314,314)(61622,61622)(61622,61623)     11      61622   2       0.149282
(402,402)(75315,75315)(75315,75316)     11      75315   2       0.156070
(388,388)(73001,73001)(73001,73002)     11      73001   2       0.164171
(389,389)(73028,73028)(73028,73029)     11      73028   2       0.165175
(455,455)(85353,85353)(85353,85354)     11      85353   2       0.184807 # Not coding, not CpG
(434,434)(81274,81274)(81274,81275)     11      81274   2       0.277235 # Not coding, not CpG
(454,454)(85342,85342)(85342,85343)     11      85342   2       0.286443 # Not coding, not CpG
(334,334)(65506,65506)(65506,65507)     11      65506   2       0.287269 # Not coding, not CpG
(343,343)(66773,66773)(66773,66774)     11      66773   2       0.364876 # Not coding, not CpG

# CELF1
# NC_044575.1:22468787-22531900	; CHR 5

# TERF2IP
# NC_044581.1:54159-56952	; CHR 11

# PRTFDC1
# NC_044572.1:17462205-17505243	; CHR 2


### These don't have any interpretable meaning, at least not yet
# ADAT1
# NC_044581.1:21238-45913	; CHR 11

# FARS2
# NC_044572.1:67051237-67283155 ; CHR 2

# MCM9
# NC_044573.1:63684526-63732610	; CHR 3