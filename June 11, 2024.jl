June 11, 2024
FST peak analyses:

# FOR 
> subset(fst, fst > 0.7)
                                                       region         chr
3936651 (187634,187634)(25267982,25267982)(25267982,25267983) NC_044575.1
5445766 (134702,134702)(16843164,16843164)(16843164,16843165) NC_044580.1
5445767 (134703,134703)(16843166,16843166)(16843166,16843167) NC_044580.1
8139248 (279030,279030)(48092710,48092710)(48092710,48092711) NC_044601.1
          midPos Nsites      fst
3936651 25267982      2 0.824873
5445766 16843164      2 0.716541
5445767 16843166      2 0.744155
8139248 48092710      2 0.701127

# Only one region showed up twice, NC_044580.1 positions 16843164 and 16843166. Both of these are within the gene: synemin SYNM (NC_044580.1:16846280-16866830)
This is an intermediate filament protein originally found in avian smooth muscles. Gene Cards says it is a "cytoskeletal protein that confers resistance to mechanical stress"








# CRA
subset(fst, fst > 0.7)

                                                           region     chr
335661      (335661,335661)(37907290,37907290)(37907290,37907291) 44571.1
335663      (335663,335663)(37907336,37907336)(37907336,37907337) 44571.1
3048343 (820574,820574)(104906414,104906414)(104906414,104906415) 44573.1
3048345 (820576,820576)(104906573,104906573)(104906573,104906574) 44573.1
3048362 (820593,820593)(104908569,104908569)(104908569,104908570) 44573.1
3048369 (820600,820600)(104909982,104909982)(104909982,104909983) 44573.1
3048373 (820604,820604)(104910683,104910683)(104910683,104910684) 44573.1
4083631     (356727,356727)(50553761,50553761)(50553761,50553762) 44575.1
4416458     (225338,225338)(27896405,27896405)(27896405,27896406) 44576.1
6735999     (548241,548241)(65558014,65558014)(65558014,65558015) 44586.1
6842675           (58407,58407)(6978972,6978972)(6978972,6978973) 44587.1
7861887     (109171,109171)(13085797,13085797)(13085797,13085798) 44601.1
           midPos Nsites      fst
335661   37907290      2 0.724117
335663   37907336      2 0.715765
3048343 104906414      2 0.769188
3048345 104906573      2 0.738477
3048362 104908569      2 0.729066
3048369 104909982      2 0.865019
3048373 104910683      2 0.728000
4083631  50553761      2 0.770351
4416458  27896405      2 0.720562
6735999  65558014      2 0.700035
6842675   6978972      2 0.716074
7861887  13085797      2 0.789113

# NC_044571.1   pos: 37907290 and 37907336
# this is not within a gene, but is just after:
NC_044571.1:37901151-37906295 testis-expressed protein 30 isoform X1	TEX30
# NCBI says: "Predicted to enable hydrolase activity."
# Gene cards says: "Diseases associated with TEX30 include Spermatogenic Failure 25 and Hypotonia-Cystinuria Syndrome"

# and before: 
NC_044571.1:37937016-37945300 LOW QUALITY PROTEIN: protein-lysine methyltransferase METTL21C	METTL21C
# NCBI says: "Enables heat shock protein binding activity and protein-lysine N-methyltransferase activity. Involved in protein methylation. Located in nucleus. Part of protein-containing complex."
# Gene cards says: "Diseases associated with METTL21C include Myostatin-Related Muscle Hypertrophy and Severe Congenital Neutropenia 6."

# 44573.1       pos: 104906414, 104906573, 104908569, 104909982, 104910683
Between these two:
NC_044573.1:104387453-104741091 glutamate receptor ionotropic, kainate 2	GRIK2
# NCBI says: "Glutamate receptors are the predominant excitatory neurotransmitter receptors in the mammalian brain and are activated in a variety of normal neurophysiologic processes. "
# gene cards says: "Diseases associated with GRIK2 include Neurodevelopmental Disorder With Impaired Language And Ataxia And With Or Without Seizures and Intellectual Developmental Disorder, Autosomal Recessive 6. Among its related pathways are Presynaptic function of Kainate receptors and Transmission across Chemical Synapses."
NC_044573.1:105878099-105927569 E3 ubiquitin-protein ligase HACE1 isoform X1	HACE1
# NCBI says: "The encoded protein is involved in specific tagging of target proteins, leading to their subcellular localization or proteasomal degradation. The protein is a potential tumor suppressor and is involved in the pathophysiology of several tumors, including Wilm's tumor."
# Gene cards says: "Diseases associated with HACE1 include Spastic Paraplegia And Psychomotor Retardation With Or Without Seizures and Charcot-Marie-Tooth Disease, Axonal, Type 2W. Among its related pathways are Class I MHC mediated antigen processing and presentation and Innate Immune System."




# PAR

# just some quick looks at the output:
fst <- read.csv('slidingwindow_singlesnps_fst.txt', sep ='\t')
max(fst$fst)
subset(fst, fst > 0.7)

                                                       region     chr   midPos
5125132 (229970,229970)(27663929,27663929)(27663929,27663930) 44578.1 27663929
5125133 (229971,229971)(27663934,27663934)(27663934,27663935) 44578.1 27663934
        Nsites      fst
5125132      2 0.710722
5125133      2 0.701462

subset(fst, fst > 0.7)

# NC_044578.1:27663929-27663934
Between the following two genes:
NC_044578.1:27643755-27644801 LOW QUALITY PROTEIN: immediate early response gene 5 protein	IER5
# NCBI says: "This gene encodes a protein that is similar to other immediate early response proteins. In the mouse, a similar gene may play an important role in mediating the cellular response to mitogenic signals. Studies in rats found the expression of a similar gene to be increased after waking and sleep deprivation."
# Gene cards doesn't say anything useful

NC_044578.1:27705862-27829414 voltage-dependent R-type calcium channel subunit alpha-1E isoform X1	CACNA1E
# NCBI says: "Voltage-dependent calcium channels are multisubunit complexes consisting of alpha-1, alpha-2, beta, and delta subunits in a 1:1:1:1 ratio. These channels mediate the entry of calcium ions into excitable cells, and are also involved in a variety of calcium-dependent processes, including muscle contraction, hormone or neurotransmitter release, gene expression, cell motility, cell division and cell death. This gene encodes the alpha-1E subunit of the R-type calcium channels, which belong to the 'high-voltage activated' group that maybe involved in the modulation of firing patterns of neurons important for information processing."
# Gene cards says "Diseases associated with CACNA1E include Developmental And Epileptic Encephalopathy 69 and Van Der Woude Syndrome 1. Among its related pathways are DREAM Repression and Dynorphin Expression and TCR Signaling (Qiagen)."



> subset(fst, fst > 0.65)

                                                           region     chr
757361      (757361,757361)(87171218,87171218)(87171218,87171219) 44571.1
757362      (757362,757362)(87171219,87171219)(87171219,87171220) 44571.1
2568818     (304218,304218)(38140696,38140696)(38140696,38140697) 44573.1
3046746     (782146,782146)(98996648,98996648)(98996648,98996649) 44573.1
3128069 (863469,863469)(108909526,108909526)(108909526,108909527) 44573.1
3979790     (189914,189914)(25219440,25219440)(25219440,25219441) 44575.1
4494346     (229325,229325)(27901632,27901632)(27901632,27901633) 44576.1
5125132     (229970,229970)(27663929,27663929)(27663929,27663930) 44578.1
5125133     (229971,229971)(27663934,27663934)(27663934,27663935) 44578.1
5349271     (196554,196554)(23290519,23290519)(23290519,23290520) 44579.1
           midPos Nsites      fst
757361   87171218      2 0.667568
757362   87171219      2 0.668538
2568818  38140696      2 0.672893
3046746  98996648      2 0.692131
3128069 108909526      2 0.662883
3979790  25219440      2 0.679678
4494346  27901632      2 0.675304
5125132  27663929      2 0.710722
5125133  27663934      2 0.701462
5349271  23290519      2 0.678239
# additional region of:
# NC_044571.1:87171218-87171219
# between these two genes:
NC_044571.1:87095605-87147397 transcriptional enhancer factor TEF-3	TEAD4
# NCBI says: "This gene product is a member of the transcriptional enhancer factor (TEF) family of transcription factors, which contain the TEA/ATTS DNA-binding domain. It is preferentially expressed in the skeletal muscle, and binds to the M-CAT regulatory element found in promoters of muscle-specific genes to direct their gene expression."
# Gene cards says: "Diseases associated with TEAD4 include Sveinsson Chorioretinal Atrophy. Among its related pathways are Gene expression (Transcription) and ERK Signaling."

NC_044571.1:87153063-87322660 tetraspanin-9 isoform X1	TSPAN9
# NCBI says: "The protein encoded by this gene is a member of the transmembrane 4 superfamily, also known as the tetraspanin family. Most of these members are cell-surface proteins that are characterized by the presence of four hydrophobic domains. The proteins mediate signal transduction events that play a role in the regulation of cell development, activation, growth and motility. "
# Gene cards didn't say anything useful
