# Independent literature review and scoring rationale

The review used non-genomic natural-history observations, parasite-removal experiments, behavioral observations, antibody assays, and an ecological review. No genomic, transcriptomic, or gene-expression findings were used. One such study appeared incidentally in broad search results and was excluded. Each GO term was judged independently; there was no target score distribution. The input order and identifiers were preserved.

## Host injury and physiological compensation

Adult flies are free living; their larvae parasitize birds. Early instars occupy nestling nostrils and other tissues, while later instars shelter in nest material and feed externally, predominantly at night. Blood removal and tissue consumption cause anemia, wounds, impaired growth, nasal/beak deformation, and potentially death. Parasite reduction improves fledging success. These observations make replacement of blood cells, iron handling, hemostasis, oxygen delivery, epithelial repair, extracellular-matrix repair, and vascular repair strong physiological predictions. They do not demonstrate that every one of those compensatory pathways has been measured in these birds.

Key sources:
- Fessl, Sinclair and Kleindorfer (2006), *The life-cycle of Philornis downsi (Diptera: Muscidae) parasitizing Darwin's finches and its impacts on nestling survival*. https://pubmed.ncbi.nlm.nih.gov/16899139/
- Koop et al. (2011), *Experimental demonstration of the fitness consequences of an introduced parasite of Darwin's finches*. https://pmc.ncbi.nlm.nih.gov/articles/PMC3092749/

## Resistance and immune defense

An experimental study found increased fly-binding antibodies in parasitized medium ground finch mothers, but not fathers or young nestlings. The inverse association between maternal antibodies and parasite abundance was marginal; antibodies did not rescue reproductive success. Thus adaptive immunity is directly relevant, while its individual cellular components are strong predictions rather than established effective protection. The study discusses inflammation, itching, defensive peptides and parasite-feeding disruption as potential mechanisms. Innate responses, cytokines, complement and leukocyte functions therefore remain plausible; pathogen-specific antiviral machinery is much less specifically connected. Wounds also provide a plausible secondary bacterial exposure route, without establishing bacterial disease in every host.

Key source:
- Koop et al. (2013), *Experimental demonstration of a parasite-induced immune response in wild birds: Darwin's finches and introduced nest flies*. https://pmc.ncbi.nlm.nih.gov/articles/PMC3930052/

## Tolerance and behavior

Mockingbird experiments support tolerance through increased nestling begging and parental provisioning, compensating for parasite damage without necessarily reducing parasite density. Tolerance varies with environmental resources and is poorer in drought years. Nutrient sensing, energy balance, carbohydrate metabolism and resource allocation are consequently relevant; generic metabolism alone provides weaker specificity.

Preening, picking at nestling nares, nest sanitation and avoidance movements have been observed, but evidence of effective behavioral resistance is limited. Females may reduce brooding and stand more, potentially altering nestling thermal conditions. Proposed leaf rubbing/self-medication and experimentally supplied insecticide-treated nesting material provide behavioral context, not evidence that intrinsic xenobiotic detoxification is a central defense. General neural development is too nonspecific to receive a strong score merely because behavior requires a nervous system.

Key sources:
- Knutie et al. (2016), *Galapagos mockingbirds tolerate introduced parasites that affect Darwin's finches*. https://pubmed.ncbi.nlm.nih.gov/28792593/
- McNew et al. (2019), *Annual environmental variation influences host tolerance to parasites*. https://pmc.ncbi.nlm.nih.gov/articles/PMC6408884/
- McNew and Clayton (2018), *Alien Invasion: Biology of Philornis Flies Highlighting Philornis downsi, an Introduced Parasite of Galapagos Birds*. https://doi.org/10.1146/annurev-ento-020117-043103

## Applying the scale

Score 5 was reserved for broad immune mechanisms directly captured by experimentally documented antibody responses. Score 4 identifies their core cellular components or particularly strong predictions from blood loss, wounds, and compensatory nutrition. Score 3 indicates supporting physiology, signaling, repair or resource use. Score 2 indicates generic cellular dependencies or indirect links. Score 1 indicates no specific connection, including many detailed embryonic, reproductive, antiviral, neural-development and nucleic-acid-maintenance processes. Regulation in either direction can be relevant; a negative-regulation term was not automatically penalized when balancing the mechanism could matter. Avian thrombocytes rather than mammalian platelets mediate hemostasis, so platelet GO labels were interpreted as conserved hemostatic functions rather than proof of mammalian cell anatomy.

GO:0160240 was supplied without a name. Its identity was checked as RNA polymerase II transcription initiation surveillance and scored as a specific housekeeping process without a particular parasite link: https://amigo.geneontology.org/amigo/term/GO%3A0160240 .
