CRISPR-Cas9 editing
===


# Table of contents




# Material required


## Plasmids

The construct is split into one module vector set (pCBC-DT1T2) contaning the guide RNAs (gRNAs) and a binary vector based on pCAMBIA (pHHE401E). The cloning is based on a primary PCR to add the proper target sequences on the gRNA on the module vector, *BsaI* restriction sites are used to ligate the module vector set to the binary vector and remove meanwhile its spectinomycin (SpR) resistance cassette [Xing et al 2014](https://bmcplantbiol.biomedcentral.com/articles/10.1186/s12870-014-0327-y).

*Cas9* expression is driven by the promoter of the egg cell-specific *EC1.1* ([AT1G76750](https://www.arabidopsis.org/servlets/TairObject?id=29908&type=locus)) gene and the enhancer  *EC1.2* gene ([AT2G21740](https://www.arabidopsis.org/servlets/TairObject?accession=locus:2052536)). This tissue-specific expression of *Cas9* allows to obtain T1 homozygous or biallelic plants instead of mosaic plants. The terminator of the *Pisum sativum rbcS E9* gene was tested according to previous observation ([Sarrion-Perdigones et al 2013](http://www.plantphysiol.org/content/162/3/1618.short)) and found more efficient than coupled with a NOS terminator. See reference [Wang et al. 2015](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0715-0).


### Module vector pCBC-DT1T2

More information on https://www.addgene.org/50590/

This vector contains 2 single guide RNAs (sgRNAs) target sequences (Ns). More targets can be added by adding primers with specific overhanging sequence for Golden Gate assembly. See more details about Golden Gate assembly on [NEB website](https://international.neb.com/applications/cloning-and-synthetic-biology/dna-assembly-and-cloning/golden-gate-assembly).

The main principle to understand is that a type IIS restriction enzyme such as *BsaI*, used here, is able to cut besides its actual recognition site, allowing to create non-palindromic overhangs which allow directional cloning.

*BsaI* restriction site:

```
5'...GGTCTC(N)1...3'
3'...CCAGAG(N)5...5'
```

The enzyme cuts 1 nucleotide directly upstream of the recognition side on the + strand but 5 nucleotides upstream on the - strand (highlighted in yellow in figure below). The 5 nucleotides overhang can be chosen and included in any construct to allow ligation. T4 ligase is used to ligate the corresponding cohesive ends (it is also able to ligate blunt ends but with much lower efficiency).


![](images/BsaI_sites.PNG) 

Note that the vector pCBC-DT1T2 does not contain the promoter of the first sgRNA and the second sgRNA contains neither its scaffolding RNA part nor its terminator. These sequences are actually in the pHHE401E binary vector and will be integrated in place after the digestion with *BsaI* enzyme and the ligation (GoldenGate cloning).


### pHHE401E binary vector

pHHE401E is based on pCambia backbone for Agrobacterium-mediated transformation. It contains the complementary parts of the pCBC-DT1T2 which ensures the proper expression of the 2 sgRNAs plus other features needed for Agrobacterium transformation and selection.

More information on https://www.addgene.org/71287/ 

## sgRNAs

The vector pCBC-DT1T2 can contain two target sites. These should be no more than 100 bp apart in the genome.
To generate 2 sgRNAs with each their target site (gRNA1 and gRNA2 indicated by blue arrows and highlighted in yellow in image below), 4 primers need to be designed. Two primers for the target site 1 and 2 primers for the target size 2 (see image below):

First sgRNA (green on picture below)

* DT1-BsF
* DT1-F0

Second sgRNA (red on picture below)

* DT2-BsR
* DT2-R0

These target sites should be 19 bp long be adjacent to a NGG site (N being any nucleotide), which is not in the primer sequence. This site is called the protospacer adjacent motif (PAM). Cas9 enzymes cuts theoretically 3 bp upstream of the PAM motif (always within the target region). The target sites are represented by Ns in the figure below. The PAM should always be originally located at the 5' end of the sgRNA target sequence (highlighted in yellow). 

For instance for the sequence `TCGAGAGAGAGCGTATTTCGGG`, the kmer `GGG` is the PAM motif located 5' of the sequence and the cut will therefore take place here: `TCGAGAGAGAGCGTATT|TTCGGG`. The sequence to include in the 2 primers will be therefore the 19-mer `TCGAGAGAGAGCGTATTTTC` (note the PAM is not included).

![](images/sgRNAs.PNG)

The 2 most outward primers indicated in red (DT1-BsF) and red (DT2-BsR) allow to add specific overhanging sites including a specific BsaI site and its corresponding overhang sequence present in the pHHE401E binary vector and allow directional cloning. The second 2 primers more inward (DT1-F0 and DT2-R0) contain a part of the vector pCBC-DT1T2. All primers contain their respective target sites (Ns). The PCR with all primers allow therefore to add in a first step the target sites with the primers DT1-F0 and DT2-R0 and in a second step the BsaI restriction site for GoldenGate cloning (DT1-BsF and DT2-BsR).

Both should contain the 19 bp of the target site (highlighted in yellow) plus flanking regions that are either used for Gateway cloning or for the merging to the vector pCBC-DT1T2.


Note: If the PAM is located on the reverse strand, the sequence used in the primer represents the reverse complement sequence of the + strand. See example below.

 ![](images/PAM_minus_strand.PNG)
 
 Note that the PAM should always be in 5' position of the target site in the pCBC-DT1T2 plasmid.

# Authors:

* [Johan Zicola](https://github.com/johanzi)

* [Emmanuel Tergemina](https://github.com/EmmanuelTergemina)
