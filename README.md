# fast pangenomes
This repository currently relates to two piblications:
* Hernandez-Salmeron JE, Irani T, Moreno-Hagelsieb G (2023) Fast genome-based delimitation of Enterobacterales species. PLoS One 18: e0291492. DOI: [10.1371/journal.pone.0291492](https://dx.doi.org/10.1371/journal.pone.0291492)
* Hernandez-Salmeron JE, Moreno-Hagelsieb G (2022) FastANI, Mash and Dashing equally differentiate between Klebsiella species. PeerJ 10: e13784. DOI: [10.7717/peerj.13784](https://dx.doi.org/10.7717/peerj.13784)


This is a set of programs, scripts, etc, needed to quickly organize prokaryotic genomes into "species"-level groups (pangenomes, or, plainly, to try and accomodate new genomes into the hierarchies existing elsewhere. So, this needs a scripted pipeline, and an explanation of what the pipeline would accomplish. We will get there at some point.

## Scripts/programs

Note: except for obvious shell scripts, like ".zsh" or ".sh" or ".bash", you should make these programs and scripts executables (under any unix, like linux or darwin):

`chmod +x *.pl *.R`

Yes, the ".R" are R scripts, but they run as executables.

## Needs:
* mash:

  Though it's possible to install mash with homebrew in mac, the version, as of this writing, is outdated. Thus, we installed it from: https://github.com/marbl/Mash

* R with packages:
  
  "cluster", "MCMCpack", "ape","reshape2", "fastcluster" for clustering; "tidyverse" and "cutpointr" for optimising your cutoffs for species and/or genus level groups. The package "fastcluster" is only needed if you want to cluster using a method other than diana (we recommend diana)

  
* genomes to cluster:

  This means the DNA sequences of such genomes (normally ending with ".fna." Our programs and scripts will work with other formats (or so we hope), but we think it's better to keep things organized. Thus we tend to use ".fna" for genome, DNA, sequences. We normally take them from NCBI, actually, we have downloaded all of the ones under RefSeq for quite a while, but that's becoming too much. Either way, if what you want is to check where your newly sequenced genomes belong, and have already some suspitions, download relevant genomes close by and around what you want to check. For example, if your genome is within cyanobacteria, maybe download at least a few of each species represented from a relevant genome database.

## Instructions:
1. Put the genomes into a single directory. Our scripts work with directories using this format:

   fna-[group]

   where "group can mean anything you like. For example, we've used "fna-Klebsiella", "fna-Enterobacterales", "fna-Cyanobacteriota"

2. The script "masher.zsh" needs you to specify the "last name" of such a directory:

   `zsh masher.zsh Klebsiella`

   The script will run mash to produce the necessary sketches and then the whole all _vs_ all comparison of the genomes in fna-Klebsiella. It will store the results, the all _vs_ all mash distances, at a subdirectory called "Mash". In the Klebsiella case, the file will be called "mash-Klebsiella.tbl.bz2"

3. Run clusterGenomes.pl

   The program automatically recognizes the format of the mash file, fixes it to run properly with R, then calls R to cluster and produce the hierarchy. The hierarchy looks like a phylogenetic tree, and it's the best way to check where your genomes belong in the hierarchy of the organisms you suspected. After clustering, the program saves several useful files. The main ones you might like to check are: a ".tree" file that can be read by almost any program for phylogenetics/trees, also readable by some packages in R, and a ".Rds" which is an R 'object' that retains all the properties of the hierarchy produced by R, and reads very quickly and efficiently into R for further examination, plotting, etc.

4. We will add, later, the programs we use to decide on cutoffs, and to plot the results.
