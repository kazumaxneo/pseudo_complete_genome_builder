######################################################################################## 
pseudo_complete_genome_builder 0.6

A Perl scripts to create pseudo complete genome from fragmented sequences(contigs).   

pseudo_complete_genome_builder: 

Kazuma Uesaka, and Kunio Ihara  



Input: 
  Complete reference genome(fasta) and contig or scaffolds(fasta)
Outnput:	
  complete genome sequeces that are ordered and oriented by the reference genome sequences

Usage:  
  perl pseudo_complete_genome_builder_v0.6.pl


########################################################################################

    
## Requirements  
- SAMTools  (version >= 1.3.1)   
- BWA (version >= 0.7.17)  
- circos (v0.67. only required if drawing indel map))  


Install Anaconda3 or Miniconda3. Then,  
```
conda env create -f=env.yml
conda activate pseudo_complete_genome_builder_v0.6
```
    

## Source
```
cd $HOME 
git clone git@github.com:kazumaxneo/pseudo_complete_genome_builder.git
cd pseudo_complete_genome_builder/
perl pseudo_complete_genome_builder_v0.6.pl
```
    


## create pseundo complete genome. If you uses spades assembler for de novo assembly, type
```
perl pseudo_complete_genome_builder_v0.6.pl -f reference.fasta -s scaffolds.fasta -max 50 -min 40 -n 1 -p 3000 -l 3000 -t 150 -a spades
```

##If you uses unicycler assembler for de novo assembly, type
```
perl pseudo_complete_genome_builder_v0.6.pl -f reference.fasta -s assembly.fasta -max 1.1 -min 0.9 -n 1 -p 3000 -l 3000 -t 150 -a unicycler
```

## Test run
```

```  

First, pseudo complete genome was created, then thier pairwise alignment with was performed with minimap2. Finall, their sequence difference, including SNVs, small indels, and large indels, was called using paftools.


## Licence ##

GPL v3.


