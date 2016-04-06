# scripts for atherosclerosis 16S rRNA gene sequencing project

### make phylogenetic tree

This code is adopted from https://github.com/ggloor/miseq_bin/blob/master/macqiime_R_plots.txt

Align sequences with muscle (http://www.drive5.com/muscle/manual)

```
nohup ./muscle -in data/OTU_seed_seqs.fa -out data/all_seed_OTUs_bad.mfa > muscle_nohup.out 2>&1&
```

Fix the formatting

``` 
awk '/^>/{gsub(/lcl\|[0-9]+\|num\|[0-9]+\|OTU\|/, "")}; 1' all_seed_OTUs_bad.mfa > all_seed_OTUs.mfa
rm all_seed_OTUs_bad.mfa
```

Build tree with FastTree

```
FastTree -nt all_seed_OTUs.mfa > OTU_seqs.tre
```

## biplot, bargraph dendrogram, ALDEx2

## UniFrac PCoA


