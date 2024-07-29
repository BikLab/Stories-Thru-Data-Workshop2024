00-create-phyloseq-object
================
desantiago
2024-06-04

### Load required libraries

    ## Loading required package: phyloseq

    ## Loading required package: tidyr

    ## Warning: package 'tidyr' was built under R version 4.3.2

    ## Loading required package: devtools

    ## Loading required package: usethis

    ## Warning: package 'usethis' was built under R version 4.3.2

    ## Loading required package: qiime2R

    ## Warning: package 'ggplot2' was built under R version 4.3.2

### Read artifact files (.qza) using qiime2R

``` r
otus <- read_qza("data/03-dada2-independent-table.qza")
```

    ## Warning in asMethod(object): sparse->dense coercion: allocating vector of size
    ## 6.0 GiB

``` r
taxonomy <- read_qza("data/04-taxonomy-independent-classification.qza")
tree <- read_qza("data/05-phylogeny-independent-tree-midrooted.qza")
metadata <- read.delim2("data/metadata_SONGS_CABight_DDT_07-19-24 (1).txt", sep = "\t", row.names = 1)
```

### View QZA files

QZA files are zipped folders with many different pieces of information
including data provenance, format, version, and the data. We need to
extract out data from this file type

``` r
otu_df <- otus$data # we can view and save the actual data by specifying '$'
taxonomy_df <- taxonomy$data # save taxonomy info as a dataframe
phylo_tree <- tree$data # tree data is stored as a phyloseq object
```

### View taxonomy file

The taxonomy file is not formatted correctly. All the taxonomy is in one
column. Each taxonomic level has to be in its own column

``` r
head(taxonomy_df) 
```

    ##                         Feature.ID
    ## 1 00000ee244fdf5afdf2ffd124c23a320
    ## 2 0000a9d6f7dd5ee6cb6f2bed94784a39
    ## 3 00018b38663c5d40aa3c8e045f4d24aa
    ## 4 0002afd8dee5d314be3adc4cda9e46d2
    ## 5 0003175263e1b9530f1252de4f03e5ee
    ## 6 00033af273b8e729cf00555541979ad9
    ##                                                                                                                                                                        Taxon
    ## 1 D_0__Eukaryota;D_1__SAR;D_2__Alveolata;D_3__Ciliophora;D_4__Intramacronucleata;D_5__Conthreep;D_6__Plagiopylea;D_7__Odontostomatida;D_8__Epalxella;D_9__uncultured ciliate
    ## 2                                                D_0__Eukaryota;D_1__SAR;D_2__Rhizaria;D_3__Cercozoa;D_4__Vampyrellidae;D_5__uncultured;D_6__uncultured freshwater eukaryote
    ## 3  D_0__Eukaryota;D_1__SAR;D_2__Alveolata;D_3__Ciliophora;D_4__Intramacronucleata;D_5__Spirotrichea;D_6__Oligotrichia;D_7__Sinistrostrombidium;D_8__Strombidium paracalkinsi
    ## 4                                                                                                                                                                 Unassigned
    ## 5                                                            D_0__Eukaryota;D_1__SAR;D_2__Rhizaria;D_3__Cercozoa;D_4__Thecofilosea;D_5__uncultured;D_6__uncultured eukaryote
    ## 6                    D_0__Eukaryota;D_1__SAR;D_2__Stramenopiles;D_3__Ochrophyta;D_4__Diatomea;D_5__Bacillariophytina;D_6__Bacillariophyceae;D_7__uncultured marine eukaryote
    ##   Consensus
    ## 1         1
    ## 2         1
    ## 3         1
    ## 4         1
    ## 5         1
    ## 6         1

### Split taxonomy into different columns and replace NA’s with Unassigned

``` r
taxonomy_fixed_df <- taxonomy_df %>% separate_wider_delim(Taxon, delim = ";", names_sep = "", too_few = "align_start")
taxonomy_fixed_df[is.na(taxonomy_fixed_df)] <- "Unassigned" # rename NAs into unassigned
head(taxonomy_fixed_df)
```

    ## # A tibble: 6 × 25
    ##   Feature.ID      Taxon1 Taxon2 Taxon3 Taxon4 Taxon5 Taxon6 Taxon7 Taxon8 Taxon9
    ##   <chr>           <chr>  <chr>  <chr>  <chr>  <chr>  <chr>  <chr>  <chr>  <chr> 
    ## 1 00000ee244fdf5… D_0__… D_1__… D_2__… D_3__… D_4__… D_5__… D_6__… D_7__… D_8__…
    ## 2 0000a9d6f7dd5e… D_0__… D_1__… D_2__… D_3__… D_4__… D_5__… D_6__… Unass… Unass…
    ## 3 00018b38663c5d… D_0__… D_1__… D_2__… D_3__… D_4__… D_5__… D_6__… D_7__… D_8__…
    ## 4 0002afd8dee5d3… Unass… Unass… Unass… Unass… Unass… Unass… Unass… Unass… Unass…
    ## 5 0003175263e1b9… D_0__… D_1__… D_2__… D_3__… D_4__… D_5__… D_6__… Unass… Unass…
    ## 6 00033af273b8e7… D_0__… D_1__… D_2__… D_3__… D_4__… D_5__… D_6__… D_7__… Unass…
    ## # ℹ 15 more variables: Taxon10 <chr>, Taxon11 <chr>, Taxon12 <chr>,
    ## #   Taxon13 <chr>, Taxon14 <chr>, Taxon15 <chr>, Taxon16 <chr>, Taxon17 <chr>,
    ## #   Taxon18 <chr>, Taxon19 <chr>, Taxon20 <chr>, Taxon21 <chr>, Taxon22 <chr>,
    ## #   Taxon23 <chr>, Consensus <dbl>

### Merge files into a phyloseq object

convert your otu table, taxonomy files, and tree into phylseq object

``` r
# fix taxonomy format 
taxonomy_fixed_df <- as.data.frame(taxonomy_fixed_df) # force into a dataframe
row.names(taxonomy_fixed_df) <- taxonomy_fixed_df$Feature.ID # make first column into row names
taxonomy_fixed_df$Feature.ID <- NULL # remove the first column
taxonomy_matrix <- as.matrix(taxonomy_fixed_df) # convert to a matrix 

physeq_otu <- otu_table(otu_df, taxa_are_rows = T) # convert into phyloseq object
physeq_tax <- tax_table(taxonomy_matrix) # convert into phyloseq object
physeq_meta <- sample_data(metadata) # convert into phyloseq object

phylo_object <- phyloseq(physeq_otu, physeq_tax, physeq_meta) # merge into phyloseq object
phylo_object_tree <- merge_phyloseq(phylo_object, phylo_tree) # add tree into phyloseq object
```

### See phyloseq summary

``` r
phylo_object_tree
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 190829 taxa and 4212 samples ]
    ## sample_data() Sample Data:       [ 4212 samples by 53 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 190829 taxa by 24 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 190829 tips and 190350 internal nodes ]
