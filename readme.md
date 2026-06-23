# Bulk RNA-seq of Rheumatoid Arthritis 

## Description
Rheumatoid arthritis (RA), a chronic inflammatory autoimmune disorder, is characterised by persistent synovial inflammation, erosion of bones and cartilage, leading to joint destruction.
Objective: to identify differentially expressed genes (DEGs) in early RA patients and explore biological processes, cellular components and molecular pathways associated with early disease development.

## Data 
Acession number : GSE89408
SRA project : SRP092408
Organism : Homo sapiens
Comparison : Ealy RA versus normal

## Differential expression analysis was performed using DESeq2 

## Dependencies
library(DESeq2)
library(apeglm)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(pathview)
library(ReactomePA)
