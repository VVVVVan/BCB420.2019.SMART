# `BCB420.2019.SMART`

#### (**SMART** data annotation of human genes)

&nbsp;

###### [Fan Shen](https://orcid.org/0000-0001-8720-5874), Undergraduate of Bioinformatics and Computaional Biology, University of Toronto, Canada. &lt;van.shen@mail.utoronto.ca&gt;

----

**If any of this information is ambiguous, inaccurate, outdated, or incomplete, please check the [most recent version](https://github.com/VVVVVan/BCB420.2019.SMART) of the package on GitHub and [file an issue](https://github.com/VVVVVan/BCB420.2019.SMART/issues).**

----

<!-- TOCbelow -->
1. About this package:<br/>
&nbsp;&nbsp;&nbsp;&nbsp;1.1. Description<br/>
&nbsp;&nbsp;&nbsp;&nbsp;1.1.0 Some explanations<br/>
&nbsp;&nbsp;&nbsp;&nbsp;1.2 Package structure<br/>
2. SMART data<br/>
&nbsp;&nbsp;&nbsp;&nbsp;2.1. Data information<br/>
3. Data download and Cleanup<br/>
4. Mapping HGNC symbols to SMART domain IDs<br/>
&nbsp;&nbsp;&nbsp;&nbsp;4.0 Preparations<br/>
&nbsp;&nbsp;&nbsp;&nbsp;4.1 Domain IDs from database<br/>
&nbsp;&nbsp;&nbsp;&nbsp;4.2 Retrive data from biomaRt<br/>
&nbsp;&nbsp;&nbsp;&nbsp;4.3 Create ID-mapping tool<br/>
&nbsp;&nbsp;&nbsp;&nbsp;4.4 Outdated Symbol<br/>
&nbsp;&nbsp;&nbsp;&nbsp;4.5 Final Validation<br/>
5. Annotating gene sets with SMART Data<br/>
&nbsp;&nbsp;&nbsp;&nbsp;5.1 Get the data from file<br/>
&nbsp;&nbsp;&nbsp;&nbsp;5.2 Load smart2sym map and simple domain statistics<br/>
&nbsp;&nbsp;&nbsp;&nbsp;5.3 Load sym2smart map and simple gene statistics<br/>
6. Annotating example gene sets with SMART Data <br/>
7. Validate import process<br/>
8. Notes<br/>
9. Reference<br/>
10. Acknowledgements<br/>
<!-- TOCabove -->

----


# 1 About this package:

### 1.1 Description
The package describes how to get the domains data from [the SMART database](http://smart.embl-heidelberg.de), how to get the domains information from [HGNC symbols](https://www.genenames.org), how to annotate the example gene set, and it reports the data statistics and validation of import process.

&nbsp;

#### 1.1.0 Some explanations
In generally, we are not searching protein domains by domain ID. Given the fact that SMART contains protein domains data, this package is going to map the HGNC symbols to domain ID and then get the domain information (a reverse way compared to others in BCB420).

&nbsp;

### 1.2 Package structure
In this project:
```text
 --BCB420.2019.SMART/
   |__.gitignore
   |__.Rbuildignore
   |__BCB420.2019.SMART.Rproj
   |__DESCRIPTION
   |__inst/
      |__dev/
         |__rptTwee.R
         |__toBrowser.R
      |__extdata/
         |__results.tsv
         |__smart2sym.RData
         |__sym2smart.RData
      |__img/       # Store images
      |__scripts/
         |__scriptTemplate.R
   |__LICENSE
   |__man/
   |__NAMESPACE
   |__R/
      |__zzz.R
   |__README.md
   |__tests/
      |__testthat.R
      |__testthat/
         |__helper-functions.R
```

&nbsp;

----

# 2 SMART data
SMART is a database of protein domains. SMART domains are found in non-redundant protein database system. For each gene, name of confidently predicted domains are shown and are presented in different shapes for each domain on sequences. All SMART data is available under a CC-BY 4.0 license.

This document describes work with [SMART version 8](http://smart.embl-heidelberg.de) [(Letunic & Bork. 2018)](https://doi.org/10.1093/nar/gkx922).

&nbsp;

### 2.1 Data information
SMART domains are detected by multiple sequence alignments of representative family members and these alignments are optimised manually by contructing hidden Markov model (HMM). [(Schultz _et al._ 2000)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC102444/)

Now SMART is a growing collection for shuffled extracellular domains and provides powerful web-based interface with visualization tools. (Letunic & Bork. 2018)

&nbsp;

----

# 3 Data download and Cleanup

To download the source data from [SMART](http://smart.embl-heidelberg.de):

1. Navigate to [SMART Website](http://smart.embl-heidelberg.de)
2. At lower left corner, follow the link to [download](http://smart.embl-heidelberg.de/smart/descriptions.pl) the domain description:

* `descriptions.pl` (344 kb) domain descriptions

3. Put the file into a sister directory of your working directory which is called data. (It should be reachable with file.path("..", "data"))

&nbsp;

----

# 4 Mapping HGNC symbols to SMART domain IDs
SMART domains have domain ID. They cannot be map to HGNC symbol directly, but could do it in opposite way. Also, one gene (HGNC symbol) may contains lots of different domains (SMART IDs) or no annotated SMART domains. To provide the best possible interpretation, we need to build a map of HGNC symbols to domain IDs. Please aware that there may exist HGNC symbols with lots of domain IDs or no domain IDs. **The usability of the dataset for annotation depends on the quality of this mapping.**

&nbsp;

### 4.0 Preparations
Install required packages:
1. `biomaRt` biomaRt is a Bioconductor package that implements the RESTful API of biomart, the annotation framwork for model organism genomes at the EBI. It is a Bioconductor package, and as such it needs to be loaded via the BiocManager (Steipe, 2019).

```R
if (! requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (! requireNamespace("biomaRt", quietly = TRUE)) {
  BiocManager::install("biomaRt")
}
```

2. `httr` httr is a package that provides a wrapper for the curl package, customised to the demands of modern web APIs

```R
if (! requireNamespace("httr", quietly = TRUE)) {
  install.packages("httr")
}
```

&nbsp;

Get HGNC data from GitHub (Steipe, 2019):

```R
myURL <- paste0("https://github.com/hyginn/",
                "BCB420-2019-resources/blob/master/HGNC.RData?raw=true")
load(url(myURL))  # loads HGNC data frame
```

&nbsp;

### 4.1 Domain IDs from database
We can check how many domain IDs we have from database.

&nbsp;

```R
# Read the domain description: a file with domain name, ID, definition and
# description.
tmp <- read.delim(file.path("../data", "descriptions.pl"), 
                         header=TRUE, 
                         sep="\t", 
                         quote="", 
                         stringsAsFactors = FALSE)

# Remove the header of the file
tmp <- tmp[-1,]

# Define the row names
colnames(tmp) <- c("name", "ID", "Definition", "Description")
rownames(tmp) <- c(seq(1,nrow(tmp)))

# The data looks like:
head(tmp)
#           name      ID                                          Definition
#1        14_3_3 SM00101                                   14-3-3 homologues
#2        35EXOc SM00474                                   3'-5' exonuclease
#3          4.1m SM00294         putative band 4.1 homologues' binding motif
#4        53EXOc SM00475                                   5'-3' exonuclease
#5          6PGD SM01350 6-phosphogluconate dehydrogenase, C-terminal domain
#6 7TM_GPCR_Srsx SM01381         Serpentine type 7TM GPCR chemoreceptor Srsx

#         Description
#1        14-3-3 homologues mediates signal transduction by binding to phosphoserine-containing proteins. They are involved in growth factor signalling and also interact with MEK kinases.
#2        3\\' -5' exonuclease proofreading domain present in DNA polymerase I, Werner syndrome helicase, RNase D and other enzymes
#3        
#4        
#5        This family represents the C-terminal all-alpha domain of 6-phosphogluconate dehydrogenase. The domain contains two structural repeats of 5 helices each.
#6 Chemoreception is mediated in Caenorhabditis elegans by members of the seven-transmembrane G-protein-coupled receptor class (7TM GPCRs) of proteins which are of the serpentine type (PMID:7585938). Srsx is a solo family amongst the superfamilies of chemoreceptors. Chemoperception is one of the central senses of soil nematodes like C. elegans which are otherwise blind and deaf (PMID:18050473)

# Note: Some of the domains do not have proper Description and Definition.

# Do all smybol start with "SM"?
all(grepl("^SM", tmp$ID)) # TRUE

# How many unique IDs do we have?
uSMART <- unique(tmp$ID) # 1299, all are unique
```

&nbsp;

### 4.2 Retrive data from biomaRt
To get the mapping, we use biomaRt to retrive data, including HGNC symbol, SMART domain IDs, APPRIS annotation, Protein stable IDs. HGNC symbol and SMART domain IDs are mapping to each other; APPRIS annotations tell the principal trancript; Protein stable IDs gives identity of a gene.

&nbsp;

#### 4.2.1 Get Data
```R
myMart <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")

# Example code to get the list of filters, attributes and look for correct ones
filters <- listFilters(myMart)
filters[grep("SMART", filters$description), ]
#          name                        description
#356 with_smart           With SMART domains ID(s)
#374      smart SMART domains ID(s) [e.g. SM00002]

attributes <- listAttributes(myMart)
attributes[grep("HGNC symbol", attributes$description, ignore.case=TRUE), ]
#          name description         page
#59 hgnc_symbol HGNC symbol feature_page

tmp <- biomaRt::getBM(filters = "with_smart",
                      attributes = c("hgnc_symbol",
                                     "smart",
                                     "transcript_appris",
                                     "ensembl_peptide_id"),
                      values = TRUE,
                      mart = myMart)
head(tmp)
#  hgnc_symbol   smart transcript_appris ensembl_peptide_id
#1   ARHGAP11A SM00324                      ENSP00000479117
#2       NLRP7 SM01289                      ENSP00000482551
#3       NLRP7 SM00368                      ENSP00000482551
#4       NLRP7 SM01289                      ENSP00000481809
#5       NLRP7 SM00368                      ENSP00000481809
#6             SM00008                      ENSP00000484387

nrow(tmp) # 59215
```

&nbsp;

The number of HGNC symbol to domain IDs is much larger than the number of domains in database. There are three possible reasons:
&nbsp;

&nbsp;

**(1)** Not all transcript are principal transcript. Lots of them have no APPRIS annotation, return `""`. Some rows are not start with `"principal"` in APPRIS anotation. We just focus on "principal" since they are "Transcript(s) expected to code for the main functional isoform based solely on the core modules in the APPRIS" (from [APPRIS help page](http://appris.bioinfo.cnio.es/#/help/scores)).

```R
nrow(tmp[tmp$transcript_appris != "",]) # 31750 rows have APPRIS annotation
nrow(tmp[grepl("^principal", tmp$transcript_appris),]) # 21921 rows is principal1
```
&nbsp;

**(2)** There also occurs some situation where no hgnc symbol is assigned to pepetide id, such that the `hgnc_symbol` returns `""`.

```R
nrow(tmp[tmp$hgnc_symbol == "",]) # 516
```

&nbsp;

**(3)** There are duplicated hgnc symbols. 
a) Some of hgnc symbols have more than one domains or one domain for many times, which is resonable.
b) Some of hgnc symbols have more than one pepeitde IDs.

```R
sum(duplicated(tmp$hgnc_symbol)) # 49277 duplications

uni <- unique(tmp[c('hgnc_symbol', 'smart')])
length(unique(uni$hgnc_symbol[duplicated(uni$hgnc_symbol)]))
# 3950 hgnc symbols have different domains

dupe <- tmp[,c('hgnc_symbol','smart')] # select columns to check duplicates
dupe <- dupe[dupe$hgnc_symbol != "",]
length(unique(dupe$hgnc_symbol[duplicated(dupe) | duplicated(dupe, fromLast=TRUE)]))
# 7500 hgnc symbols have same domain for many times

uni <- unique(tmp[c('hgnc_symbol', 'ensembl_peptide_id')])
length(unique(uni$hgnc_symbol[duplicated(uni$hgnc_symbol)]))
# 7505 hgnc symbols have differnt peptide IDs

# Clean up the environment if necessary
rm(uni)
rm(dupe)
```

&nbsp;


#### 4.2.2 Cleanup and validation of data
We need to clean up the data from previous section.

&nbsp;

**(1)** Remove all rows that are not `principal1` for APPRIS annotation.

```R
tmp <- tmp[grepl("^principal", tmp$transcript_appris),]

# Comfirm
all(grepl("^principal", tmp$transcript_appris)) # TRUE
```

&nbsp;

**(2)** Remove all empty hgnc symbols, that return `""`.

```R
tmp <- tmp[tmp$hgnc_symbol != "",]

# Comfiem
all(!("" %in% tmp$hgnc_symbol)) # TRUE
```

**(3)** Deal with it when constructing the ID-mapping tool

&nbsp;

### 4.3 Create ID-mapping tool (store in list with different length of hgnc symbols)
We need to create a ID-mapping tool. We create two different mapping. One from SMART ID to symbol because each domain could be in many different genes. One from symbol to SMART IS becasue each HGNC symbols may have different number of domains. We store the data in a list. For each HGNC symbol, we just store the unique domain rather than duplicated domain. This will also remove HGNC symbol with same domain but different pepetide IDs.

&nbsp;

```R
# smart2sym and sym2smart list
smart2sym <- list()
sym2smart <- list()
for (i in seq_along(tmp$hgnc_symbol)) {
  # smart to sym
  if (is.null(smart2sym[[tmp$smart[i]]])) { # not present in mapping
    smart2sym[[tmp$smart[i]]] <- c(tmp$hgnc_symbol[i])
  } else { # the symbol is in the mapping
    # Record the hgnc_symbol if it is not in the mapping yet
    if  (! (tmp$hgnc_symbol[i] %in% smart2sym[[tmp$smart[i]]])) { 
      smart2sym[[tmp$smart[i]]] <- 
                c(smart2sym[[tmp$smart[i]]], tmp$hgnc_symbol[i])
    }
  }
  # sym to smart
  if (is.null(sym2smart[[tmp$hgnc_symbol[i]]])) { # not present in mapping
    sym2smart[[tmp$hgnc_symbol[i]]] <- c(tmp$smart[i])
  } else { # the smart is in the mapping
    # Record the smart if it is not in the mapping yet
    if  (! (tmp$smart[i] %in% sym2smart[[tmp$hgnc_symbol[i]]])) { 
      sym2smart[[tmp$hgnc_symbol[i]]] <- 
                c(sym2smart[[tmp$hgnc_symbol[i]]], tmp$smart[i])
    }
  }
}

head(smart2sym)
#$SM00409
#  [1] "MYBPC2"    "F11R"      "CEACAM5"   "PAPLN"     "TRIO"      "LINGO4"    "LRIT3"    
#  [8] "CEACAM21"  "A1BG"      "SIGLEC11"  "HAPLN4"    "MXRA5"     "CD300C"    "SEMA4D" 
# .... more not shown
#
#$SM00408
#  [1] "MYBPC2"   "F11R"     "CEACAM5"  "PAPLN"    "TRIO"     "LINGO4"   "LRIT3"   
#  [8] "CEACAM21" "A1BG"     "SIGLEC11" "HAPLN4"   "MXRA5"    "SEMA4D"   "MYOM2"     
# .... more not shown
#
#$SM00060
#  [1] "MYBPC2"   "CMYA5"    "PTPRB"    "IL31RA"   "NDNF"     "MYOM2"    "AXL"     
#  [8] "TYRO3"    "PRLR"     "EPHA7"    "FSD1"     "IGSF22"   "PTPRT"    "IGF1R"  
# .... more not shown
#
#$SM00054
#  [1] "RCN3"     "RASEF"    "GPD2"     "RHOT1"    "FKBP10"   "CRNN"     "S100A6"  
#  [8] "S100A11"  "EFCAB2"   "CABP1"    "PLCH1"    "ACTN1"    "TESC"     "SPTAN1"    
# .... more not shown
#
#$SM00222
# [1] "ARFGEF3" "IQSEC3"  "PSD2"    "CYTH3"   "ARFGEF1" "CYTH1"   "CYTH4"   "IQSEC2" 
# [9] "CYTH2"   "PSD4"    "IQSEC1"  "FBXO8"   "PSD3"    "ARFGEF2" "GBF1"    "PSD"     
#
#$SM00248
#  [1] "TRPV5"           "SOWAHA"          "ANKRD33"         "TRPV4"          
#  [5] "ABTB2"           "MTPN"            "ANKRD37"         "ANKRD52" 
# .... more not shown

head(sym2smart)
#$MYBPC2
#[1] "SM00409" "SM00408" "SM00060"
#
#$RCN3
#[1] "SM00054"
#
#$ARFGEF3
#[1] "SM00222"
#
#$TRPV5
#[1] "SM00248"
#
#$OR2M4
#[1] "SM01381"
#
#$ZNF574
#[1] "SM00355"
```

&nbsp;

#### 4.3.1 Additional symbols
Now we can check if all domain IDs in my database is in this list. If not, check if we could get the data from biomaRt by the SMART domain IDs

&nbsp;

```R
sel <- (! (uSMART %in% names(smart2sym)))

test <- biomaRt::getBM(filters = "smart",
                      attributes = c("hgnc_symbol",
                                     "smart",
                                     "transcript_appris",
                                     "ensembl_peptide_id"),
                      values = uSMART[sel],
                      mart = myMart)
                      
# Check if they are start with `principal`
length(test$smart[grepl("^principal", test$transcript_appris)]) # 0
# No additional symbols

sum(sel) # 302 domain ID do not map to any HGNC symbols.

# Store the symbol as NA in smart2sym
for (item in uSMART[sel]) {
  smart2sym[[item]] <- NA
}

# Check if all domain ID from database is in the map now.
all((uSMART %in% names(smart2sym))) # TRUE

# Check if all ID in map could be found in database.
all(names(smart2sym) %in% uSMART) # FALSE 
sum(! names(smart2sym) %in% uSMART) # 23 not found in database.

# Store these domain ID as NotMapped in sym2smart
sym2smart[["NotMapped"]] <- uSMART[sel]
```

&nbsp;

Thus, there are no additional symbol for my list and 382 of domain IDs from my database do not have any HGNC symbols. TODO: Ask for help if there is any other way to map the domain IDs to HGNC symbols.

&nbsp;

### 4.4 Outdated Symbol
We now validate our mapping tool by the HGNC reference from Prof Steipe, loaded before.

&nbsp;

```R
# The code is based on Prof Steipe's BCB420.2019.STRING outdataed symbol section
# It is hard to get the symbols from list, so get the outdated symbol from tmp
sel <- (  !  (tmp$hgnc_symbol %in% HGNC$sym))
length(       tmp$hgnc_symbol[sel]) # 131 unkown
length(unique(tmp$hgnc_symbol[sel])) # 76 are unique

unkSym <- data.frame(unk = unique(tmp$hgnc_symbol[sel]),
                     new = NA,
                     stringsAsFactors = FALSE)

# We will leave those that not map to any HGNC symbol in resource to NA

# grep() for the presence of the symbols in either HGNC$prev or
# HGNC$synonym. If either is found, that symbol replaces NA in unkSym$new
for (i in seq_len(nrow(unkSym))) {
  iPrev <- grep(unkSym$unk[i], HGNC$prev)[1] # take No. 1 if there are several
  if (length(iPrev) == 1) {
    unkSym$new[i] <- HGNC$sym[iPrev]
  } else {
    iSynonym <- which(grep(unkSym$unk[i], HGNC$synonym))[1]
    if (length(iSynonym) == 1) {
      unkSym$new[i] <- HGNC$sym[iSynonym]
    }
  }
}
  
sum(! is.na(unkSym$new)) # 7 symbol has map to symbol in HGNC resource

# We add the contents of unkSym$new back into smart2sym and sym2smart. This
# way, the newly mapped symbols are updated, and the old symbols that did not
# map are delete in smart2sym and set to NA in sym2smart.

# smart2sym
for (i in seq_along(unkSym$unk)) {
  IDs <- tmp$smart[tmp$hgnc_symbol == unkSym$unk[i]] # Get the ID
  for (ID in IDs) {
    if (length(smart2sym[[ID]]) == 1) { # Replace the only item
      smart2sym[[ID]] <- c(unkSym$new[i])
    } else {
      if (! is.na(match(unkSym$unk[i], smart2sym[[ID]]))) { # if there is a match
        smart2sym[[ID]] <- smart2sym[[ID]][-match(unkSym$unk[i], smart2sym[[ID]])]
        # If the item is not NA, add to list
        # If the item is NA and length of list is 0, add NA to the list
        if (! is.na(unkSym$new[i])) {
          smart2sym[[ID]] <- c(smart2sym[[ID]], unkSym$new[i])
        } else if (length(smart2sym[[ID]]) == 0) {
          smart2sym[[ID]] <- NA
        }
      }
    }
  }
}

# sym2smart
for (i in seq_along(unkSym$unk)) {
   names(sym2smart)[match(unkSym$unk[i], names(sym2smart))] <- unkSym$new[i]
}
```

&nbsp;

### 4.5 Final validation
Finally, we validate my mapping tool.

&nbsp;

```R
# 1. Are all domain ID in the map?
# Check if all domain ID from database is in the map now.
all((uSMART %in% names(smart2sym))) # TRUE

# Check if all ID in map could be found in database.
all(names(smart2sym) %in% uSMART) # FALSE 
sum(! names(smart2sym) %in% uSMART) # 23 not found in database.

# 2. How many symbol are there?
sum(! (is.na(names(sym2smart)) | names(sym2smart) == "NotMapped" )) # 9821

# 4. How many percent of symbol is mapped?
sum(! (is.na(names(sym2smart)) | names(sym2smart) == "NotMapped" )) * 100 /
                                (length(names(sym2smart)) -1) # 99.30%

# 5. Are all symbols now in reference table?
all(names(sym2smart)[! (is.na(names(sym2smart)) | 
                     names(sym2smart) == "NotMapped")] %in% HGNC$sym) # TRUE

# Done.
# Save the map:
save(smart2sym, file = file.path("inst", "extdata", "smart2sym.RData"))
save(sym2smart, file = file.path("inst", "extdata", "sym2smart.RData"))

# From an RStudio project, the file can be loaded with:
load(file = file.path("inst", "extdata", "smart2sym.RData"))
load(file = file.path("inst", "extdata", "sym2smart.RData"))
```

&nbsp;

----

# 5 Annotating gene sets with SMART Data
Given our mapping tool, we can now annotate gene sets with SMART data. First, we analyze the entire STRING data. 

&nbsp;

#### 5.1 Get the data from file
```R
# Read the domain description: a file with domain name, ID, definition and
# description.
tmp <- read.delim(file.path("../data", "descriptions.pl"), 
                         header=TRUE, 
                         sep="\t", 
                         quote="", 
                         stringsAsFactors = FALSE)

# Remove the header of the file
tmp <- tmp[-1,]

# Define the row names
colnames(tmp) <- c("name", "ID", "Definition", "Description")
rownames(tmp) <- c(seq(1,nrow(tmp)))

storedDefinition <- tmp
save(storedDefinition, file = file.path("inst", "extdata", "storedDefinition.RData"))

# There are some domains do not have definition or description
sum(tmp$Definition == "") # 202
sum(tmp$Description == "") # 337
sum(tmp$Description == "" & tmp$Definition == "") # 31

length(tmp$ID)
a <- data.frame(LackDefinition = c(sum(tmp$Definition == ""), 
                round(sum(tmp$Definition == "") / length(tmp$ID) * 100, 2)),
               LackDescription = c(sum(tmp$Description == ""),
                round(sum(tmp$Description == "") / length(tmp$ID) * 100, 2)),
               LackBoth = c(sum(tmp$Description == "" & tmp$Definition == ""),
                round(sum(tmp$Description == "" & tmp$Definition == "") / 
                       length(tmp$ID), 2)),
               stringsAsFactors = FALSE)
rownames(a) <- c("Counts", "Percentage %")

library(gridExtra)
library(grid)
grid.newpage()
grid.table(a)
```
![](./inst/img/SMART.jpeg?sanitize=true "SMART")
&nbsp;

#### 5.2 Load smart2sym map and simple domain statistics

```R
# Load the map
load(file = file.path("inst", "extdata", "smart2sym.RData"))

# First, how many genes does one domain located.
lengths <- as.integer(unlist(lapply(smart2sym, function(x) length(x))))
max(lengths) # 527
min(lengths) # 1
hist(lengths, 
     xlim = c(min(lengths), max(lengths)),
     main = "SMART domain mapped gene numbers",
     col = rainbow(4, s= 0.5),
     xlab = "Number of genes",
     ylab = "Counts")
abline(v = mean(lengths), lwd = 0.5, col = "black")
```
![](./inst/img/SMART_domain.jpeg?sanitize=true "SMART domain")

&nbsp;

There are too much small numbers in the map, so it is hard to tell what happened in right part of the plot. We could zoom in to know more.

&nbsp;

```R
a <- hist(lengths, 
     xlim = c(min(lengths) + 1, max(lengths) + 50),
     ylim = c(0, 100),
     breaks = c(seq(0,max(lengths) + 50, 20)),
     main = "SMART domain mapped gene numbers",
     col = heat.colors(11),
     xlab = "Number of genes",
     ylab = "Counts",
     labels = TRUE)
abline(v = mean(lengths), lwd = 0.5, col = "black")
```
![](./inst/img/SMART_zoom.jpeg?sanitize=true "SMART domain zoom in")

&nbsp;

We can tell that most domains are mapped to 0-20 genes, while the number of mapping decrease dramatically after 20 genes.

&nbsp;

#### 5.3 Load sym2smart map and simple gene statistics
```R
# Load the map
load(file = file.path("inst", "extdata", "sym2smart.RData"))

# Clean up NA and NotMapped in sym2smart
symDomains <- sym2smart
indexes <- which(is.na(names(symDomains)))
for(i in indexes) {
  symDomains[[i]] <- NULL
}
symDomains[["NotMapped"]] <- NULL

# Save it
save(symDomains, file = file.path("inst", "extdata", "symDomains.RData"))

# Coverage of human protein genes
length(names(symDomains)) * 100 / sum(HGNC$type == "protein") # 51.10%

# Then, how many domains does one gene has.
lengths <- as.integer(unlist(lapply(symDomains, function(x) length(x))))
largest <- max(lengths) # 10
smallest <- min(lengths) # 1

# Plot the hist
a <- hist(lengths, 
           xlim = c(smallest, largest),
           breaks = c(seq(1, largest,1)),
           main = "SMART gene mapped domain numbers",
           col = heat.colors(10),
           xlab = "Number of domains",
           ylab = "Counts",
           labels = TRUE)
abline(v = mean(lengths[lengths != max(lengths)]), lwd = 1, col = "black")
```
![](./inst/img/SMART_gene.jpeg?sanitize=true "SMART gene")

&nbsp;
We can tell that most genes have one domain. Very few genes have more than one domain from the database.
&nbsp;

----

# 6 Annotating example gene sets with SMART Data
Next, we use data to analyze the domains for our example gene set and store it.

```R
# The code is based on Prof Steipe's BCB420.2019.STRING project

# The specification of the sample set is copy-paste from the 
# BCB420 resources project.

xSet <- c("AMBRA1", "ATG14", "ATP2A1", "ATP2A2", "ATP2A3", "BECN1", "BECN2",
          "BIRC6", "BLOC1S1", "BLOC1S2", "BORCS5", "BORCS6", "BORCS7",
          "BORCS8", "CACNA1A", "CALCOCO2", "CTTN", "DCTN1", "EPG5", "GABARAP",
          "GABARAPL1", "GABARAPL2", "HDAC6", "HSPB8", "INPP5E", "IRGM",
          "KXD1", "LAMP1", "LAMP2", "LAMP3", "LAMP5", "MAP1LC3A", "MAP1LC3B",
          "MAP1LC3C", "MGRN1", "MYO1C", "MYO6", "NAPA", "NSF", "OPTN",
          "OSBPL1A", "PI4K2A", "PIK3C3", "PLEKHM1", "PSEN1", "RAB20", "RAB21",
          "RAB29", "RAB34", "RAB39A", "RAB7A", "RAB7B", "RPTOR", "RUBCN",
          "RUBCNL", "SNAP29", "SNAP47", "SNAPIN", "SPG11", "STX17", "STX6",
          "SYT7", "TARDBP", "TFEB", "TGM2", "TIFA", "TMEM175", "TOM1",
          "TPCN1", "TPCN2", "TPPP", "TXNIP", "UVRAG", "VAMP3", "VAMP7",
          "VAMP8", "VAPA", "VPS11", "VPS16", "VPS18", "VPS33A", "VPS39",
          "VPS41", "VTI1B", "YKT6")


# Which gene has no domain annotated?
x <- which( ! (xSet %in% names(symDomains)))
length(x) # 41 genes do not have domains
cat(sprintf("\t%s\t(%s, %s)\n", HGNC[xSet[x], "sym"], HGNC[xSet[x], "name"],HGNC[xSet[x], "UniProtID"]))
#	  ATG14	(autophagy related 14, Q6ZNE5)
# 	BECN1	(beclin 1, Q14457)
# 	BECN2	(beclin 2, A8MW95)
# 	BLOC1S1	(biogenesis of lysosomal organelles complex 1 subunit 1, P78537)
# 	BLOC1S2	(biogenesis of lysosomal organelles complex 1 subunit 2, Q6QNY1)
# 	BORCS5	(BLOC-1 related complex subunit 5, Q969J3)
# 	BORCS6	(BLOC-1 related complex subunit 6, Q96GS4)
# 	BORCS7	(BLOC-1 related complex subunit 7, Q96B45)
# 	BORCS8	(BLOC-1 related complex subunit 8, Q96FH0)
# 	CALCOCO2	(calcium binding and coiled-coil domain 2, Q13137)
# 	EPG5	(ectopic P-granules autophagy protein 5 homolog, Q9HCE0)
# 	GABARAP	(GABA type A receptor-associated protein, O95166)
# 	GABARAPL1	(GABA type A receptor associated protein like 1, Q9H0R8)
# 	GABARAPL2	(GABA type A receptor associated protein like 2, P60520)
# 	HSPB8	(heat shock protein family B (small) member 8, Q9UJY1)
# 	IRGM	(immunity related GTPase M, A1A4Y4)
# 	KXD1	(KxDL motif containing 1, Q9BQD3)
# 	LAMP1	(lysosomal associated membrane protein 1, P11279)
# 	LAMP2	(lysosomal associated membrane protein 2, P13473)
# 	LAMP3	(lysosomal associated membrane protein 3, Q9UQV4)
# 	LAMP5	(lysosomal associated membrane protein family member 5, Q9UJQ1)
# 	MAP1LC3A	(microtubule associated protein 1 light chain 3 alpha, Q9H492)
# 	MAP1LC3B	(microtubule associated protein 1 light chain 3 beta, Q9GZQ8)
# 	MAP1LC3C	(microtubule associated protein 1 light chain 3 gamma, Q9BXW4)
# 	NAPA	(NSF attachment protein alpha, P54920)
# 	OPTN	(optineurin, Q96CV9)
# 	PI4K2A	(phosphatidylinositol 4-kinase type 2 alpha, Q9BTU6)
# 	SNAP47	(synaptosome associated protein 47, Q5SQN1)
# 	SNAPIN	(SNAP associated protein, O95295)
# 	SPG11	(SPG11, spatacsin vesicle trafficking associated, Q96JI7)
# 	TIFA	(TRAF interacting protein with forkhead associated domain, Q96CG3)
# 	TMEM175	(transmembrane protein 175, Q9BSA9)
# 	TPCN1	(two pore segment channel 1, Q9ULQ1)
# 	TPCN2	(two pore segment channel 2, Q8NHX9)
# 	TPPP	(tubulin polymerization promoting protein, O94811)
# 	VAMP3	(vesicle associated membrane protein 3, Q15836)
# 	VAMP8	(vesicle associated membrane protein 8, Q9BV40)
# 	VAPA	(VAMP associated protein A, Q9P0L0)
# 	VPS16	(VPS16, CORVET/HOPS core subunit, Q9H269)
# 	VPS18	(VPS18, CORVET/HOPS core subunit, Q9P253)
# 	VPS33A	(VPS33A, CORVET/HOPS core subunit, Q96AX1)

# Most type of not annoated proteins is transmembrane protein. 
# Check some of the genes in SMART, http://smart.embl-heidelberg.de/#.
# They don't have domain annotated in SMART.
# It seems like SMART do not assign ID/name to transmembrane domains.


# For our annotation, we just store the domain information for those have 
# domain in example set:
library(httr)
sel <- which(names(symDomains) %in% xSet)

store <- list()

load(file = file.path("inst", "extdata", "storedDefinition.RData"))

for (sym in sel) {
  domains <- symDomains[[sym]]
  for (domain in domains) {
    store[[names(symDomains)[sym]]] <- c(store[[names(symDomains)[sym]]],
                 storedDefinition$Definition[storedDefinition$ID == domain],
                 storedDefinition$Description[storedDefinition$ID == domain])
  }
}

domainsOut <- data.frame(symbol = character(),
                         Definition = character(),
                         Description = character(),
                         stringsAsFactors = FALSE)
for (i in seq_along(store)) {
  for (j in seq_along(store[[i]])) {
    if (j%%2 == 1) {
      add <- data.frame(symbol = names(store)[i],
                        Definition = store[[i]][j],
                        Description = store[[i]][j+1],
                        stringsAsFactors = FALSE)
      domainsOut <- rbind(domainsOut,add)
    }
  }
}


writeLines(c("symobl\tDefinition\tDescription",
             sprintf("%s\t%s\t%s\t", domainsOut[,1], domainsOut[,2], domainsOut[,3])),
           con = "inst/extdata/domains.tsv")
           
# The data set can be read back in again (in an RStudio session) with
myXset <- read.delim(file.path("inst", "extdata", "domains.tsv"),
                     stringsAsFactors = FALSE)

# From an installed package, the command would be:
myXset <- read.delim(system.file("extdata",
                                  "domains.tsv",
                                  package = "BCB420.2019.SMART"),
                     stringsAsFactors = FALSE)
```
&nbsp;

----

# 7 Validate import process
To validate the import process, we get the domain, symbol and uniprot ID from HGNC and go to SMART website to check if the gene has the domain.

```R
myMart <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")

# Example code to get the list of filters, attributes and look for correct ones
filters <- biomaRt::listFilters(myMart)
filters[grep("SMART", filters$description), ]

attributes <- biomaRt::listAttributes(myMart)
attributes[grep("sequence", attributes$page, ignore.case=TRUE), ]

tmp <- biomaRt::getBM(filters = "with_smart",
                      attributes = c("hgnc_symbol",
                                     "smart",
                                     "smart_start",
                                     "smart_end",
                                     "transcript_appris",
                                     "ensembl_peptide_id"),
                      values = TRUE,
                      mart = myMart)
                      
# Clean up the data as before
tmp <- tmp[grepl("^principal", tmp$transcript_appris),]
tmp <- tmp[tmp$hgnc_symbol != "",]

set.seed(410420)
testSet <- sample(tmp$hgnc_symbol, 10)
for (i in seq_along(testSet)) {
  if (sum(tmp$hgnc_symbol == testSet[i]) >20) {
    testSet[i] <- NA
  }
}
for(gene in testSet){
  if (! is.na(gene)) {
    print(unique(tmp[tmp$hgnc_symbol == gene,c(2,3,4,6)]))
  }
}
```
&nbsp;

The output is:

```text
        smart smart_start smart_end ensembl_peptide_id
57961 SM00112         376       458    ENSP00000231484
57962 SM00112         157       242    ENSP00000231484
57963 SM00112         621       704    ENSP00000231484
57964 SM00112         266       350    ENSP00000231484
57965 SM00112         482       563    ENSP00000231484
57966 SM00112          36       133    ENSP00000231484
         smart smart_start smart_end ensembl_peptide_id
122570 SM00295        1795      2012    ENSP00000386461
122571 SM00295        1193      1412    ENSP00000386461
122572 SM00139        1644      1793    ENSP00000386461
122573 SM00139         989      1192    ENSP00000386461
122574 SM00015         785       807    ENSP00000386461
122575 SM00015         762       784    ENSP00000386461
122576 SM00015         854       876    ENSP00000386461
122577 SM00015         831       853    ENSP00000386461
122578 SM00326        1504      1566    ENSP00000386461
122579 SM00242          59       761    ENSP00000386461
122584 SM00295        1795      2012    ENSP00000415090
122585 SM00295        1193      1412    ENSP00000415090
122586 SM00139        1644      1793    ENSP00000415090
122587 SM00139         989      1192    ENSP00000415090
122588 SM00015         785       807    ENSP00000415090
122589 SM00015         762       784    ENSP00000415090
122590 SM00015         854       876    ENSP00000415090
122591 SM00015         831       853    ENSP00000415090
122592 SM00326        1504      1566    ENSP00000415090
122593 SM00242          59       761    ENSP00000415090
        smart smart_start smart_end ensembl_peptide_id
57435 SM00054         819       847    ENSP00000317997
57436 SM00054         855       883    ENSP00000317997
57437 SM00054         819       847    ENSP00000373689
57438 SM00054         855       883    ENSP00000373689
        smart smart_start smart_end ensembl_peptide_id
48161 SM00349           4        64    ENSP00000395733
48162 SM00355         341       363    ENSP00000395733
48163 SM00355         313       335    ENSP00000395733
48164 SM00355         369       391    ENSP00000395733
48165 SM00355         509       531    ENSP00000395733
48166 SM00355         285       307    ENSP00000395733
48167 SM00355         425       447    ENSP00000395733
48168 SM00355         397       419    ENSP00000395733
48169 SM00355         201       223    ENSP00000395733
48170 SM00355         453       475    ENSP00000395733
48171 SM00355         481       503    ENSP00000395733
48172 SM00355         257       279    ENSP00000395733
48173 SM00355         173       195    ENSP00000395733
48174 SM00355         229       251    ENSP00000395733
```
&nbsp;

Then go to [SMART](http://smart.embl-heidelberg.de/) check if the gene with uniprot ID have the same domain description and start and end.

For example search for first gene with peptide id `ENSP00000231484`:
![](./inst/img/exampleSmart.png?sanitize=true "exammple smart screenshot")
We can tell that the start and end of the domains are same as what we get from biomaRt and in the page, we can tell the ID of the domain name is same as `SM00112`.

&nbsp;

I also check other symbols, the start and end of the domains are same. However, generally, SMART website will produce more domains including `low complexity`.

&nbsp;

To further check if the sequence are correct. I get the transcript sequence from biomaRt and check it to SMART website.

```R
#        smart smart_start smart_end ensembl_peptide_id
#57961 SM00112         376       458    ENSP00000231484
#57962 SM00112         157       242    ENSP00000231484
#57963 SM00112         621       704    ENSP00000231484
#57964 SM00112         266       350    ENSP00000231484
#57965 SM00112         482       563    ENSP00000231484
#57966 SM00112          36       133    ENSP00000231484

# Use the above symbol as example, we can tell the peptide id, the start and end
test <- biomaRt::getBM(filters = "ensembl_peptide_id",
                      attributes = c("peptide",
                                     "ensembl_peptide_id"),
                      values = "ENSP00000231484",
                      mart = myMart)
                      
seq <- test$peptide
substring(seq, 36, 133)
# [1] "VSEEVPSGTVIGKLSQELGREERRRQAGAAFQVLQLPQALPIQVDSEEGLLSTGRRLDREQLCRQWDPCLVSFDVLATGDLALIHVEIQVLDINDHQP"

# Then we go to the smart website page again and serch for the peptide id
# The domain seq from 36 to 133:
# VSEEVPSGTVIGKLSQELGREERRRQAGAAFQVLQLPQALPIQVDSEEGLLSTGRRLDRE
#QLCRQWDPCLVSFDVLATGDLALIHVEIQVLDINDHQP

# Since the seq is short, we can tell by eye that they are same.

# More test:
substring(seq, 157, 242)
#[1] "DRALDPDTGPNTLHTYTLSPSEHFALDVIVGPDETKHAELIVVKELDREIHSFFDLVLTAYDNGNPPKSGTSLVKVNVLDSNDNSP"

# DRALDPDTGPNTLHTYTLSPSEHFALDVIVGPDETKHAELIVVKELDREIHSFFDLVLTA
#YDNGNPPKSGTSLVKVNVLDSNDNSP


substring(seq, 376, 458)
#[1] "VMADDLDSGHNGLVHCWLSQELGHFRLKRTNGNTYMLLTNATLDREQWPKYTLTLLAQDQGLQPLSAKKQLSIQISDINDNAP"

# VMADDLDSGHNGLVHCWLSQELGHFRLKRTNGNTYMLLTNATLDREQWPKYTLTLLAQDQ
#GLQPLSAKKQLSIQISDINDNAP

```
THE END of Validation!!!
&nbsp;

----

# 8 Notes

Some useful keyboard shortcuts for package authoring:

Build and Reload Package: `Cmd + Shift + B`<br/>
Update Documentation: `Cmd + Shift + D` or `devtools::document()`<br/>
Test Package: `Cmd + Shift + T`<br/>
Check Package: `Cmd + Shift + E` or `devtools::check()`<br/>

&nbsp;
----

# 9 Reference
This package is written according to the example provided by [Prof Steipe, BCB420.2019.STRING](https://github.com/hyginn/BCB420.2019.STRING).

The HGNC symbol data and example genes are get from [Prof Steipe, BCB420.2019.Resource](https://github.com/hyginn/BCB420-2019-resources).

* Letunic, I., & Bork, P. (2018). 20 years of the SMART protein domain annotation resource. [_Nucleic Acids Research_, D1, D493–D496](https://doi.org/10.1093/nar/gkx922).

* Schultz, J., Copley, R. R., Doerks, T., Ponting, C. P., & Bork, P. (2000). SMART: a web-based tool for the study of genetically mobile domains. [_Nucleic acids research_, 28(1), 231-4](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC102444/).

&nbsp;
----

# 10 Acknowledgements

&nbsp;
----
&nbsp;

<!-- END -->
