\name{Dataset}
\alias{soynam}
\alias{data.check}
\alias{data.line}
\alias{gen.raw}
\alias{soybase}
\alias{data.check.qa}
\alias{data.line.qa}
\alias{gen.qa}
\alias{data.check.in}
\alias{data.line.in}
\alias{gen.in}
\title{
Datasets
}
\description{
Genotypes and phenotypes from quality assured (soybase) or original (soynam) datasets. An additional dataset containing yield components (soyin) collected at Purdue University is also available. See the section "details" for the description of data objects.

The SoyNAM population (soynam.org) is a nested association mapping panel that comprises more than 5000 recombinant inbred lines (RILs), including determinate, indeterminate, and semi-determinate genotypes from maturity groups (MG) ranging from late MG II to early MG IV, derived from 40 biparental populations, where progenies were not exposed to selection. Each biparental population approximately contains 140 individuals and all families share the cultivar IA3023 as the standard parent. From the other 40 founder parents, 17 lines are elite public germplasm from different regions, 15 have diverse ancestry and 8 are plant introductions. The SoyNAM population was designed to dissect the genetic architecture of complex traits and to map yield-associated quantitative trait loci (QTL) using a diverse panel.

Parental lines were sequenced to derive the SNP allele calls. A total of 5303 SNP loci were selected with the criterion of maximizing the number of families segregating for those loci. The SNPs were used to build the SoyNAM 6K BreadChip SNP array using the Illumina Infinium HD Assay platform (Song et al. 2017). Among those SNPs, a subset of 4312 markers were selected by the SoyNAM group as quality-assured based on proportion of missing loci and correct segregation patterns. Both raw and quality assured genotypes are available in the R package SoyNAM.

}
\usage{
data(soybase)
data(soynam)
}
\details{
Datasets of the SoyNAM project, original and quality assured (QA) versions. Data was downloaded on November 16th 2015 from \url{soynam.org}. The data collected in Indiana was collected in 2013-2014 and made available on January 2018. Studies performed on the entire dataset with additional detail about the experimental settings include Diers et al. (2018) and Xavier et al. (2018).

Genotypic matrices are named "gen.raw" and "gen.qa" for the raw and QA versions, respectively. In each dataset, phenotypes are allocated into two objects, one with the lines ("data.line") and one with checks and parents ("data.checks"). Information on data objects include year, location, environment (combination of year and location), strain, family, set (set in each environment), spot (combination of set and environment), height (in centimeters), R8 (number of days to maturity), planting date (501 represents may 1), flowering (701 represents july 1), maturity (901 represents september 1), lodging (score from 1 to 5), yield (in Kg/ha), moisture, protein (percentage in the seed), oil (percentage in the seed), fiber (percentage in the seed), seed size (in grams of 100 seeds).

The dataset including yield components collected at Purdue University (West Lafayette, Indiana) was used to investigate genomic prediction (Xavier et al. 2016) and interaction among traits (Xavier et al. 2017). This dataset contains the genotypic information in the matrix "gen.in", with missing values imputed using the software MaCH (Li et al. 2010). Similar to the datasets previously described, phenotypes are allocated into two objects, lines ("data.line.in") and checks ("data.checks.in"). Information on these data objects include year, location, environment (combination of year and location), strain, family, set (set in each environment), spot (combination of set and environment), the spatial coordinates of the field plots (BLOCK, ROW and COLUMN), plant height (in centimeters), R1 (number of days to flowering), R8 (number of days to maturity), lodging (score from 1 to 5), yield (in bu/ac), leaf shape (ratio length:width), number of nodes in the main stem, number of pods in the main stem, number of pods per node, average canopy coverage, rate of canopy coverage, growing degree day to flowering (GDD_R1), growing degree day to maturity (GDD_R1), and length of reproductive period in terms of growing degree day (GDD_REP).

}
\author{
Alencar Xavier
}
\references{

Diers, B. W., Specht, J., Rainey, K. M., Cregan, P., Song, Q., Ramasubramanian, V., ... & Shannon, G. (2018). Genetic Architecture of Soybean Yield and Agronomic Traits. G3: Genes, Genomes, Genetics, g3-200332.

Li, Y., Willer, C. J., Ding, J., Scheet, P., & Abecasis, G. R. (2010). MaCH: using sequence and genotype data to estimate haplotypes and unobserved genotypes. Genetic epidemiology, 34(8), 816-834.

Song, Q., Yan, L., Quigley, C., Jordan, B.D., Fickus, E., Schroeder, S., Song, B., An, Y. Q. C., Hyten, D., Nelson, R., Rainey, K. M., Beavis, W. D., Specht, J. E., Diers, B. W., Cregan, P. (2017). Development and Genetic Characterization of the Soybean Nested Association Mapping (NAM) Population. Plant Genome. 10(2):1-14.

Xavier, A., Muir, W. M., & Rainey, K. M. (2016). Assessing predictive properties of genome-wide selection in soybeans. G3: Genes, Genomes, Genetics, 6(8), 2611-2616.

Xavier, A., Hall, B., Casteel, S., Muir, W., & Rainey, K. M. (2017). Using unsupervised learning techniques to assess interactions among complex traits in soybeans. Euphytica, 213(8), 200.

Xavier, A., Jarquin, D., Howard, R., Ramasubramanian, V., Specht, J. E., Graef, G. L., ... & Nelson, R. (2017). Genome-Wide analysis of grain yield stability and environmental interactions in a multiparental soybean population. G3: Genes, Genomes, Genetics, g3-300300.

}
\examples{
data(soynam)
data(soybase)
data(soyin)
}