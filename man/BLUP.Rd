\name{BLUP}
\alias{BLUP}
\alias{ENV}
\alias{gwnam}
\alias{mrr}
\title{
Best Linear Unbias Predictor
}
\description{
Genetic values for a given trait computed by REML.
}
\usage{
BLUP(trait="yield",family="all",env="all",dereg=FALSE,
     MAF=0.05,use.check=TRUE,impute="FM",rm.rep=TRUE)
}
\arguments{
  \item{trait}{
Character. Trait of interest. The options are: "yield" (grain yield in Kg/ha), "maturity" (days to maturity), "height" (plant height in cm), "lodging" (lodging score from 1 to 5), "protein" (protein percentage in the grain), "oil" (oil percentage in the grain), "size" (seed size = mass of 100 seeds in grams) and "fiber" (fiber percentage in the grain).
}
  \item{family}{
Numberic vector or "all". Which SoyNAM families to use.
}
  \item{env}{
Numberic vector or "all". Which environments to use. The environments are coded as follows: 1 (IA_2012), 2 (IA_2013), 3 (IL_2011), 4 (IL_2012), 5 (IL_2013), 6 (IN_2012), 7 (IN_2013), 8 (KS_2012), 9 (KS_2013), 10 (MI_2012), 11 (MO_2012), 12 (MO_2013), 13 (NE_2011), 14 (NE_2012), 15 (OHmc_2012), 16 (OHmc_2013), 17 (OHmi_2012) and 18 (OHmi_2013).
}
  \item{dereg}{
Logic. If true, deregress BLUPs (remove shrinkage).
}
  \item{MAF}{
Numeric. Minor allele frequency threshold for the markers.
}
  \item{use.check}{
Logical. If TRUE, it includes a control term as fixed effect in the model.
}

  \item{impute}{
NULL, 'EXP' of 'FM'. If 'EXP', it imputes missing genotypes using expectation (allele frequency). If 'FM' is imputes missing genotypes using a forward Markov Chain algorithm, filling missing loci with the most likely genotype based on the previous marker.
}

  \item{rm.rep}{
Logical. If TRUE, it removes replicated genotypes. Genotypes are treated as identical when the genotypes are more than 95 percent identical. This argument requires imputed genotypes.
}
}
\details{
This function uses the raw dataset (\eqn{data(soynam)}), allowing user-defined data quality control for genotypes and BLUPs of genetic values.

The algorithm start from selecting the chosen families and environment that will be used for the best linear unbias predictor (BLUP). The BLUP values are calculates based on the following model: 

\eqn{Trait = Control + Environment + Genotype})

Where control is a covariate set as fixed effect based on the  checks of each set (microenvironment); Environment is a random effect that represents the combination of location and year; and Genotype is the random effect associated to the lines. The BLUP values are the regression coefficients corresponding to the Genotype effect. The BLUP is calculated using the R package lme4 (Bates 2010) using REML.

If checks are used as covariate (use.check=TRUE), then the best linear unbias estimator (BLUE) of the check effects is assigned to each set as a micro-environmental control. Each set had between one and five controls, including the SoyNAM parents and five other cultivars. These genotypes are normalized by environment and the BLUE of each set is calculated. All genotypes in a same set will have the same check effect.
}
\value{
This function returns a list with four objects. A numeric vector with the  BLUP solution of the phenotyes ("Phen"); the corresponding genotypes ("Gen"); a vector with the respective family ("Fam"); and a numeric vector with the number of SNPs per cromosome ("Chrom"). The output of this fuction has the exact input format for the NAM package (Xavier et al. 2015) to perform genome-wide association analysis.
}
\references{
Bates, D. M. (2010). lme4: Mixed-effects modeling with R. URL http://lme4.r-forge.r-project.org/book.

Xavier, A., Xu, S., Muir, W. M., & Rainey, K. M. (2015). NAM: association studies in multiple populations. Bioinformatics, 31(23), 3862-3864.
}
\author{
Alencar Xavier
}
\examples{
Test=BLUP(trait="yield",family=2:3,env=1:2)
}