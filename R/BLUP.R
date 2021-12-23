BLUP=function(trait="yield",family="all",env="all",dereg=FALSE,
              MAF=0.05,use.check=TRUE,impute="FM",rm.rep=TRUE){
    
    # Line added for debugging purpose
    gen.raw=data.check=data.line=matrix(NA,2,2)
  
    # Load data
    data(soynam,envir=environment(),package="SoyNAM")
        
    # Genotypic matrix of lines
    geno = gen.raw[grep('DS1',rownames(gen.raw)),]
        
    # FAM
    fam=rownames(geno)
    fam=gsub('DS1.-','',fam)
    fam=gsub('...$','',fam,perl = T)
    fam=as.numeric(fam)
    
    # CHR
    chr=rep(NA,20)
    for(i in 1:20) chr[i]=length(grep(
      paste("Gm",sprintf("%02d",i),sep=''),
      colnames(geno)));rm(i)
    
    # Subsetting
    if(is.numeric(family)) data.line = data.line[data.line$family%in%family,]
    
    if(is.numeric(env))  if(length(env)==1) stop("At least two environments where the trait was measured are required")
    
    if(is.numeric(env)){
      E1 = as.numeric(data.line$environ)
      data.line = data.line[E1%in%env,]
      E2 = as.numeric(data.check$environ)
      data.check = data.check[E2%in%env,]
    }
    
    # Check function
    CHECK=function(trait){ test=dcast(data.check,environ+spot~strain,value.var=trait,mean)
                           rownames(test)=test[,2];E=test[,1];test=test[,-c(1,2)];test=data.matrix(test);test[is.nan(test)]=NA;
                           X=function(X) unlist(tapply(X,E,FUN=function(x){m=mean(x,na.rm=T);SD=sd(x,na.rm=T);return((x-m)/SD)}))
                           MEAN=apply(test,2,X);C=rowMeans(MEAN,na.rm=T);names(C)=rownames(test);C[is.nan(C)]=0;return(C)}
    
    # Model terms
    Y = data.line[,trait]
    G = data.line[,"strain"]
    E = data.line[,"environ"]
    
    # BLUP
    if(use.check){
      cat('solving BLUE of checks\n')
      check = CHECK(trait);set = as.character(data.line[,"spot"])
      C = check[set]
      cat('solving BLUP of phenotypes\n')
      blup=lmer(Y~C+(1|E)+(1|G))
    }else{
      cat('solving BLUP of phenotypes\n')
      blup=lmer(Y~(1|E)+(1|G))}
    
    BV = rowMeans(ranef(blup)$G)
    ge = intersect(names(BV),rownames(geno))
    BV = BV[ge]
    BV = BV+mean(Y,na.rm=T)
    
    # Check variance components
    n = c(table(blup@frame$G))[ge]
    vc = data.frame(VarCorr(blup))
    vg = vc$vcov[1]
    ve = vc$vcov[3]
    r2 = vg/(vg+ve/n)
    cat('Broad-sense H2 of',trait,'is',round(mean(r2),2),'\n')
    # Deregression
    if(dereg){BV = BV/r2}
    
    if(is.numeric(family)){
      BV = BV[fam%in%family]
      geno = geno[fam%in%family,]
      fam = fam[fam%in%family]
    }
    
    geno = snpQC(geno,MAF=MAF)
    
    cat('Removing markers with more than 50% missing values\n')
    mi=apply(geno,2,function(q)mean(is.na(q)))
    geno = geno[,mi<.5]
    
    for(i in 1:20) chr[i]=length(grep(
      paste("Gm",sprintf("%02d",i),sep=''),
      colnames(geno)));rm(i)
    
    if(anyNA(geno)){
      
      if(impute=="EXP"){
        cat('Imputing NA using allele expectation\n')
        geno = IMP(geno)
        if(anyNA(geno)){
          geno[is.na(geno)] = 1
        }
      }
      
      if(impute=="FM"){
        cat('Imputing NA using transition probabilities with neighbor SNPs \n')
        geno = markov(geno,chr)
      }
      
    }
        
    # CHR
    for(i in 1:20) chr[i]=length(grep(
      paste("Gm",sprintf("%02d",i),sep=''),
      colnames(geno)));rm(i)
    
    if(!anyNA(geno)){
      if(rm.rep){
        cat('removing duplicates\n')
        w = duplicated(c(tcrossprod(runif(ncol(geno)),geno)))
        if(any(w)){
          cat(sum(w),'individuals with identical genotypic profile\n')
          BV=BV[!w]
          geno=geno[!w,]
          fam=fam[!w]
        }
      }
    }
  
    LIST = list('Phen'=BV,'Gen'=geno,'Chrom'=chr,'Fam'=fam,'r2'=r2,'nReps'=n)
    
    return(LIST)
  }

ENV=function(trait="yield"){
  data.line.qa = 0
  gen.qa=0
  data(soynam,envir=environment(),package="SoyNAM")
  data(soybase,envir=environment(),package="SoyNAM")
  test=dcast(data.line.qa,strain~environ,value.var=trait,mean)
  Strain=as.character(test[,1])
  test=test[,-1]
  test=data.matrix(test)
  rownames(test)=Strain
  test=test[intersect(rownames(gen.qa),rownames(test)),]
  test[is.nan(test)]=NA
  test = test[,-which(apply(test,2,function(x)mean(is.na(x)))==1)]
  gen = gen.qa[rownames(test),]
  fam=as.numeric(substr(rownames(gen),6,7))
  chr=rep(NA,20)
  for(i in 1:20) chr[i]=length(grep(
    paste("Gm",sprintf("%02d",i),sep=''),
    colnames(gen)));rm(i)
  gen = markov(gen,chr)
  h = list(Y=test,gen=gen,fam=fam,chr=chr)
  return(h)}

# Extended Flexible Gwas
gwnam = function(pheno,geno,pop){
  tag = 0; numSnps = ncol(geno)
  pb = txtProgressBar(style = 3)
  sma = apply(geno,2,function(x,y,pop){
    # Print marker under evaluation
    tag <<- tag+1
    TmpDta = data.frame(y=y,f=factor(pop),x=x)
    lvl = levels(TmpDta$f)
    TmpDta = droplevels.data.frame(TmpDta[rowMeans(!is.na(TmpDta))==1,])
    Vy = c(var(TmpDta$y))
    # Null model
    fit0 = lmer(y~(1|f),TmpDta)
    ll0 = logLik(fit0)
    ## Model 1 - Within-family effect
    fit1 = suppressMessages(lmer(y~(1|f)+(1|f):x,TmpDta))
    eff1 = ranef(fit1)$f[,2]
    names(eff1) = rownames(ranef(fit1)$f)
    eff1 = eff1[lvl]
    ll1 = logLik(fit1)
    LRT1 = ll1-ll0
    PVAL1 = -log10(1-pchisq(LRT1,1))
    ## Model 2 - Across-family effect
    eff2 = suppressMessages(lmer(y~x+(1|f),TmpDta))@beta[2]
    ## Coeff of determination
    R2 = 1-c(fit0@devcomp$cmp['sigmaREML'],
             fit1@devcomp$cmp['sigmaREML'])/Vy
    names(R2)=paste0('R2.model',0:1)
    ## Output
    NumObs = nrow(TmpDta)
    out = c(NumObs,P1=PVAL1,R2,Fxd=eff2,eff1)
    setTxtProgressBar(pb, tag/numSnps)
    return(out)},y=pheno,pop=pop)
  rownames(sma) = c('NumObs','MinusLogPvalue','NullModel_R2','AltrModel_R2','OverallSnpEffect',
                    paste('SnpEffPop',sort(unique(pop)),sep=''))
  close(pb); sma[is.na(sma)] = 0
  return(data.frame(t(sma)))}