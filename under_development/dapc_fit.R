

dat <- data_4pops
dat, sampCol='SAMPLE', locusCol='LOCUS', genoCol='GT', popCol='POP', num.cores=1

pcPreds <- 3

colnames(dat)[match(c(sampCol, locusCol, popCol, genoCol),colnames(dat))] <- c(
  'SAMPLE','LOCUS','POP','GT'
)

FUN_snp_da_contrib <- function(x){
  temp <- sum(x*x)
  if(temp < 1e-12) return(rep(0, length(x)))
  return(x*x / temp)
}


PCA <- pca_genos(dat, scaling=scaling, popCol=popCol)

pops <- PCA$pops %>%  as.factor()
k <- length(unique(pops))

X <- PCA$x[, 1:pcPreds] %>% as.data.frame()

DA <- lda(X, pops, prior=rep(1,k)/k, tol=1e-30)

snp.da.load <- as.matrix(PCA$rotation[, 1:pcPreds]) %*% as.matrix(DA$scaling)
snp.da.contr <- apply(snp.da.load, 2, FUN_snp_da_contrib) %>%
  as.data.frame() %>%
  rownames_to_column(., 'LOCUS') %>%
  as.data.table

DA.tab <- da2DT(DA, sampVec=rownames(X), obsPops=pops)
PCA.tab <- pca2DT(PCA) %>%
  .[AXIS %in% paste0('PC', 1:pcPreds)]

list(da.fit=DA, da.tab=DA.tab, pca.fit=PCA, pca.tab=PCA.tab, snp.contrib=snp.da.contr)

# Leave-one-out CV
samps <- dat$SAMPLE %>% unique
pops.uniq <- dat$POP %>% unique

i=1

predTab <- lapply(1:length(samps), function(i){
  PCA.i <- pca_genos(dat[SAMPLE!=samps[i],], scaling=scaling, popCol='POP')

  DA.i <- lda(
    as.data.frame(PCA.i$x[, 1:pcPreds]),
    PCA.i$pops,
    prior=rep(1,k)/k,
    tol=1e-30)

  dat.i <- dat[SAMPLE==samps[i],]

  Xnew.i <- DT2Mat_genos(dat.i) %*% PCA.i$rotation %>%
    .[, 1:pcPreds] %>%
    matrix(., ncol=pcPreds, nrow=1) %>%
    as.data.frame() %>%
    setnames(., new=paste0('PC', 1:pcPreds))

  data.table(
    POP=dat.i$POP[1],
    SAMPLE=samps[i],
    POP.PRED=predict(DA.i, newdata = Xnew.i[1,])$class
  )
}) %>%
  do.call('rbind',.)


# Training-testing partitioning
train.samps <- dat %>%
  .[, sample(SAMPLE, round(length(unique(SAMPLE)))), by=POP] %>%
  .[['V1']]

test.samps <- dat[!SAMPLE%in%train.samps]

DAPC.train <- dat[SAMPLE %in% train.samps] %>%
  adegenet_DT2genlight(
    ., sampCol=sampCol, locusCol=locusCol, genoCol=genoCol, popCol=popCol
  ) %>%
  dapc(., n.pca=pcPreds, n.da=length(genlObj$pop)-1)

predict.dapc(DAPC.train, newddata=DAPC.train)

DAPC.train$prior


# CV statistics
predTab[, sum(POP==POP.PRED)/length(SAMPLE)]

popComps <- CJ(POP=pops.uniq, POP.PRED=pops.uniq)

predPairsLong <- lapply(1:nrow(popComps), function(i){
  pop.obs <- popComps$POP[i]
  pop.pred <- popComps$POP.PRED[i]
  assign <- nrow(predTab[POP==pop.obs & POP.PRED==pop.pred])/nrow(predTab[POP==pop.obs])
  data.table(POP=pop.obs, POP.PRED=pop.pred, ASSIGN=assign)
}) %>%
  do.call('rbind', .)

predPairsWide <- predPairsLong %>%
  dcast(., POP~POP.PRED, value.var='ASSIGN')

