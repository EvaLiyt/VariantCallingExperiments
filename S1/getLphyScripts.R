## Libs
library(glue)
library(ape)
#install.packages("phangorn", repos = "https://cloud.r-project.org", type = "source")
#install.packages("phytools", repos = "https://cloud.r-project.org")
#install.packages("geiger", repos = "https://cloud.r-project.org")
#install.packages("TreeSim", repos = "https://cloud.r-project.org")
library(geiger)
library(TreeSim)

home_dir <- path.expand("~")

set.seed(42)
## Chosen Truth
tumourTreeLow <- 5
tumourTreeHigh <- 8

normalTreeLow <- 50
normalTreeHigh <- 90

sequence_len <- 10000
ntaxa <- 40

## define functions
get_lognormal_params <- function(lower, upper) {
  # Check inputs
  if (lower <= 0 || upper <= 0 || lower >= upper) {
    stop("Bounds must be positive and lower < upper")
  }
  
  # z-scores for 2.5% and 97.5% quantiles
  z_lower <- qnorm(0.025)
  z_upper <- qnorm(0.975)
  
  # Log-transform the bounds
  log_lower <- log(lower)
  log_upper <- log(upper)
  
  # Estimate sdlog
  sdlog <- (log_upper - log_lower) / (z_upper - z_lower)
  
  # Estimate meanlog
  meanlog <- (log_lower + log_upper) / 2
  
  return(list(meanlog = meanlog, sdlog = sdlog))
}

estimate_beta_params <- function(L, U, p_L=0.025, p_U=0.975) {
  mean_est <- (L + U) / 2
  var_est <- ((U - L) / 4)^2
  conc <- (mean_est * (1 - mean_est)) / var_est - 1
  alpha_start <- mean_est * conc
  beta_start <- (1 - mean_est) * conc
  
  objective <- function(params) {
    a <- params[1]
    b <- params[2]
    if(a <= 0 || b <= 0) return(Inf)
    q_low <- qbeta(p_L, a, b)
    q_high <- qbeta(p_U, a, b)
    (q_low - L)^2 + (q_high - U)^2
  }
  
  res <- optim(c(alpha_start, beta_start), objective, method="L-BFGS-B",
               lower = c(0.005, 0.005))
  
  return(list(alpha = res$par[1], beta = res$par[2]))
}

generate_yule_tree_length <- function(root_age, n_tips, lambda, n_sim = 100) {
  set.seed(42)
  tree_lengths <- numeric(n_sim)
  
  for (i in 1:n_sim) {
    # Simulate one Yule tree per iteration
    trees <- sim.bd.taxa.age(n = n_tips, numbsim = 1, lambda = lambda, mu = 1, age = root_age)
    
    # Handle failure cases
    if (length(trees) == 0 || is.null(trees[[1]])) {
      warning(paste("Simulation", i, "failed. Skipping."))
      tree_lengths[i] <- NA
    } else {
      tree_lengths[i] <- sum(trees[[1]]$edge.length)
    }
  }
  
  # Remove failed simulations (if any)
  tree_lengths <- na.omit(tree_lengths)
  
  return(mean(tree_lengths))
}

#### Calibrations for Normal and Tumour Root
normalRoot_params <- get_lognormal_params(normalTreeLow, normalTreeHigh)
mean_normalRoot <- normalRoot_params$meanlog
sdlog_normalRoot <- normalRoot_params$sdlog

tumourRoot_params <- get_lognormal_params(tumourTreeLow, tumourTreeHigh)
mean_tumourRoot <- tumourRoot_params$meanlog
sdlog_tumourRoot <- tumourRoot_params$sdlog


## Mutation Rate Calculation
### Calculate Mutation Rate for Colon Epithelium Normal Cells
# The data is from "https://www.pnas.org/doi/full/10.1073/pnas.1607794113"
patient_ages <- c(10, 30, 65)
ressearch_mutation_rates_low <- c(1.8, 5.5, 11.0) * 10^-7

normal_mut_model <- lm(ressearch_mutation_rates_low ~ patient_ages)
model_summary <- summary(normal_mut_model)
# slope is mutation rate per site per year
normalRate_wholeGenome <- model_summary $ coefficients[2, 1]

### Calculate Mutation Rate for Colorectal Caner
### The rate is from "https://www.nature.com/articles/nature05610/tables/1"
### assume early stage to be 3 years from paper "https://ascopubs.org/doi/10.1200/JCO.2021.39.15_suppl.e18584"
## tumour rate per site per year
tumourMuts_wholeGenome <- 1.21 * 10^-6 / 3
scale_tumour <- tumourMuts_wholeGenome/normalRate_wholeGenome

### Calculate the lambda for tumour and normal
lambda_high_normal <- sum(1 / (2:ntaxa)) / normalTreeLow
lambda_low_normal <- sum(1 / (2:ntaxa)) / normalTreeHigh
normal_lambda_params <- get_lognormal_params(lambda_low_normal, lambda_high_normal)
meanlog_lambda_normal <- normal_lambda_params$meanlog
sdlog_lambda_normal <- normal_lambda_params$sdlog

lambda_high_tumour <- sum(1 / (2:ntaxa)) / tumourTreeLow
lambda_low_tumour <- sum(1 / (2:ntaxa)) / tumourTreeHigh
tumour_lambda_params <- get_lognormal_params(lambda_low_tumour, lambda_high_tumour)
meanlog_lambda_tumour <- tumour_lambda_params$meanlog
sdlog_lambda_tumour <- tumour_lambda_params$sdlog


### Scale Branch Rates
## calculate tree length
normalTreeLength <- generate_yule_tree_length((normalTreeHigh + normalTreeLow)/2, ntaxa, (lambda_high_normal+ lambda_low_normal)/2) -
  generate_yule_tree_length((tumourTreeHigh + tumourTreeLow)/2, 20, (lambda_high_tumour + lambda_low_tumour)/2)
normalTreeLength_high <- generate_yule_tree_length(normalTreeHigh, ntaxa, lambda_low_normal) - generate_yule_tree_length(tumourTreeLow, 20, lambda_high_tumour)
normalTreeLength_low <- generate_yule_tree_length(normalTreeLow, ntaxa, lambda_high_normal) - generate_yule_tree_length(tumourTreeHigh, 20, lambda_low_tumour)

tumourTreeLength <- generate_yule_tree_length((tumourTreeHigh + tumourTreeLow)/2 , 20, lambda_high_tumour) + (normalTreeHigh + normalTreeLow)/2 - (tumourTreeHigh + tumourTreeLow)/2
tumourTreeLength_high <- generate_yule_tree_length(tumourTreeHigh, 20, lambda_low_tumour) + normalTreeHigh - tumourTreeHigh
tumourTreeLength_low <- generate_yule_tree_length(tumourTreeLow, 20, lambda_low_tumour) + normalTreeLow - tumourTreeLow

## scale the mutation rates
expect_normal_mut <- 0.01 # per site
normalRate_high <- expect_normal_mut / normalTreeLength_low
normalRate_low <- expect_normal_mut / normalTreeLength_high

tumourRate_high <- normalRate_high * scale_tumour
tumourRate_low <- normalRate_low * scale_tumour


### Generate LogNormal Distribution for Branch Rates
normalRate_params <- get_lognormal_params(normalRate_low, normalRate_high)
mean_normalRate <- normalRate_params$meanlog
sdlog_normalRate <- normalRate_params$sdlog
meanRate_normal <- exp(mean_normalRate + sdlog_normalRate / 2)

tumourRate_params <- get_lognormal_params(tumourRate_low, tumourRate_high)
mean_tumourRate <- tumourRate_params$meanlog
sdlog_tumourRate <- tumourRate_params$sdlog
meanRate_tumour <- exp(mean_tumourRate + sdlog_tumourRate / 2)


# ADO rate
delta_high <- 0.3
delta_low <- 0.1
params_delta <- estimate_beta_params(delta_low, delta_high)

alpha_delta <- params_delta$alpha
beta_delta <- params_delta$beta

# t (coverage) and v (variance)
cov_high <- 6
cov_low <- 5
cov_params <- get_lognormal_params(cov_low, cov_high)
sdlog_cov <- cov_params$sdlog
mean_cov <- cov_params$meanlog

v_high <- 25
v_low <- 15
v_params <- get_lognormal_params(v_low, v_high)
sdlog_v <- v_params$sdlog
mean_v <- v_params$meanlog

# epsilon (effective sequencing rate)
epsilon_high <- 0.0025
epsilon_low <- 0.00076

params_epsilon <- estimate_beta_params(epsilon_low, epsilon_high)

alpha_epsilon <- params_epsilon$alpha
beta_epsilon <- params_epsilon$beta

# overdispersion
w1_high <- 55
w1_low <- 50
w1_params <- get_lognormal_params(w1_low, w1_high)
mean_w1 <- w1_params$meanlog
sdlog_w1 <- 0.33

w2_low <- 2.5
w2_high <- 2.55
w2_params <- get_lognormal_params(w2_low, w2_high)
mean_w2 <- w2_params$meanlog
sdlog_w2 <- 0.33

## Simplest Model (SM)
### The simplest model is the Yule tree contains 40 colorectal tumour cells only.
# mutation rate for 40 taxa tumour tree
tumourTreeLength_high_SM <- generate_yule_tree_length(tumourTreeHigh, ntaxa, lambda_low_tumour)
tumourTreeLength_low_SM <- generate_yule_tree_length(tumourTreeLow, ntaxa, lambda_low_tumour)

## scale the mutation rates
mu_high <- expect_normal_mut* scale_tumour / tumourTreeLength_low_SM
mu_low <- expect_normal_mut* scale_tumour / tumourTreeLength_high_SM

mu_params <- get_lognormal_params(mu_low, mu_high)
mean_mu <- mu_params$meanlog
sdlog_mu <- mu_params$sdlog
meanRate_mu <- exp(mean_mu + sdlog_mu / 2)

# Create output folder if it doesn't exist
if (!dir.exists("./simpleModel/1-simulation/lphy_scripts")) {
  dir.create("./simpleModel")
  dir.create("./simpleModel/1-simulation")
  dir.create("./simpleModel/1-simulation/lphy_scripts")
}
for (i in 0:99) {
  script <- glue('
data {{
    RD = readFasta(file = \"{home_dir}/data/10K.fna\");
    L = RD.nchar();
    n = {ntaxa};
}}
model{{
    snp ~ SNPSampler(alignment=RD, p=0.001);
    D = injectSNP(alignment=RD, snp=snp);
    FD = fasta(haploid(D));

    lambda ~ LogNormal(meanlog = {meanlog_lambda_tumour}, sdlog = {sdlog_lambda_tumour});
    mu ~ LogNormal(meanlog ={mean_mu}, sdlog = {sdlog_mu});

    ψ ~ Yule(lambda=lambda, n=n);

    // gt16 model
    π ~ Dirichlet(conc=[3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0]);
    rates ~ Dirichlet(conc=[1.0, 2.0, 1.0, 1.0, 2.0, 1.0]);
    Q = gt16(freq=π, rates=rates);
    A ~ PhyloCTMC(L=L, Q=Q, tree=ψ, dataType=phasedGenotype(), mu=mu, root=D);

    conc = rep(8.0, {ntaxa});
    s ~ Dirichlet(conc=conc, sum={ntaxa});
    
    delta ~ Beta(alpha = {alpha_delta}, beta = {beta_delta});
    t ~ LogNormal(meanlog = {mean_cov}, sdlog = {sdlog_cov});
    v ~ LogNormal(meanlog = {mean_v}, sdlog = {sdlog_v});
    epsilon ~ Beta(alpha = {alpha_epsilon}, beta = {beta_epsilon});
    w1 ~ LogNormal(meanlog = {mean_w1}, sdlog = {sdlog_w1});
    w2 ~ LogNormal(meanlog = {mean_w2}, sdlog = {sdlog_w2});
    
    alpha ~ Ploidy(l= L, n= n, delta= delta);
    cov ~ CoverageModel(alpha= alpha, t= t, v= v, s= s); // mimic CRC28 dataset
    r ~ ReadCountModel(D=A, epsilon=epsilon, alpha=alpha, coverage=cov, w1=w1, w2=w2);
}}
')
  
  file_name <- sprintf("simpleModel/1-simulation/lphy_scripts/sequenceGeneration_r%d.lphy", i)
  writeLines(script, file_name)
}

## inference model
# Create output folder if it doesn't exist
if (!dir.exists("./simpleModel/3-xml/TrueErrModelScripts")) {
  dir.create("./simpleModel/3-xml")
  dir.create("./simpleModel/3-xml/TrueErrModelScripts")
}
for (i in 0:99) {
  script <- glue('
data {{
    // alignment
    RA = readNexus(file = "../../1-simulation/nexus/sequenceGeneration_r{i}_A.nexus");
    var = variableSites(alignment=RA);
    A = copySites(ids=var, alignment=RA);
    L = A.nchar();
    effectiveMeanRate={meanRate_mu}*{sequence_len}/L;
    logMean = log(effectiveMeanRate);
}}
model{{
    lambda ~ LogNormal(meanlog = {meanlog_lambda_tumour}, sdlog = {sdlog_lambda_tumour});
    mu ~ LogNormal(meanlog = logMean, sdlog=1);
    tree ~ Yule(lambda=lambda, n={ntaxa});

    // gt16 model
    π ~ Dirichlet(conc=[3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0]);
    rates ~ Dirichlet(conc=[1.0, 2.0, 1.0, 1.0, 2.0, 1.0]);
    Q = gt16(freq=π, rates=rates);
    D ~ PhyloCTMC(L=L, Q=Q, tree=tree, mu=mu, dataType=phasedGenotype());
    epsilon ~ Beta(alpha=2, beta=18);
    delta ~ Beta(alpha=1.5, beta=4.5);
    A ~ GT16ErrorModel(delta=delta, epsilon=epsilon, alignment=D);
}}
')
  
  file_name <- sprintf("simpleModel/3-xml/TrueErrModelScripts/TrueErrModel_r%d.lphy", i)
  writeLines(script, file_name)
}
# Create output folder if it doesn't exist
if (!dir.exists("./simpleModel/3-xml/TruewtErrModelScripts")) {
  dir.create("./simpleModel/3-xml/TruewtErrModelScripts")
}
for (i in 0:99) {
  script <- glue('
data {{
    // alignment
    RA = readNexus(file = "../../1-simulation/nexus/sequenceGeneration_r{i}_A.nexus");
         
    var = variableSites(alignment=RA);
    A = copySites(ids=var, alignment=RA);
    L = A.nchar();
    effectiveMeanRate={meanRate_mu}*{sequence_len}/L;
    logMean = log(effectiveMeanRate);
}}
model{{
    lambda ~ LogNormal(meanlog = {meanlog_lambda_tumour}, sdlog = {sdlog_lambda_tumour});
    mu ~ LogNormal(meanlog = logMean, sdlog=1);
    tree ~ Yule(lambda=lambda, n={ntaxa});

    // gt16 model
    π ~ Dirichlet(conc=[3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0]);
    rates ~ Dirichlet(conc=[1.0, 2.0, 1.0, 1.0, 2.0, 1.0]);
    Q = gt16(freq=π, rates=rates);
    A ~ PhyloCTMC(L=L, Q=Q, tree=tree, mu=mu, dataType=phasedGenotype());
}}
')
  
  file_name <- sprintf("simpleModel/3-xml/TruewtErrModelScripts/TruewtErrModel_r%d.lphy", i)
  writeLines(script, file_name)
}


## Calibrated Yule Model (CM)
## The CalibratedYule model is a tree with two calibration nodes: the root and the mrca of tumour cells. The tree contain 20 normal cells and 20 tumour cells.
### Simulation
#### Print Lines
tumour_names_line <- paste0('"', paste0("clade_", 0:19), '"', collapse = ", ")
if (!dir.exists("./calibratedModel/1-simulation/")) {
  dir.create("./calibratedModel")
  dir.create("./calibratedModel/1-simulation/")
}
calibration_model_script <- glue('
data {{
    RD = readFasta(file = "{home_dir}/data/10K.fna");
    L = RD.nchar();
    n = {ntaxa};
}}
model {{
    snp ~ SNPSampler(alignment=RD, p=0.001);
    D = injectSNP(alignment=RD, snp=snp);

    lambda ~ LogNormal(meanlog = {meanlog_lambda_normal}, sdlog ={sdlog_lambda_normal});
    rootAge ~ LogNormal(meanlog = {mean_normalRoot}, sdlog = {sdlog_normalRoot});
    tumourAge ~ LogNormal(meanlog = {mean_tumourRoot}, sdlog = {sdlog_tumourRoot});
    rootRate ~ LogNormal(meanlog = {mean_normalRate}, sdlog = {sdlog_normalRate});
    tumourRate ~ LogNormal(meanlog = {mean_tumourRate}, sdlog = {sdlog_tumourRate});

    tumourNames = [{tumour_names_line}];

    tumourTaxa = taxa(names = tumourNames);

    tree ~ CalibratedYule(lambda=lambda, n=n, cladeTaxa=tumourTaxa, cladeMRCAAge=[tumourAge], rootAge=rootAge);

    tumourRoot = mrca(tree=tree, taxa=tumourNames);
    branchRates = localClock(tree=tree, clades=[tumourRoot], cladeRates=[tumourRate], rootRate=rootRate, includeStem=true);

    // gt16 model
    π ~ Dirichlet(conc=[3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0]);
    rates ~ Dirichlet(conc=[1.0, 2.0, 1.0, 1.0, 2.0, 1.0]);
    Q = gt16(freq=π, rates=rates);

    A ~ PhyloCTMC(L=L, Q=Q, tree=tree, dataType=phasedGenotype(), branchRates = branchRates, root=D);
}}
')

writeLines(calibration_model_script, "./calibratedModel/1-simulation/calibrationSimulation.lphy")

# Create output folder if it doesn't exist
if (!dir.exists("./calibratedModel/1-simulation/lphy_scripts")) {
  dir.create("./calibratedModel/1-simulation/lphy_scripts")
}
for (i in 0:99) {
  script <- glue('
data {{
    treeList = readTrees(file="../trueTrees/tree_r{i}.trees");
    tree=treeList[0];
    RD = readFasta(file = "{home_dir}/data/10K.fna");
    L = RD.nchar();
    n = {ntaxa};
}}
model {{
    snp ~ SNPSampler(alignment=RD, p=0.001);
    D = injectSNP(alignment=RD, snp=snp);
    FD = fasta(haploid(D));

    lambda ~ LogNormal(meanlog = {meanlog_lambda_normal}, sdlog ={sdlog_lambda_normal});
    rootAge ~ LogNormal(meanlog = {mean_normalRoot}, sdlog = {sdlog_normalRoot});
    tumourAge ~ LogNormal(meanlog = {mean_tumourRoot}, sdlog = {sdlog_tumourRoot});
    rootRate ~ LogNormal(meanlog = {mean_normalRate}, sdlog = {sdlog_normalRate});
    tumourRate ~ LogNormal(meanlog = {mean_tumourRate}, sdlog = {sdlog_tumourRate});

    tumourNames = [{tumour_names_line}];

    tumourTaxa = taxa(names = tumourNames);

    tree ~ CalibratedYule(lambda=lambda, n=n, cladeTaxa=tumourTaxa, cladeMRCAAge=[tumourAge], rootAge=rootAge);

    tumourRoot = mrca(tree=tree, taxa=tumourNames);
    branchRates = localClock(tree=tree, clades=[tumourRoot], cladeRates=[tumourRate], rootRate=rootRate, includeStem=true);

    // gt16 model
    π ~ Dirichlet(conc=[3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0]);
    rates ~ Dirichlet(conc=[1.0, 2.0, 1.0, 1.0, 2.0, 1.0]);
    Q = gt16(freq=π, rates=rates);

    A ~ PhyloCTMC(L=L, Q=Q, tree=tree, branchRates=branchRates, dataType=phasedGenotype(), root=D);
    conc = rep(8.0, {ntaxa});
    s ~ Dirichlet(conc=conc, sum={ntaxa});
    
    delta ~ Beta(alpha = {alpha_delta}, beta = {beta_delta});
    t ~ LogNormal(meanlog = {mean_cov}, sdlog = {sdlog_cov});
    v ~ LogNormal(meanlog = {mean_v}, sdlog = {sdlog_v});
    epsilon ~ Beta(alpha = {alpha_epsilon}, beta = {beta_epsilon});
    w1 ~ LogNormal(meanlog = {mean_w1}, sdlog = {sdlog_w1});
    w2 ~ LogNormal(meanlog = {mean_w2}, sdlog = {sdlog_w2});
    
    alpha ~ Ploidy(l= L, n= n, delta= delta);
    cov ~ CoverageModel(alpha= alpha, t= t, v= v, s= s); // mimic CRC28 dataset
    r ~ ReadCountModel(D=A, epsilon=epsilon, alpha=alpha, coverage=cov, w1=w1, w2=w2);
}}
')
  
  file_name <- sprintf("calibratedModel/1-simulation/lphy_scripts/sequenceGeneration_r%d.lphy", i)
  writeLines(script, file_name)
}

## inference model
# Create output folder if it doesn't exist
if (!dir.exists("./calibratedModel/3-xml/TrueErrModelScripts")) {
  dir.create("./calibratedModel/3-xml")
  dir.create("./calibratedModel/3-xml/TrueErrModelScripts")
}
for (i in 0:99) {
  script <- glue('
data {{
    // alignment
    RA = readNexus(file = "../../1-simulation/nexus/sequenceGeneration_r{i}_A.nexus");
    var = variableSites(alignment=RA);
    A = copySites(ids=var, alignment=RA);
    L = A.nchar();
    effectiveNormalRate={meanRate_normal}*{sequence_len}/L;
    logMeanNormal = log(effectiveNormalRate);
    effectivetumourRate={meanRate_tumour}*{sequence_len}/L;
    logMeanTumour = log(effectivetumourRate);
}}
model{{
    lambda ~ LogNormal(meanlog = {meanlog_lambda_normal}, sdlog = {sdlog_lambda_normal});
    rootAge ~ LogNormal(meanlog = {mean_normalRoot}, sdlog = {sdlog_normalRoot});
    tumourAge ~ LogNormal(meanlog = {mean_tumourRoot}, sdlog = {sdlog_tumourRoot});
    rootRate ~ LogNormal(meanlog = logMeanNormal, sdlog=1);
    tumourRate ~ LogNormal(meanlog = logMeanTumour, sdlog=1);

    tumourNames = [{tumour_names_line}];

    tumourTaxa = taxa(names = tumourNames);

    tree ~ CalibratedYule(lambda=lambda, n={ntaxa}, cladeTaxa=tumourTaxa, cladeMRCAAge=[tumourAge], rootAge=rootAge);

    tumourRoot = mrca(tree=tree, taxa=tumourNames);
    branchRates = localClock(tree=tree, clades=[tumourRoot], cladeRates=[tumourRate], rootRate=rootRate, includeStem=true);

    // gt16 model
    π ~ Dirichlet(conc=[3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0]);
    rates ~ Dirichlet(conc=[1.0, 2.0, 1.0, 1.0, 2.0, 1.0]);
    Q = gt16(freq=π, rates=rates);

    D ~ PhyloCTMC(L=L, Q=Q, tree=tree, branchRates=branchRates, dataType=phasedGenotype());
    
    epsilon ~ Beta(alpha=2, beta=18);
    delta ~ Beta(alpha=1.5, beta=4.5);
    A ~ GT16ErrorModel(delta=delta, epsilon=epsilon, alignment=D);
}}
')
  
  file_name <- sprintf("calibratedModel/3-xml/TrueErrModelScripts/TrueErrModel_r%d.lphy", i)
  writeLines(script, file_name)
}
# Create output folder if it doesn't exist
if (!dir.exists("./calibratedModel/3-xml/TruewtErrModelScripts")) {
  dir.create("./calibratedModel/3-xml/TruewtErrModelScripts")
}
for (i in 0:99) {
  script <- glue('
data {{
    // alignment
    RA = readNexus(file = "../../1-simulation/nexus/sequenceGeneration_r{i}_A.nexus");
         
    var = variableSites(alignment=RA);
    A = copySites(ids=var, alignment=RA);
     L = A.nchar();
    effectiveNormalRate={meanRate_normal}*{sequence_len}/L;
    logMeanNormal = log(effectiveNormalRate);
    effectivetumourRate={meanRate_tumour}*{sequence_len}/L;
    logMeanTumour = log(effectivetumourRate);
}}
model{{
    lambda ~ LogNormal(meanlog = {meanlog_lambda_normal}, sdlog = {sdlog_lambda_normal});
    rootAge ~ LogNormal(meanlog = {mean_normalRoot}, sdlog = {sdlog_normalRoot});
    tumourAge ~ LogNormal(meanlog = {mean_tumourRoot}, sdlog = {sdlog_tumourRoot});
    rootRate ~ LogNormal(meanlog = logMeanNormal, sdlog=1);
    tumourRate ~ LogNormal(meanlog = logMeanTumour, sdlog=1);

    tumourNames = [{tumour_names_line}];

    tumourTaxa = taxa(names = tumourNames);

    tree ~ CalibratedYule(lambda=lambda, n={ntaxa}, cladeTaxa=tumourTaxa, cladeMRCAAge=[tumourAge], rootAge=rootAge);

    tumourRoot = mrca(tree=tree, taxa=tumourNames);
    branchRates = localClock(tree=tree, clades=[tumourRoot], cladeRates=[tumourRate], rootRate=rootRate, includeStem=true);

    // gt16 model
    π ~ Dirichlet(conc=[3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0]);
    rates ~ Dirichlet(conc=[1.0, 2.0, 1.0, 1.0, 2.0, 1.0]);
    Q = gt16(freq=π, rates=rates);

    A ~ PhyloCTMC(L=L, Q=Q, tree=tree, branchRates=branchRates, dataType=phasedGenotype());
}}
')
  
  file_name <- sprintf("calibratedModel/3-xml/TruewtErrModelScripts/TruewtErrModel_r%d.lphy", i)
  writeLines(script, file_name)
}

## Realistic Model (RM)
### Have a 10,000 taxa tree as healthy tree with normal ages and a 10,000 taxa tree as tumour tree with tumour age.
### Sample a branch to substitute on healthy tree with the tumour tree.
### Sample 20 normal cells and 20 tumour cells to get subsampled tree.
### Simulation
#### Lambda
bigN <- 10000
lambda_high_realistic_normal <- sum(1 / (2:bigN)) / normalTreeLow
lambda_low_realistic_normal <- sum(1 / (2:bigN)) / normalTreeHigh
bigTreeLambdaParams_normal <- get_lognormal_params(lambda_low_realistic_normal, lambda_high_realistic_normal)
mean_lambda_RM_normal <- bigTreeLambdaParams_normal$meanlog
sdlog_lambda_RM_normal <- bigTreeLambdaParams_normal$sdlog

lambda_high_realistic_tumour <- sum(1 / (2:bigN)) / tumourTreeLow
lambda_low_realistic_tumour <- sum(1 / (2:bigN)) / tumourTreeHigh
bigTreeLambdaParams_tumour <- get_lognormal_params(lambda_low_realistic_tumour, lambda_high_realistic_tumour)
mean_lambda_RM_tumour <- bigTreeLambdaParams_tumour$meanlog
sdlog_lambda_RM_tumour <- bigTreeLambdaParams_tumour$sdlog
#### Print Lines
# Create output folder if it doesn't exist
if (!dir.exists("./complexModel/1-simulation/lphy_scripts")) {
  dir.create("./complexModel")
  dir.create("./complexModel/1-simulation")
  dir.create("./complexModel/1-simulation/lphy_scripts")
}
for (i in 0:99) {
  script <- glue('
data {{
    RD = readFasta(file = \"{home_dir}/data/10K.fna\");
    L = RD.nchar();
    n = {ntaxa};
    bigN = {bigN};
}}
model{{
    snp ~ SNPSampler(alignment=RD, p=0.001);
    D = injectSNP(alignment=RD, snp=snp);

    lambda ~ LogNormal(meanlog = {mean_lambda_RM_normal}, sdlog ={sdlog_lambda_RM_normal});
    tumourAge ~ LogNormal(meanlog = {mean_tumourRoot}, sdlog = {sdlog_tumourRoot});
    rootRate ~ LogNormal(meanlog = {mean_normalRate}, sdlog = {sdlog_normalRate});
    tumourRate ~ LogNormal(meanlog = {mean_tumourRate}, sdlog = {sdlog_tumourRate});
    lambdaT ~ LogNormal(meanlog = {mean_lambda_RM_tumour}, sdlog ={sdlog_lambda_RM_tumour});
    
    // healthy Yule tree
    ψhealthy ~ Yule(lambda=lambda, n=bigN);

    // sample a branch at the time
    node ~ SampleBranch(tree=ψhealthy, age=tumourAge);

    // tumour Yule tree with known root age
    ψtumour ~ Yule(lambda=lambdaT, n=bigN, rootAge=tumourAge);

    // generate the full tree
    nodeLabel = "tumourRoot";
    ψfull = substituteClade(baseTree=ψhealthy, cladeTree=ψtumour, node=node, nodeLabel=nodeLabel);
    tumourRoot = ψfull.getLabeledNode(nodeLabel);

    // get the names for healthy and tumour cells
    allTaxa = ψfull.getTaxaNames();
    tumourTaxa = tumourRoot.getLeafNames();
    healthyTaxa = setDifference(mainSet=allTaxa, excludeSet=tumourTaxa);

    sampledTumourTaxa ~ sample(array=tumourTaxa, size=20);
    sampledHealthyTaxa ~ sample(array=healthyTaxa, size=20);
    scsTaxa = setUnion(firstSet=sampledTumourTaxa, secondSet=sampledHealthyTaxa);

    otherHealthy = setDifference(mainSet=healthyTaxa, excludeSet=sampledHealthyTaxa);
    bulkHealthyTaxa ~ sample(array=otherHealthy, size=980);
    tips = setUnion(firstSet=scsTaxa, secondSet=bulkHealthyTaxa);
    
    // sample the tree
    ψ ~ SubsampledTree(tree=ψfull, taxa=[tips], sampleFraction=[1.0]);
    sampledTumourRoot = mrca(tree=ψ , taxa=sampledTumourTaxa);
    ψlabelled = labelClade(tree=ψ, taxa=sampledTumourTaxa, label=nodeLabel);
   
    // gt16 model
    conc = rep(3.0, 16);
    π ~ Dirichlet(conc=conc);
    rates ~ Dirichlet(conc=[1.0, 2.0, 1.0, 1.0, 2.0, 1.0]);
    Q = gt16(freq=π, rates=rates);
    
     // generate branch rates
    branchRates = localClock(tree=ψ, clades=[sampledTumourRoot], cladeRates=[tumourRate], rootRate=rootRate, includeStem=true);

    A ~ PhyloCTMC(L=L, Q=Q, tree=ψ, dataType=phasedGenotype(), branchRates = branchRates, root=D);
    scA = extractAlignment(alignment=A, taxa=scsTaxa);
    bulkNames=setDifference(mainSet=tips, excludeSet=sampledTumourTaxa);
    bulkA = extractAlignment(alignment=A,taxa=bulkNames);
    Fbulk=fasta(haploid(bulkA));
     
    conc = rep(8.0, {ntaxa});
    s ~ Dirichlet(conc=conc, sum={ntaxa});
    
    delta ~ Beta(alpha = {alpha_delta}, beta = {beta_delta});
    t ~ LogNormal(meanlog = {mean_cov}, sdlog = {sdlog_cov});
    v ~ LogNormal(meanlog = {mean_v}, sdlog = {sdlog_v});
    epsilon ~ Beta(alpha = {alpha_epsilon}, beta = {beta_epsilon});
    w1 ~ LogNormal(meanlog = {mean_w1}, sdlog = {sdlog_w1});
    w2 ~ LogNormal(meanlog = {mean_w2}, sdlog = {sdlog_w2});
    
    alpha ~ Ploidy(l= L, n= n, delta= delta);
    cov ~ CoverageModel(alpha= alpha, t= t, v= v, s= s); // mimic CRC28 dataset
    r ~ ReadCountModel(D=scA, epsilon=epsilon, alpha=alpha, coverage=cov, w1=w1, w2=w2);
}}
')
  
  file_name <- sprintf("complexModel/1-simulation/lphy_scripts/sequenceGeneration_r%d.lphy", i)
  writeLines(script, file_name)
}

## inference model
# Create output folder if it doesn't exist
if (!dir.exists("./complexModel/3-xml/TrueErrModelScripts")) {
  dir.create("./complexModel/3-xml")
  dir.create("./complexModel/3-xml/TrueErrModelScripts")
}
for (i in 0:99) {
  script <- glue('
data {{
    // getting information from tree tumour root and taxa
    treeList = readTrees(file = "../../1-simulation/trees/sequenceGeneration_r{i}_psilabelled.trees");
    trueTree = treeList[0];

    n = trueTree.leafCount();
    nodeLabel = "tumourRoot";
    tumourRoot_scs = trueTree.getLabeledNode(nodeLabel);
    tumourLeaves = tumourRoot_scs.getLeafNames();
    normalLeaves = setDifference(trueTree.getTaxaNames(), tumourLeaves);
    
    // alignment
    RA = readNexus(file = "../../1-simulation/nexus/sequenceGeneration_r{i}_scA.nexus");
         
    var = variableSites(alignment=RA);
    A = copySites(ids=var, alignment=RA);
    L = A.nchar();
    effectiveNormalRate={meanRate_normal}*{sequence_len}/L;
    logMeanNormal = log(effectiveNormalRate);
    effectivetumourRate={meanRate_tumour}*{sequence_len}/L;
    logMeanTumour = log(effectivetumourRate);
}}
model{{
    lambda ~ LogNormal(meanlog = {meanlog_lambda_normal}, sdlog = {sdlog_lambda_normal});
    tumourAge ~ LogNormal(meanlog = {mean_tumourRoot}, sdlog = {sdlog_tumourRoot});
    rootRate ~ LogNormal(meanlog = logMeanNormal, sdlog=1);
    tumourRate ~ LogNormal(meanlog = logMeanTumour, sdlog=1);
    
    tumourTaxa=taxa(names=tumourLeaves);
    normalTaxa=taxa(names=normalLeaves);
    tree ~ CalibratedYule(lambda=lambda, n=n, cladeTaxa=tumourTaxa, cladeMRCAAge=[tumourAge], otherTaxa=normalTaxa);
    tumourRoot = mrca(tree=tree, taxa=tumourLeaves);
    branchRates = localClock(tree=tree, clades=[tumourRoot], cladeRates=[tumourRate], rootRate=rootRate, includeStem=true);

    // gt16 model
    π ~ Dirichlet(conc=[3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0]);
    rates ~ Dirichlet(conc=[1.0, 2.0, 1.0, 1.0, 2.0, 1.0]);
    Q = gt16(freq=π, rates=rates);

    D ~ PhyloCTMC(L=L, Q=Q, tree=tree, branchRates=branchRates, dataType=phasedGenotype());
    
    epsilon ~ Beta(alpha=2, beta=18);
    delta ~ Beta(alpha=1.5, beta=4.5);
    A ~ GT16ErrorModel(delta=delta, epsilon=epsilon, alignment=D);
}}
')
  
  file_name <- sprintf("complexModel/3-xml/TrueErrModelScripts/TrueErrModel_r%d.lphy", i)
  writeLines(script, file_name)
}
# Create output folder if it doesn't exist
if (!dir.exists("./complexModel/3-xml/TruewtErrModelScripts")) {
  dir.create("./complexModel/3-xml/TruewtErrModelScripts")
}
for (i in 0:99) {
  script <- glue('
data {{
     // getting information from tree tumour root and taxa
    treeList = readTrees(file = "../../1-simulation/trees/sequenceGeneration_r{i}_psilabelled.trees");
    trueTree = treeList[0];

    n = trueTree.leafCount();
    nodeLabel = "tumourRoot";
    tumourRoot_scs = trueTree.getLabeledNode(nodeLabel);
    tumourLeaves = tumourRoot_scs.getLeafNames();
    normalLeaves = setDifference(trueTree.getTaxaNames(), tumourLeaves);
    
    // alignment
    RA = readNexus(file = "../../1-simulation/nexus/sequenceGeneration_r{i}_scA.nexus");
    
    var = variableSites(alignment=RA);
    A = copySites(ids=var, alignment=RA);
    L = A.nchar();
    effectiveNormalRate={meanRate_normal}*{sequence_len}/L;
    logMeanNormal = log(effectiveNormalRate);
    effectivetumourRate={meanRate_tumour}*{sequence_len}/L;
    logMeanTumour = log(effectivetumourRate);
}}
model{{
    lambda ~ LogNormal(meanlog = {meanlog_lambda_normal}, sdlog = {sdlog_lambda_normal});
    tumourAge ~ LogNormal(meanlog = {mean_tumourRoot}, sdlog = {sdlog_tumourRoot});
    rootRate ~ LogNormal(meanlog = logMeanNormal, sdlog=1);
    tumourRate ~ LogNormal(meanlog = logMeanTumour, sdlog=1);
    
    tumourTaxa=taxa(names=tumourLeaves);
    normalTaxa=taxa(names=normalLeaves);
    tree ~ CalibratedYule(lambda=lambda, n=n, cladeTaxa=tumourTaxa, cladeMRCAAge=[tumourAge], otherTaxa=normalTaxa);
    tumourRoot = mrca(tree=tree, taxa=tumourLeaves);
    branchRates = localClock(tree=tree, clades=[tumourRoot], cladeRates=[tumourRate], rootRate=rootRate, includeStem=true);

    // gt16 model
    π ~ Dirichlet(conc=[3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0]);
    rates ~ Dirichlet(conc=[1.0, 2.0, 1.0, 1.0, 2.0, 1.0]);
    Q = gt16(freq=π, rates=rates);

    A ~ PhyloCTMC(L=L, Q=Q, tree=tree, branchRates=branchRates, dataType=phasedGenotype());
}}
')
  
  file_name <- sprintf("complexModel/3-xml/TruewtErrModelScripts/TruewtErrModel_r%d.lphy", i)
  writeLines(script, file_name)
}

## DxVx inference models
dxvx_pairs <- list(
  c(1,1), c(2,1), c(2,2), c(3,1), c(3,2), c(3,3),
  c(4,1), c(4,2), c(4,3), c(4,4),
  c(5,1), c(5,2), c(5,3), c(5,4), c(5,5)
)

## Simple Model DxVx
for (pair in dxvx_pairs) {
  d <- pair[1]; v <- pair[2]
  DV_folder <- sprintf("D%dV%d", d, v)
  align_suffix <- sprintf("_d%d_v%d", d, v)

  # ErrModel
  err_dir <- sprintf("./simpleModel/3-xml/%s/ErrModelScripts", DV_folder)
  if (!dir.exists(err_dir)) dir.create(err_dir, recursive = TRUE)
  for (i in 0:99) {
    script <- glue('
data {{
   // alignment
    RA = readFasta(file = "../../../2-variantCalling/{DV_folder}/alignment/alignment_r{i}{align_suffix}.fasta", sequenceType=phasedGenotype());
    var = variableSites(alignment=RA);
    A = copySites(ids=var, alignment=RA);
    L = A.nchar();
    effectiveMeanRate={meanRate_mu}*{sequence_len}/L;
    logMean = log(effectiveMeanRate);
}}
model{{
    lambda ~ LogNormal(meanlog = {meanlog_lambda_tumour}, sdlog = {sdlog_lambda_tumour});
    mu ~ LogNormal(meanlog = logMean, sdlog=1);
    tree ~ Yule(lambda=lambda, n={ntaxa});

    // gt16 model
    π ~ Dirichlet(conc=[3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0]);
    rates ~ Dirichlet(conc=[1.0, 2.0, 1.0, 1.0, 2.0, 1.0]);
    Q = gt16(freq=π, rates=rates);
    D ~ PhyloCTMC(L=L, Q=Q, tree=tree, dataType=phasedGenotype(), mu=mu);

    epsilon ~ Beta(alpha=2, beta=18);
    delta ~ Beta(alpha=1.5, beta=4.5);
    A ~ GT16ErrorModel(delta=delta, epsilon=epsilon, alignment=D);
}}
')
    file_name <- sprintf("%s/ErrModel_r%d%s.lphy", err_dir, i, align_suffix)
    writeLines(script, file_name)
  }

  # wtErrModel
  wt_dir <- sprintf("./simpleModel/3-xml/%s/wtErrModelScripts", DV_folder)
  if (!dir.exists(wt_dir)) dir.create(wt_dir, recursive = TRUE)
  for (i in 0:99) {
    script <- glue('
data {{
   // alignment
    RA = readFasta(file = "../../../2-variantCalling/{DV_folder}/alignment/alignment_r{i}{align_suffix}.fasta", sequenceType=phasedGenotype());
    var = variableSites(alignment=RA);
    A = copySites(ids=var, alignment=RA);
    L = A.nchar();
    effectiveMeanRate={meanRate_mu}*{sequence_len}/L;
    logMean = log(effectiveMeanRate);
}}
model{{
    lambda ~ LogNormal(meanlog = {meanlog_lambda_tumour}, sdlog = {sdlog_lambda_tumour});
    mu ~ LogNormal(meanlog = logMean, sdlog=1);
    tree ~ Yule(lambda=lambda, n={ntaxa});

    // gt16 model
    π ~ Dirichlet(conc=[3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0]);
    rates ~ Dirichlet(conc=[1.0, 2.0, 1.0, 1.0, 2.0, 1.0]);
    Q = gt16(freq=π, rates=rates);
    A ~ PhyloCTMC(L=L, Q=Q, tree=tree, dataType=phasedGenotype(), mu=mu);
}}
')
    file_name <- sprintf("%s/wtErrModel_r%d%s.lphy", wt_dir, i, align_suffix)
    writeLines(script, file_name)
  }
}

## Calibrated Model DxVx
for (pair in dxvx_pairs) {
  d <- pair[1]; v <- pair[2]
  DV_folder <- sprintf("D%dV%d", d, v)
  align_suffix <- sprintf("_d%d_v%d", d, v)

  # ErrModel
  err_dir <- sprintf("./calibratedModel/3-xml/%s/ErrModelScripts", DV_folder)
  if (!dir.exists(err_dir)) dir.create(err_dir, recursive = TRUE)
  for (i in 0:99) {
    script <- glue('
data {{
   // alignment
    RA = readFasta(file = "../../../2-variantCalling/{DV_folder}/alignment/alignment_r{i}{align_suffix}.fasta", sequenceType=phasedGenotype());
    var = variableSites(alignment=RA);
    A = copySites(ids=var, alignment=RA);
    L = A.nchar();
    effectiveNormalRate={meanRate_normal}*{sequence_len}/L;
    logMeanNormal = log(effectiveNormalRate);
    effectivetumourRate={meanRate_tumour}*{sequence_len}/L;
    logMeanTumour = log(effectivetumourRate);
}}
model{{
    lambda ~ LogNormal(meanlog = {meanlog_lambda_normal}, sdlog = {sdlog_lambda_normal});
    rootAge ~ LogNormal(meanlog = {mean_normalRoot}, sdlog = {sdlog_normalRoot});
    tumourAge ~ LogNormal(meanlog = {mean_tumourRoot}, sdlog = {sdlog_tumourRoot});
    rootRate ~ LogNormal(meanlog = logMeanNormal, sdlog=1);
    tumourRate ~ LogNormal(meanlog = logMeanTumour, sdlog=1);

    tumourNames = [{tumour_names_line}];

    tumourTaxa=taxa(names=tumourNames);
    tree ~ CalibratedYule(lambda=lambda, n={ntaxa}, cladeTaxa=tumourTaxa, cladeMRCAAge=[tumourAge], rootAge=rootAge);
    tumourRoot = mrca(tree=tree, taxa=tumourNames);
    branchRates = localClock(tree=tree, clades=[tumourRoot], cladeRates=[tumourRate], rootRate=rootRate, includeStem=true);

    // gt16 model
    π ~ Dirichlet(conc=[3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0]);
    rates ~ Dirichlet(conc=[1.0, 2.0, 1.0, 1.0, 2.0, 1.0]);
    Q = gt16(freq=π, rates=rates);

    D ~ PhyloCTMC(L=L, Q=Q, tree=tree, branchRates=branchRates, dataType=phasedGenotype());

    epsilon ~ Beta(alpha=2, beta=18);
    delta ~ Beta(alpha=1.5, beta=4.5);
    A ~ GT16ErrorModel(delta=delta, epsilon=epsilon, alignment=D);
}}
')
    file_name <- sprintf("%s/ErrModel_r%d%s.lphy", err_dir, i, align_suffix)
    writeLines(script, file_name)
  }

  # wtErrModel
  wt_dir <- sprintf("./calibratedModel/3-xml/%s/wtErrModelScripts", DV_folder)
  if (!dir.exists(wt_dir)) dir.create(wt_dir, recursive = TRUE)
  for (i in 0:99) {
    script <- glue('
data {{
    // alignment
    RA = readFasta(file = "../../../2-variantCalling/{DV_folder}/alignment/alignment_r{i}{align_suffix}.fasta", sequenceType=phasedGenotype());
    var = variableSites(alignment=RA);
    A = copySites(ids=var, alignment=RA);
    L = A.nchar();
    effectiveNormalRate={meanRate_normal}*{sequence_len}/L;
    logMeanNormal = log(effectiveNormalRate);
    effectivetumourRate={meanRate_tumour}*{sequence_len}/L;
    logMeanTumour = log(effectivetumourRate);
}}
model{{
    lambda ~ LogNormal(meanlog = {meanlog_lambda_normal}, sdlog = {sdlog_lambda_normal});
    tumourAge ~ LogNormal(meanlog = {mean_tumourRoot}, sdlog = {sdlog_tumourRoot});
    rootAge ~ LogNormal(meanlog = {mean_normalRoot}, sdlog = {sdlog_normalRoot});
    rootRate ~ LogNormal(meanlog = logMeanNormal, sdlog=1);
    tumourRate ~ LogNormal(meanlog = logMeanTumour, sdlog=1);

    tumourNames = [{tumour_names_line}];

    tumourTaxa = taxa(names = tumourNames);

    tree ~ CalibratedYule(lambda=lambda, n={ntaxa}, cladeTaxa=tumourTaxa, cladeMRCAAge=[tumourAge], rootAge=rootAge);

    tumourRoot = mrca(tree=tree, taxa=tumourNames);
    branchRates = localClock(tree=tree, clades=[tumourRoot], cladeRates=[tumourRate], rootRate=rootRate, includeStem=true);

    // gt16 model
    π ~ Dirichlet(conc=[3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0]);
    rates ~ Dirichlet(conc=[1.0, 2.0, 1.0, 1.0, 2.0, 1.0]);
    Q = gt16(freq=π, rates=rates);

    A ~ PhyloCTMC(L=L, Q=Q, tree=tree, branchRates=branchRates, dataType=phasedGenotype());
}}
')
    file_name <- sprintf("%s/wtErrModel_r%d%s.lphy", wt_dir, i, align_suffix)
    writeLines(script, file_name)
  }
}

## Realistic Model DxVx
for (pair in dxvx_pairs) {
  d <- pair[1]; v <- pair[2]
  DV_folder <- sprintf("D%dV%d", d, v)
  align_suffix <- sprintf("_d%d_v%d", d, v)

  # ErrModel
  err_dir <- sprintf("./complexModel/3-xml/%s/ErrModelScripts", DV_folder)
  if (!dir.exists(err_dir)) dir.create(err_dir, recursive = TRUE)
  for (i in 0:99) {
    script <- glue('
data {{
    // getting information from tree tumour root and taxa
    treeList = readTrees(file = "../../../1-simulation/trees/sequenceGeneration_r{i}_psilabelled.trees");
    trueTree = treeList[0];

    n = trueTree.leafCount();
    nodeLabel = "tumourRoot";
    tumourRoot_scs = trueTree.getLabeledNode(nodeLabel);
    tumourLeaves = tumourRoot_scs.getLeafNames();

    // alignment
    RA = readFasta(file = "../../../2-variantCalling/{DV_folder}/alignment/alignment_r{i}{align_suffix}.fasta", sequenceType=phasedGenotype());
    taxa = RA.taxa();
    n = length(taxa);
    normalLeaves = setDifference(taxa.taxaNames(), tumourLeaves);
    
    var = variableSites(alignment=RA);
    A = copySites(ids=var, alignment=RA);
    L = A.nchar();
    effectiveNormalRate={meanRate_normal}*{sequence_len}/L;
    logMeanNormal = log(effectiveNormalRate);
    effectivetumourRate={meanRate_tumour}*{sequence_len}/L;
    logMeanTumour = log(effectivetumourRate);
}}
model{{
    lambda ~ LogNormal(meanlog = {meanlog_lambda_normal}, sdlog = {sdlog_lambda_normal});
    tumourAge ~ LogNormal(meanlog = {mean_tumourRoot}, sdlog = {sdlog_tumourRoot});
    rootRate ~ LogNormal(meanlog = logMeanNormal, sdlog=1);
    tumourRate ~ LogNormal(meanlog = logMeanTumour, sdlog=1);

    tumourTaxa=taxa(names=tumourLeaves);
    normalTaxa=taxa(names=normalLeaves);
    tree ~ CalibratedYule(lambda=lambda, n=n, cladeTaxa=tumourTaxa, cladeMRCAAge=[tumourAge], otherTaxa=normalTaxa);
    tumourRoot = mrca(tree=tree, taxa=tumourLeaves);
    branchRates = localClock(tree=tree, clades=[tumourRoot], cladeRates=[tumourRate], rootRate=rootRate, includeStem=true);

    // gt16 model
    π ~ Dirichlet(conc=[3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0]);
    rates ~ Dirichlet(conc=[1.0, 2.0, 1.0, 1.0, 2.0, 1.0]);
    Q = gt16(freq=π, rates=rates);

    D ~ PhyloCTMC(L=L, Q=Q, tree=tree, branchRates=branchRates, dataType=phasedGenotype());

    epsilon ~ Beta(alpha=2, beta=18);
    delta ~ Beta(alpha=1.5, beta=4.5);
    A ~ GT16ErrorModel(delta=delta, epsilon=epsilon, alignment=D);
}}
')
    file_name <- sprintf("%s/ErrModel_r%d%s.lphy", err_dir, i, align_suffix)
    writeLines(script, file_name)
  }

  # wtErrModel
  wt_dir <- sprintf("./complexModel/3-xml/%s/wtErrModelScripts", DV_folder)
  if (!dir.exists(wt_dir)) dir.create(wt_dir, recursive = TRUE)
  for (i in 0:99) {
    script <- glue('
data {{
    // getting information from tree tumour root and taxa
    treeList = readTrees(file = "../../../1-simulation/trees/sequenceGeneration_r{i}_psilabelled.trees");
    trueTree = treeList[0];

    n = trueTree.leafCount();
    nodeLabel = "tumourRoot";
    tumourRoot_scs = trueTree.getLabeledNode(nodeLabel);
    tumourLeaves = tumourRoot_scs.getLeafNames();

    // alignment
    RA = readFasta(file = "../../../2-variantCalling/{DV_folder}/alignment/alignment_r{i}{align_suffix}.fasta", sequenceType=phasedGenotype());
    taxa = RA.taxa();
    n = length(taxa);
    normalLeaves = setDifference(taxa.taxaNames(), tumourLeaves);
    
    var = variableSites(alignment=RA);
    A = copySites(ids=var, alignment=RA);
    L = A.nchar();
    effectiveNormalRate={meanRate_normal}*{sequence_len}/L;
    logMeanNormal = log(effectiveNormalRate);
    effectivetumourRate={meanRate_tumour}*{sequence_len}/L;
    logMeanTumour = log(effectivetumourRate);
}}
model{{
    lambda ~ LogNormal(meanlog = {meanlog_lambda_normal}, sdlog = {sdlog_lambda_normal});
    tumourAge ~ LogNormal(meanlog = {mean_tumourRoot}, sdlog = {sdlog_tumourRoot});
    rootRate ~ LogNormal(meanlog = logMeanNormal, sdlog=1);
    tumourRate ~ LogNormal(meanlog = logMeanTumour, sdlog=1);

    tumourTaxa=taxa(names=tumourLeaves);
    normalTaxa=taxa(names=normalLeaves);
    tree ~ CalibratedYule(lambda=lambda, n=n, cladeTaxa=tumourTaxa, cladeMRCAAge=[tumourAge], otherTaxa=normalTaxa);
    tumourRoot = mrca(tree=tree, taxa=tumourLeaves);
    branchRates = localClock(tree=tree, clades=[tumourRoot], cladeRates=[tumourRate], rootRate=rootRate, includeStem=true);

    // gt16 model
    π ~ Dirichlet(conc=[3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0]);
    rates ~ Dirichlet(conc=[1.0, 2.0, 1.0, 1.0, 2.0, 1.0]);
    Q = gt16(freq=π, rates=rates);

    A ~ PhyloCTMC(L=L, Q=Q, tree=tree, branchRates=branchRates, dataType=phasedGenotype());
}}
')
    file_name <- sprintf("%s/wtErrModel_r%d%s.lphy", wt_dir, i, align_suffix)
    writeLines(script, file_name)
  }
}

## Note: the realistic model's true single-cell alignment is sequenceGeneration_r{i}_scA.nexus
## (the 40 single cells), not _A.nexus (the full subsampled tree). The TrueErrModel/TruewtErrModel
## scripts above read _scA.nexus accordingly.
