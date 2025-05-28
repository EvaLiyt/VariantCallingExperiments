## Libs
library(glue)

## Chosen Truth
tumourTreeLow <- 5
tumourTreeHigh <- 8

normalTreeLow <- 50
normalTreeHigh <- 90

sequence_len <- 1*10^6

z975 <- qnorm(0.975)

## Mutation Rate Calculation
### Calculate Mutation Rate for Colon Epithelium Normal Cells
# The data is from "https://www.pnas.org/doi/full/10.1073/pnas.1607794113"

patient_ages <- c(10, 30, 50)
ressearch_mutation_rates_low <- c(1.8-0.5, 5.5-1.6, 11.0-1.5) * 10^-7
ressearch_mutation_rates_high <- c(1.8+0.5, 5.5+1.6, 11.0+1.5) * 10^-7

normal_mut_model_low <- lm(ressearch_mutation_rates_low ~ patient_ages)
normal_mut_model_high <- lm(ressearch_mutation_rates_high ~ patient_ages)

model_summary_low <- summary(normal_mut_model_low)
model_summary_high <- summary(normal_mut_model_high)

normalRate_low_wholeGenome <- model_summary_low $ coefficients[2, 1]
normalRate_high_wholeGenome <- model_summary_high $ coefficients[2, 1]

### Calculate Mutation Rate for Colorectal Caner
### The rate is from "https://www.nature.com/articles/nature05610/tables/1"
tumourMuts_wholeGenome <- 1.21 * 10^-6
tumourRate_wholeGenome_low <- tumourMuts_wholeGenome / tumourTreeHigh
tumourRate_wholeGenome_high <- tumourMuts_wholeGenome / tumourTreeLow

### Scale Branch Rates
### Reference genome is from "https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/"
wholeGenomeSize <- 3099441038
sizeScale <- wholeGenomeSize / sequence_len
normalRate_low <- normalRate_low_wholeGenome * sizeScale
normalRate_high <- normalRate_high_wholeGenome * sizeScale
tumourRate_high <- tumourRate_wholeGenome_high * sizeScale
tumourRate_low <- tumourRate_wholeGenome_low * sizeScale

### Generate LogNormal Distribution for Branch Rates
sdlog_normalRate <- (log(normalRate_high) - log(normalRate_low)) / (2 * z975)
mean_normalRate <- (log(normalRate_high) + log(normalRate_low)) / 2

sdlog_tumourRate <- (log(tumourRate_high) - log(tumourRate_low)) / (2 * z975)
mean_tumourRate <- (log(tumourRate_high) + log(tumourRate_low)) / 2

### Generate LogNormal Distribution for read count model parameters
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

# ADO rate
delta_high <- 0.3
delta_low <- 0.1
params_delta <- estimate_beta_params(delta_low, delta_high)

alpha_delta <- params_delta$alpha
beta_delta <- params_delta$beta

# t (coverage) and v (variance)
cov_high <- 6
cov_low <- 5

sdlog_cov <- (log(cov_high) - log(cov_low)) / (2 * z975)
mean_cov <- (log(cov_low) + log(cov_high)) / 2

v_high <- 25
v_low <- 15
sdlog_v <- (log(v_high) - log(v_low)) / (2 * z975)
mean_v <- (log(v_low) + log(v_high)) / 2

# epsilon (effective sequencing rate)
epsilon_high <- 0.0025
epsilon_low <- 0.00076

params_epsilon <- estimate_beta_params(epsilon_low, epsilon_high)

alpha_epsilon <- params_epsilon$alpha
beta_epsilon <- params_epsilon$beta

# overdispersion
w1_high <- 55
w1_low <- 50
mean_w1 <- (log(w1_low) + log(w1_high)) / 2
sdlog_w1 <- 0.33

w2_high <- 2.5
w2_low <- 2.55
mean_w2 <- (log(w2_low) + log(w2_high)) / 2
sdlog_w2 <- 0.33

### print the parameters
readCount_parameters <- sprintf(
  "delta ~ Beta(alpha = %.2f, beta = %.2f);
   t ~ LogNormal(meanlog = %.2f, sdlog = %.2f);
   v ~ LogNormal(meanlog = %.2f, sdlog = %.2f);
   epsilon ~ Beta(alpha = %.2f, beta = %.2f);
   w1 ~ LogNormal(meanlog = %.2f, sdlog = %.2f);
   w2 ~ LogNormal(meanlog = %.2f, sdlog = %.2f);",
  alpha_delta, beta_delta,
  mean_cov, sdlog_cov,
  mean_v, sdlog_v,
  alpha_epsilon, beta_epsilon,
  mean_w1, sdlog_w1,
  mean_w2, sdlog_w2
)

## Simplest Model (SM)
### The simplest model is the Yule tree contains 40 colorectal tumour cells only.
### Simulation
#### Lambda
## calculate the bonds for lambda
simpleTaxa <- 40 
lambda_high_SM <- sum(1 / (2:simpleTaxa)) / tumourTreeLow
lambda_low_SM <- sum(1 / (2:simpleTaxa)) / tumourTreeHigh

# calculate the parameters for log normal distribution
sdlog_lambda_SM <- (log(lambda_high_SM) - log(lambda_low_SM)) / (2 * z975)
mean_lambda_SM <- (log(lambda_low_SM) + log(lambda_high_SM)) / 2

#### Get Parameter Chunk
simulation_parameters_simple_model <- sprintf("lambda ~ LogNormal(meanlog = %.2f, sdlog = %.2f);
    mu ~ LogNormal(meanlog = %.2f, sdlog = %.2f);
", mean_lambda_SM, sdlog_lambda_SM, mean_tumourRate, sdlog_tumourRate)

### Replace Parameters
simple_model_template <- glue('
data {{
    RD = readFasta(file = \"/home/yil218/data/chr1.fna\");
    L = RD.nchar();
    n = 40;
}}
model{{
    snp ~ SNPSampler(alignment=RD, p=0.001);
    D = injectSNP(alignment=RD, snp=snp);

    {simulation_parameters_simple_model}

    ψ ~ Yule(lambda=lambda, n=n);

    // gt16 model
    π ~ Dirichlet(conc=[3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0]);
    rates ~ Dirichlet(conc=[1.0, 2.0, 1.0, 1.0, 2.0, 1.0]);
    Q = gt16(freq=π, rates=rates);
    A ~ SparsePhyloCTMC(L=L, Q=Q, tree=ψ, dataType=phasedGenotype(), mu=mu, root=D);

    conc = rep(8.0, 40);
    s ~ Dirichlet(conc=conc, mean=1);
    
    {readCount_parameters}
    
    alpha ~ Ploidy(l= L, n= n, delta= delta);
    cov ~ CoverageModel(alpha= alpha, t= t, v= v, s= s); // mimic CRC28 dataset
    r ~ ReadCountModel(D=A, epsilon=epsilon, alpha=alpha, coverage=cov, w1=w1, w2=w2);
}}')

writeLines(simple_model_template, "./simpleModel/1-simulation/sequenceGeneration.lphy")

## Calibrated Yule Model (CM)
## The CalibratedYule model is a tree with two calibration nodes: the root and the mrca of tumour cells. The tree contain 20 normal cells and 20 tumour cells.

### Simulation
#### Lambda
## calculate the bonds for lambda
calibratedTaxa <- 40 
lambda_high_CM <- sum(1 / (2:calibratedTaxa)) / normalTreeLow
lambda_low_CM <- sum(1 / (2:calibratedTaxa)) / normalTreeHigh

# calculate the parameters for log normal distribution
sdlog_lambda_CM <- (log(lambda_high_CM) - log(lambda_low_CM)) / (2 * z975)
mean_lambda_CM <- (log(lambda_low_CM) + log(lambda_high_CM)) / 2

#### Calibrations for Normal and Tumour Root
mean_normalRoot <- (log(normalTreeHigh) + log(normalTreeLow)) / 2
sdlog_normalRoot <- (log(normalTreeHigh) - log(normalTreeLow)) / (2 * z975)

mean_tumourRoot <- (log(tumourTreeHigh) + log(tumourTreeLow)) / 2
sdlog_tumourRoot <- (log(tumourTreeHigh) - log(tumourTreeLow)) / (2 * z975)

#### Print Lines
simulation_parameters_calibrated_model <- sprintf(
  "lambda ~ LogNormal(meanlog = %.2f, sdlog = %.2f);\n
    rootAge ~ LogNormal(meanlog = %.2f, sdlog = %.2f);\n
    tumourAge ~ LogNormal(meanlog = %.2f, sdlog = %.2f);\n
    rootRate ~ LogNormal(meanlog = %.2f, sdlog = %.2f);\n
    cladeRate ~ LogNormal(meanlog = %.2f, sdlog = %.2f);\n",
  mean_lambda_CM, sdlog_lambda_CM,
  mean_normalRoot, sdlog_normalRoot,
  mean_tumourRoot, sdlog_tumourRoot,
  mean_normalRate, sdlog_normalRate,
  mean_tumourRate, sdlog_tumourRate)

tumour_names_line <- paste0('"', paste0("clade_", 0:19), '"', collapse = ", ")

calibration_model_script <- glue('
data {{
    RD = readFasta(file = "/home/yil218/data/chr1.fna");
    L = RD.nchar();
    n = 40;
}}
model {{
    snp ~ SNPSampler(alignment=RD, p=0.001);
    D = injectSNP(alignment=RD, snp=snp);

    {simulation_parameters_calibrated_model}

    tumourNames = [{tumour_names_line}];

    tumourTaxa = taxa(names = tumourNames);

    tree ~ CalibratedYule(lambda=lambda, n=40, cladeTaxa=tumourTaxa, cladeMRCAAge=[tumourAge], rootAge=rootAge);

    tumourRoot = mrca(tree=tree, taxa=tumourNames);
    branchRates = localClock(tree=tree, clades=[tumourRoot], cladeRates=[cladeRate], rootRate=rootRate, includeStem=true);

    // gt16 model
    π ~ Dirichlet(conc=[3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0]);
    rates ~ Dirichlet(conc=[1.0, 2.0, 1.0, 1.0, 2.0, 1.0]);
    Q = gt16(freq=π, rates=rates);

    A ~ SparsePhyloCTMC(L=L, Q=Q, tree=tree, branchRates=branchRates, dataType=phasedGenotype(), root=D);
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
    RD = readFasta(file = "/home/yil218/data/chr1.fna");
    L = RD.nchar();
    n = 40;
}}
model {{
    snp ~ SNPSampler(alignment=RD, p=0.001);
    D = injectSNP(alignment=RD, snp=snp);

    {simulation_parameters_calibrated_model}

    tumourNames = [{tumour_names_line}];

    tumourTaxa = taxa(names = tumourNames);

    tree ~ CalibratedYule(lambda=lambda, n=40, cladeTaxa=tumourTaxa, cladeMRCAAge=[tumourAge], rootAge=rootAge);

    tumourRoot = mrca(tree=tree, taxa=tumourNames);
    branchRates = localClock(tree=tree, clades=[tumourRoot], cladeRates=[cladeRate], rootRate=rootRate, includeStem=true);

    // gt16 model
    π ~ Dirichlet(conc=[3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0]);
    rates ~ Dirichlet(conc=[1.0, 2.0, 1.0, 1.0, 2.0, 1.0]);
    Q = gt16(freq=π, rates=rates);

    A ~ SparsePhyloCTMC(L=L, Q=Q, tree=tree, branchRates=branchRates, dataType=phasedGenotype(), root=D);
    conc = rep(8.0, 40);
    s ~ Dirichlet(conc=conc, mean=1);
    
    {readCount_parameters}
    
    alpha ~ Ploidy(l= L, n= n, delta= delta);
    cov ~ CoverageModel(alpha= alpha, t= t, v= v, s= s); // mimic CRC28 dataset
    r ~ ReadCountModel(D=A, epsilon=epsilon, alpha=alpha, coverage=cov, w1=w1, w2=w2);
}}
')
  
  file_name <- sprintf("calibratedModel/1-simulation/lphy_scripts/sequenceGeneration_r%d.lphy", i)
  writeLines(script, file_name)
}

## Realistic Model (RM)
### Have a 10,000 taxa tree as healthy tree with normal ages and a 10,000 taxa tree as tumour tree with tumour age.
### Sample a branch to substitute on healthy tree with the tumour tree.
### Sample 20 normal cells and 20 tumour cells to get subsampled tree.
### Simulation
#### Lambda

## calculate the bonds for lambda
n <- 10000 # taxa number
lambda_high_RM <- sum(1 / (2:n)) / normalTreeLow
lambda_low_RM <- sum(1 / (2:n)) / normalTreeHigh

## calculate the parameters for log normal distribution
sdlog_lambda_RM <- (log(lambda_high_RM) - log(lambda_low_RM)) / (2 * z975)
mean_lambda_RM <- (log(lambda_low_RM) + log(lambda_high_RM)) / 2


#### Print Lines
simulation_parameters_realistic_model <- sprintf(
  "lambda ~ LogNormal(meanlog = %.2f, sdlog = %.2f);
tumourAge ~ LogNormal(meanlog = %.2f, sdlog = %.2f);
rootRate ~ LogNormal(meanlog = %.2f, sdlog = %.2f);
cladeRate ~ LogNormal(meanlog = %.2f, sdlog = %.2f);
",
  mean_lambda_RM, sdlog_lambda_RM,
  mean_tumourRoot, sdlog_tumourRoot,
  mean_normalRate, sdlog_normalRate,
  mean_tumourRate, sdlog_tumourRate
)

### Replace Parameters
realistic_model_template <- glue('
data {{
    RD = readFasta(file = \"/home/yil218/data/chr1.fna\");
    L = RD.nchar();
    n = 40;
    bigN = 10000;
}}
model{{
    snp ~ SNPSampler(alignment=RD, p=0.001);
    D = injectSNP(alignment=RD, snp=snp);

    {simulation_parameters_realistic_model}
    
    // healthy Yule tree
    ψhealthy ~ Yule(lambda=lambda, n=bigN);

    // sample a branch at the time
    node ~ SampleBranch(tree=ψhealthy, age=tumourAge);

    // tumour Yule tree with known root age
    ψtumour ~ Yule(lambda=lambda, n=bigN, rootAge=tumourAge);

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

    // sample the tree
    ψ ~ SubsampledTree(tree=ψfull, taxa=[scsTaxa], sampleFraction=[1.0]);
    sampledTumourRoot = mrca(tree=ψ , taxa=sampledTumourTaxa);
    ψlabelled = labelClade(tree=ψ, taxa=sampledTumourTaxa, label=nodeLabel);

    // gt16 model
    conc = rep(3.0, 16);
    π ~ Dirichlet(conc=conc);
    rates ~ Dirichlet(conc=[1.0, 2.0, 1.0, 1.0, 2.0, 1.0]);
    Q = gt16(freq=π, rates=rates);
    
     // generate branch rates
    branchRates = localClock(tree=ψ, clades=[sampledTumourRoot], cladeRates=[cladeRate], rootRate=rootRate, includeStem=true);


    A ~ SparsePhyloCTMC(L=L, Q=Q, tree=ψ, dataType=phasedGenotype(), branchRates = branchRates, root=D);

    conc = rep(8.0, 40);
    s ~ Dirichlet(conc=conc, mean=1);
    
    {readCount_parameters}
    
    alpha ~ Ploidy(l= L, n= n, delta= delta);
    cov ~ CoverageModel(alpha= alpha, t= t, v= v, s= s); // mimic CRC28 dataset
    r ~ ReadCountModel(D=A, epsilon=epsilon, alpha=alpha, coverage=cov, w1=w1, w2=w2);
}}')

writeLines(realistic_model_template, "./realisticModel/1-simulation/sequenceGeneration.lphy")
