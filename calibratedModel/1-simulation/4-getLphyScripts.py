import os
import shutil

# LPhy template
lphy_template = '''data{{
    treeList = readTrees(file="../trueTrees/tree_r{i}.trees");
    tree=treeList[0];
    RD = readFasta(file = "/home/yil218/data/chr1.fna");
    L = RD.nchar();
}}
model{{
    snp ~ SNPSampler(alignment=RD, p=0.001);
    D = injectSNP(alignment=RD, snp=snp);
    
    lambda ~ LogNormal(meanlog = -3.02, sdlog = 0.15);
    rootAge ~ LogNormal(meanlog = 4.21, sdlog = 0.15);
    tumourAge ~ LogNormal(meanlog = 1.84, sdlog = 0.12);
    rootRate ~ LogNormal(meanlog = -9.55, sdlog = 0.06);
    cladeRate ~ LogNormal(meanlog = -7.43, sdlog = 0.12);

    tumourNames = ["clade_0","clade_1","clade_2","clade_3","clade_4",
    "clade_5","clade_6","clade_7","clade_8","clade_9",
    "clade_10","clade_11","clade_12","clade_13","clade_14",
    "clade_15","clade_16","clade_17","clade_18","clade_19"];

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
    
    delta ~ LogNormal(meanlog = -1.96, sdlog = 0.18); //ADO rate 0.1-0.2
    t ~ LogNormal(meanlog = 1.70, sdlog = 0.05); //coverage 5-6
    v ~ LogNormal(meanlog = 2.96, sdlog = 0.13); // variance 15-25
    epsilon ~ LogNormal(meanlog = -6.59, sdlog = 0.30); // effective sequencing rate 7.6*10 -4 ~ 2.5*10-3
    w1 ~ LogNormal(meanlog = 3.96, sdlog = 0.33); // wildtype overdispersion 50-55
    w2 ~ LogNormal(meanlog = 0.93, sdlog = 0.33); // alternative overdispersion 2.5-2.55

    alpha ~ Ploidy(l= L, n= 40, delta= delta);
    cov ~ CoverageModel(alpha= alpha, t= t, v= v, s= s); // mimic CRC28 dataset
    r ~ ReadCountModel(D=A, epsilon=epsilon, alpha=alpha, coverage=cov, w1=w1, w2=w2);
}}
'''

# set the path
output_dir = "lphy_scripts"  # The directory to save the generated scripts

# Check and create the output directory
if os.path.exists(output_dir):
    shutil.rmtree(output_dir)  # Remove the directory and its contents
    print(f"The directory {output_dir} and all its contents were deleted.")

os.makedirs(output_dir)  # Create the output directory
print(f"The directory {output_dir} was created.")

# Process all tree files
for i in range(100):  # r0-r99
    # Save the generated script to the output folder
    file_name = f"{output_dir}/sequenceGeneration_r{i}.lphy"
    lphy_script = lphy_template.format(i=i)

    with open(file_name, "w") as f:
        f.write(lphy_script)
        print(f"Generated: {file_name}")

print("All LPhy scripts generated successfully.")
