# Variational_Inference_Of_PRS_Using_Empirical_Prior
In continuation of project: https://github.com/Sagarnandeshwar/Variational_Inference_Of_Bayesian_Linear_Regression
In this project, we will explore the usefulness of empirical prior in PRS i.e. estimated from real-world data rather than being predefined based on theoretical assumptions. Empirical priors can capture the actual distribution of genetic effects observed in the population, allowing for more accurate and data-driven inference 

## Motivation: 
Large-scale genome-wide association studies (GWASs) has driven the development of statistical methods for phenotype prediction. One of such techniques is the polygenic risk score (PRS) methods that formulate the task of polygenic prediction in terms of a multiple linear regression framework, where the goal is to infer the joint effect sizes of all genetic variants on the trait. These methods are used to estimate an individual’s genetic risk for a particular trait or disease based on the combined effect of multiple genetic variants, hence they are a valuable tool for predicting genetic risk and understanding the genetic basis of complex traits and diseases. 

## Challenges: 
PRS methods when paired with modern GWAS sample sizes, which consists of hundreds of thousands of individuals, with high dimensional data possess several computational and statistical challenges. The closed-form update equations for some of the variational parameters involve terms that relate to the LD between the focal variant and all other variants in the genome. This computationally prohibitive to compute for millions of variants and for hundreds of EM iterations. Furthermore, most individual-level GWAS data sources are protected for privacy concerns. 

## Related Work 
One of the common methods in these analyses is single nucleotide polymorphisms (SNPs), measured by either genotyping arrays or imputed using reference haplotypes. Byesian models are other PRS models that incorporate prior knowledge such as probability distributions over the genetic causal architecture of complex traits. However, a major limitation of Bayesian methods is that their scalability is extremely slow and inefficient with its inference techniques. 

## Method  
### Overview 
In this project, I implemented Variational inference of polygenic risk score (VIPRS) which is an is an efficient Bayesian inference algorithm utilizes variational inference to approximate the posterior for the effect sizes for high-dimensional multiple regression applied in the context of polygenic risk score (PRS) problem. Here, I have implemented the updates of point estimates for the annotation weights w and the point estimates for the variance explained per annotation k, updating the Empirical Bayes algorithm (Expectation-Maximization (EM)). 

### Theory: 
For P SNPs and and N individuals, let X denotes P x N genotype matrix and y ∈ $R^{N×1}$ the phenotype vector. Let $a_{j,k}$ ∈ {0, 1}
be the binary indicator for whether SNP j is in annotation k, i.e., j ∈ $C_k$, where $C_k$ is the set of SNPs that are in annotation k. Therefore, for K annotations, we have a P x K binary annotation matrix A.

Let w ∈ $R^{K×1}$ be the unknown annotation weights, $σ_k^2$ the variance explained by SNP j if it is in annotation k, and β and s be the effect size and binary indicator variables of causal SNPs, respectively.
Together, we have the following data generative process.

![theory](https://github.com/Sagarnandeshwar/Variational_Inference_Of_PRS_Using_Empirical_Prior/blob/main/images/model.png)


### Algorithm:  
I modify the the codes from https://github.com/Sagarnandeshwar/Variational_Inference_Of_Bayesian_Linear_Regression, by implementing by implementing the updates of point estimates for the annotation weights w ∈ $R^{K×1}$ and the point estimates for the variance explained per annotation k: $σ_k^2$ = $τ_k^{−1}$ for k ∈ {1, . . . ,K}. To achieve that, you will follow the Expectation-Maximization (EM) algorithm (aka Empirical Bayes) outlined in Algorithm 1. The algorithm is similar to the original VIPRS you implemented in A3 but with the changes highlighted in red to incorporate the functional priors. Equations (6) and (7) were derived by taking the partial derivative of the ELBO w.r.t. $τ_k$ and w, respectively. Make sure you know how to derive them. Because of the logistic function, there is no closed-form update for w. Instead, we perform gradient ascent in (7) with some fixed learning rate η. Equation (5) is optional. You may fix it to 1 if you have trouble of convergence.

![em](https://github.com/Sagarnandeshwar/Variational_Inference_Of_PRS_Using_Empirical_Prior/blob/main/images/em.png)

## Dataset: 

**Linkage-Disequilibrium matrices** The LD matrices can be computed from a reference panel, such as the 1000 Genome dataset (ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/) or you can download pre-computed LD matrices from the UK Biobank dataset from Zenodo: https://zenodo.org/record/7036625#.ZBeDiOzMJQI. 

**Annotations** The ‘annotations’ folder stores the annotations for SNPs on chromosome 21 and 22. In total, there are 101 columns and the last 96 columns are the annotations. The above algorithm works only for binary annotations. Some annotations are continuous (e.g., GERP.NS, MAF Adj LLD AFR, Recomb Rate 10kb, Nucleotide Diversity 10kb, etc). 

**Summary statistics data** The training and testing summary statistics data computed for human standing height on chromosome 21 and 22 from 337,205 White British individuals from the UK Biobank. The data are divided into 5 folds to facilitate evaluation described below. 

## Data Processing 
I use magenpy, a python library for loading, manipulating, and simulating with genotype data. I created GWADataLoader object for training and testing data per fold, and then creates several numpy arrays to store per SNPs data and per annotations data. The use of numpy methods has significantly speed up the algorithm. I have also limits the γ value it 0.01, 0.99), in order to stabalize the algorithm. 

## Simulation 
For each of the two LD matrix (ch21 and ch21) For each of the five fold
1. Use mergepy GWADataLoader to
   a Load the LD matrix
   (b) Load the Training marginal beta
   (c) Load the Annotation matrix
   mergy simultaneous harmonize data all the three dataset
2. Use numpy to create array for
   - per snps data
      - $μ_{β_j}^*$
      - $τ_{β_j}^*$
      - $τ_β$
      - π
      - $γ_j^*$
   - per annotation data
     - $τ_k$
     - w
3. Run the EM algorithm for 10 iterations with
   - Learning rate = 0.001
   - $τ_ϵ$ = 1,0
4. Use mergepy to load testing marginal beta
5. match the testing Marginal beta with training Marginal beta
6. Compute the R square
 
## Evaluation metric 
### R-squared
We use R-squared to evaluate how your trained model performs on the corresponding testing fold.  

R-squared measures the proportion of the variance in the dependent variable that is explained by the independent variable(s) included in the model. R-squared values range from 0 to 1, with a value of 1 indicating that the model explains 100% of the variation in the dependent variable and a value of 0 indicating that the model explains none of the variation. 

Assuming standardized phenotype, the R-squared on testing fold can be computed by summary statistics: 
![r2](https://github.com/Sagarnandeshwar/Variational_Inference_Of_PRS_Using_Empirical_Prior/blob/main/images/r2.png)

We performed 5-fold cross-validation by training and evaluating your model on the training and testing fold with the same index and repeat your experiments 5 times to compute standard error of the R-squared estimates for each method. 

### Evidence lower bound 

The evidence lower bound (ELBO) of the model is:  
![elbo](https://github.com/Sagarnandeshwar/Variational_Inference_Of_PRS_Using_Empirical_Prior/blob/main/images/elbo.png)
More specifically,  
![elbo](https://github.com/Sagarnandeshwar/Variational_Inference_Of_PRS_Using_Empirical_Prior/blob/main/images/elbos.png)
where $γ_j^*$ , $μ_j^*$ and $τ_j$ are the inferred PIP, mean and precision of the effect size for SNP j at the E-step, respectively; ◦ is the elementwise product of two vectors. 

I have also implemented the LELBO for the model, but due to computational limitation, I am unable to run it for the whole data. 

## Result 
![result](https://github.com/Sagarnandeshwar/Variational_Inference_Of_PRS_Using_Empirical_Prior/blob/main/images/1.png)
![result](https://github.com/Sagarnandeshwar/Variational_Inference_Of_PRS_Using_Empirical_Prior/blob/main/images/2.png)
![result](https://github.com/Sagarnandeshwar/Variational_Inference_Of_PRS_Using_Empirical_Prior/blob/main/images/3.png)
I was able to implement the EM algorithm for VIPRS and performed 5-fold cross-validation by training and evaluating your model on the training and testing fold by computing standard error of the R-squared estimates for each method for both chromosomes ch21 and ch22 sets in an efficient way with the help of magenpy and numpy. Within the same chromosomes group the r square values are very similar, however we see that the model performs slightly better when trained with LD of chromosomes ch21 than ch22. All the R square values observed are low. This may be due to small sample space, biased initial values or because of the hyperparameter selection. 

## Discussion 
The method can be affected from population stratification, which occurs when there are systematic differences in allele frequencies between different populations. This can lead to spurious associations between genetic variants and traits or diseases. 

Secondly, the method assumes that there is a linear relationship between the effect of individual genetic variants and the risk for the trait or disease. The actual relationship may be non-linear or depend on environmental factors or different genetic variants interactions. 

The PRS methods heavily relies on genome-wide association studies (GWAS) to identify genetic variants associated with a particular trait or disease. However, GWAS studies often have limited sample sizes, which can lead to false positives or miss important genetic variants that contribute to the trait or disease. 

## Future work 
We could try to run the model with different hyperparameters and initial values, and could also implement Evidence lower bound to observe the models performances over the time of training. 

Second, the current spike-and-slab prior assumes that all genetic variants have a uniform prior probability of being causal and that the causal SNPs have equal expected contribution to the heritability. For the future work we could explore a more general and flexible Gaussian mixture prior. 

Lastly, we could try to joint the model effect sizes from multiple ancestrally homogeneous populations within the same framework and evaluate the performance across different populations. 

## Reference 
”Fast and accurate Bayesian polygenic risk modeling with variational inference”(2022) - Shadi Zabad, Simon Gravel, and Yue Li 
