Case study: variable selection for the approximation of body fat
================
Lacey W. Heinsberg
November 06, 2024, 15:18







The purpose of this .Rmd is to perform variable selection and
multivariable modeling of body fat percentage.

Goal: Approximate the proportion of body fat that can be estimated by
simple anthropometric measures.

This file uses source code that accompanies section 3.3 of the
manuscript [“Variable selection: A review and recommendations for the
practicing
statistician”](https://onlinelibrary.wiley.com/doi/full/10.1002/bimj.201700067)
that was created by Georg Heinze, Christine Wallisch, and Daniela
Dunkler. This code has been modified and annotated for:

Heinsberg LW. Statistical Tools to Support Implementation: Variable
Selection and Post-Selection Inference in Genomic Nursing Research.
\[Expert Lecturer Abstract, Podium\]. Presented at the International
Society of Nurses in Genetics, November 2024, San Diego, California.

For questions, comments, or remarks, feel free to contact Lacey
Heinsberg (<law145@pitt.edu>).

# Load libraries

``` r
library(shrink)
library(pander)
library(corrplot)
library(dplyr)
```

# Load and prepare data

The data set used for this example are from [Johnson’s
(1996)](https://www.tandfonline.com/doi/full/10.1080/10691898.1996.11910505)
body fat study which was intended as an educational data set to teach
multivariable linear regression in the classroom. The original purpose
of the data set was to support the approximation of a costly measurement
of body density (from which proportion of body fat can be derived with
[Siri’s
(1956)](https://www.sciencedirect.com/science/article/abs/pii/B978148323110550011X)
formula) by a combination of age, height, weight, and ten simple
anthropometric circumference measures through multivariable linear
regression.

A more detailed description of the data can be found at:
<https://ww2.amstat.org/publications/jse/datasets/fat.txt>

In line with other literature, Heinze, Wallisch, and Dunkler excluded
one individual from the original data set with an implausible height
observation. The team also converted the units of some variables for the
approach detailed below.

``` r
# Load data ------------------------------------------------------
case1.bodyfat <- read.table("data/case1_bodyfat.txt", header = T, sep = ";")
n <- nrow(case1.bodyfat)
```

Before using the data, we are going to modify it a bit by adding some
simulated genetic data and nursing-centric covariates.

``` r
# Set seed for reproducibility
set.seed(123)

# Simulate SNP data (rs1 to rs6) with varying minor allele frequencies (MAFs)
case1.bodyfat <- case1.bodyfat %>%
  mutate(
    rs1 = rbinom(n, 2, 0.3),  # rs1 with MAF of 0.3
    rs2 = rbinom(n, 2, 0.4),  # rs2 with MAF of 0.4
    rs3 = rbinom(n, 2, 0.2),  # rs3 with MAF of 0.2
    rs4 = rbinom(n, 2, 0.5),  # rs4 with MAF of 0.5
    rs5 = rbinom(n, 2, 0.3),  # rs5 with MAF of 0.3
    rs6 = rbinom(n, 2, 0.4)   # rs6 with MAF of 0.4
  )

# Simulate social drivers of health (SDOH) and behavioral factorss
case1.bodyfat <- case1.bodyfat %>%
  mutate(
    num_children = pmax(0, rpois(n, 2)),  # Number of children, minimum 0
    family_mealtime_freq = pmin(21, pmax(0, round(rlnorm(n, 2, 0.5)))),  # Family mealtime frequency, 0 and 21 meals per week
    physical_environment_score = rnorm(n, 50, 10), # Physical environment score (unitless)
    stress_level = pmin(10, pmax(0, rnorm(n, 5, 2))), # Stress, 0-10 with higher values indicating higher stress
    physical_activity = pmin(16, pmax(0, rnorm(n, 1, 0.5))), # Physical activity in hours per day, 0-16 hours per day
    sleep_duration = pmin(16, pmax(0, rnorm(n, 9, 2))), # Sleep duration in hours 0 to 16 hours
    healthy_eating_index = pmin(100, pmax(0, rnorm(n, 50, 10))) # Healthy eating index, 0-100 (unitless)
  )

# Create a GxE interaction effect (for added interest) 
case1.bodyfat <- case1.bodyfat %>%
  mutate(
    gxe_effect = rs1 * sleep_duration  # Interaction term between rs1 and sleep_duration
  )

# Create a new outcome that preserves the original siri structure with 
# added SNP and interaction effects (so we can pretend we are onto something BIG!)
case1.bodyfat <- case1.bodyfat %>%
  mutate(
    siri_simulated = siri + 1.5 * rs1 - 1.2 * rs2 - 0.4 * stress_level +
                    0.5 * physical_activity - 0.3 * (sleep_duration - 7)^2 + 
                    1.2 * gxe_effect 
  )

# Double the dataset for a larger sample size
duplicated_data <- case1.bodyfat[sample(1:n, n, replace = TRUE), ]
case1.bodyfat <- rbind(case1.bodyfat, duplicated_data)

# Check the size of the new dataset and preview the first few rows
n <- nrow(case1.bodyfat)
head(case1.bodyfat)
```

    ##   case brozek siri density age weight weight_kg height height_cm neck chest
    ## 1    1   12.6 12.3  1.0708  23 154.25      70.1  67.75       172 36.2  93.1
    ## 2    2    6.9  6.1  1.0853  22 173.25      78.8  72.25       184 38.5  93.6
    ## 3    3   24.6 25.3  1.0414  22 154.00      70.0  66.25       168 34.0  95.8
    ## 4    4   10.9 10.4  1.0751  26 184.75      84.0  72.25       184 37.4 101.8
    ## 5    5   27.8 28.7  1.0340  24 184.25      83.8  71.25       181 34.4  97.3
    ## 6    6   20.6 20.9  1.0502  24 210.25      95.6  74.75       190 39.0 104.5
    ##   abdomen   hip thigh knee ankle biceps forearm wrist rs1 rs2 rs3 rs4 rs5 rs6
    ## 1    85.2  94.5  59.0 37.3  21.9   32.0    27.4  17.1   0   0   0   0   1   2
    ## 2    83.0  98.7  58.7 37.3  23.4   30.5    28.9  18.2   1   1   0   1   0   0
    ## 3    87.9  99.2  59.6 38.9  24.0   28.8    25.2  16.6   0   0   0   1   1   0
    ## 4    86.4 101.2  60.1 37.3  22.8   32.4    29.4  18.2   1   1   0   1   0   2
    ## 5   100.0 101.9  63.2 42.2  24.0   32.2    27.7  17.7   2   1   0   2   0   1
    ## 6    94.4 107.8  66.0 42.0  25.6   35.7    30.6  18.8   0   0   0   0   0   0
    ##   num_children family_mealtime_freq physical_environment_score stress_level
    ## 1            2                    7                   70.95705     3.884423
    ## 2            2                   16                   60.90176     4.363083
    ## 3            2                    3                   38.54626     4.558916
    ## 4            2                    7                   29.45192     3.840704
    ## 5            2                   21                   62.58065     5.377836
    ## 6            5                    6                   43.52461     4.297470
    ##   physical_activity sleep_duration healthy_eating_index gxe_effect
    ## 1         1.2085865      11.086203             49.11380    0.00000
    ## 2         1.9078413      14.081942             65.61617   14.08194
    ## 3         0.7157718       7.876441             59.66893    0.00000
    ## 4         0.7950628      12.488092             49.49942   12.48809
    ## 5         0.7437824       6.717081             34.27745   13.43416
    ## 6         0.7138788      10.818667             46.77544    0.00000
    ##   siri_simulated
    ## 1       6.341407
    ## 2       7.460846
    ## 3      23.603875
    ## 4      15.511214
    ## 5      44.817738
    ## 6      15.163286

Our example data set consists of 502 participants with data on body fat
percentage (variable name `siri_simulated`, our outcome of interest) and
several candidate variables hypothesized to be related to body fat
percentage.

# Define study question

Objective: The study aims to examine how genetic variation
(obesity-related SNPs), anthropometric measures, social drivers of
health (SDOH), and behavioral factors influence a continuous measure of
fat mass in a sample of adults aged 22-81 years.

Population: A sample of 502 adults, aged 22-81 years.

Outcome Variable: \* siri_simulated (continuous)

Explanatory Variables:

1.  Genetics

- rs1 to rs6: Obesity-related genetic markers (0 = no risk alleles, 1 =
  one risk allele, 2 = two risk alleles).
- Gene x Environment Effect: Interaction effect of rs1 and sleep
  duration.

2.  Demographics

- Age: Continuous (years)

3.  Anthropometrics

- Weight: Continuous (kg).
- Height: Continuous (cm).
- Circumference Values: \*\* Neck, abdomen, thigh, knee, ankle, biceps,
  forearm, wrist: Continuous (cm).

4.  Family dynamics

- Number of Children in Household: Continuous, the total number of
  children living in the household.
- Family Mealtime Frequency: Continuous, number of family meals per week
  (0–21, right-skewed).

5.  Environmental factors

- Physical Environment Score: Composite continuous score, higher scores
  indicate a safer and more exercise-friendly environment.

6.  Psychosocial factors

- Stress Level: Continuous, with higher scores indicating more stress.

7.  Behavioral factors

- Physical Activity Level: Continuous, hours per day spent on physical
  activity.
- Sleep Duration: Continuous, hours of sleep per night.
- Healthy Eating Index: Continuous, with higher scores indicating better
  eating habits.

# Analyses

## Events-per-variable (EPV)

The ratio between sample size and the number of independent variables is
termed “events-per-variable” (EPV). EPV quantifies the balance between
the amount of information provided by the data and the number of unknown
parameters that should be estimated.

With a limited sample size, it is not possible to accurately estimate
many regression coefficients. Therefore, EPV “rules of thumb” (such as
that recommend in Table 3 of [Heinze, Wallisch, and Dunkler
(2018)](https://onlinelibrary.wiley.com/doi/full/10.1002/bimj.201700067))
can be considered when developing your analytical approach.

For an expanded commentary on EPV, see section 2.1.3 of [Heinze,
Wallisch, and Dunkler
(2018)](https://onlinelibrary.wiley.com/doi/full/10.1002/bimj.201700067).

Let’s calculate the EPV for our example data set.

``` r
# EPV --------------------------------------------------------------
# Store predictors of interest
pred <- c("age", "weight_kg", "height_cm", "neck",  "abdomen",  
          "thigh", "knee", "ankle", "biceps", "forearm", "wrist", "rs1", "rs2", "rs3", "rs4", "rs5", "rs6", "num_children", "family_mealtime_freq", "physical_environment_score", "stress_level", "physical_activity", "sleep_duration", "healthy_eating_index", "gxe_effect" )

# Compute EPV 
epv <- dim(case1.bodyfat)[1] / length(pred)
epv
```

    ## [1] 20.08

In our example data set, we have an acceptable EPV of 20.08.

## Correlation structure

Note that the EPV ratio above often oversimplifies the analytical
approach because - beyond the number of independent variables - many
other quantities such as the correlation structure of a data set may
influence accuracy [(Courvoisier, Combescure, Agoritsas, Gayet-Ageron, &
Perneger, 2011)](https://pubmed.ncbi.nlm.nih.gov/21411281/). The
recommended EPV limits should be adapted to the situation (e.g., raised
if correlations between candidate independent variables are particularly
strong; lowered if candidate variables are all independent of each
other) [(Heinze, Wallisch, and Dunkler,
2018)](https://onlinelibrary.wiley.com/doi/full/10.1002/bimj.201700067).

Further, strong correlations impose some challenges in model development
and interpretation. In particular, interpretation of regression
coefficients as adjusted effects in the global model, or, if variable
selection is applied, interpretation of non-selected variables as
“non-predictive” can be problematic. As such, understanding the
correlation structure of our data set is very important.

Let’s examine it now.

``` r
# Correlation structure --------------------------------------------
pander(cor(case1.bodyfat[, pred]))
```

|                                |    age    | weight_kg | height_cm |   neck    |
|:------------------------------:|:---------:|:---------:|:---------:|:---------:|
|            **age**             |     1     | -0.01247  |  -0.2645  |  0.1228   |
|         **weight_kg**          | -0.01247  |     1     |  0.4687   |  0.8046   |
|         **height_cm**          |  -0.2645  |  0.4687   |     1     |   0.264   |
|            **neck**            |  0.1228   |  0.8046   |   0.264   |     1     |
|          **abdomen**           |  0.2343   |  0.8839   |  0.1426   |  0.7328   |
|           **thigh**            |  -0.1779  |   0.863   |  0.3357   |  0.6487   |
|            **knee**            | 0.008152  |  0.8352   |  0.5158   |   0.629   |
|           **ankle**            |  -0.1235  |  0.6609   |  0.4197   |  0.4869   |
|           **biceps**           | -0.01943  |  0.8078   |  0.2773   |  0.7315   |
|          **forearm**           |  -0.1262  |  0.6209   |  0.2644   |  0.6174   |
|           **wrist**            |  0.2168   |  0.7238   |  0.3819   |  0.7346   |
|            **rs1**             | -0.01434  | 0.007234  |  0.09656  |  0.02961  |
|            **rs2**             |  -0.0263  |  -0.041   |  0.04244  | -0.009341 |
|            **rs3**             |  0.1448   |  -0.1237  |  -0.1345  | -0.03345  |
|            **rs4**             | -0.009383 | -0.02824  |  0.04312  | 0.004635  |
|            **rs5**             | -0.01368  |  -0.1076  | 0.003148  |  -0.0907  |
|            **rs6**             |  0.06632  | -0.002413 | -0.07218  |  0.04036  |
|        **num_children**        |  0.01667  | 0.008047  |  0.05035  | -0.07331  |
|    **family_mealtime_freq**    |  0.0355   | -0.02631  | -0.03879  |  0.05389  |
| **physical_environment_score** | -0.02363  | -0.04044  |  -0.1022  | -0.001796 |
|        **stress_level**        | -0.06532  | -0.08934  |  -0.147   | -0.08747  |
|     **physical_activity**      |  -0.1253  |  0.1012   |  0.06056  |  0.06895  |
|       **sleep_duration**       | -0.04537  |  0.05496  | -0.009391 |  0.1458   |
|    **healthy_eating_index**    |  0.01366  | -0.001971 | 0.009323  | -0.008544 |
|         **gxe_effect**         | -0.03343  |  0.03397  |  0.1271   |  0.0835   |

Table continues below

|                                |  abdomen  |   thigh   |   knee    |   ankle   |
|:------------------------------:|:---------:|:---------:|:---------:|:---------:|
|            **age**             |  0.2343   |  -0.1779  | 0.008152  |  -0.1235  |
|         **weight_kg**          |  0.8839   |   0.863   |  0.8352   |  0.6609   |
|         **height_cm**          |  0.1426   |  0.3357   |  0.5158   |  0.4197   |
|            **neck**            |  0.7328   |  0.6487   |   0.629   |  0.4869   |
|          **abdomen**           |     1     |  0.7507   |  0.6875   |  0.4838   |
|           **thigh**            |  0.7507   |     1     |  0.7941   |  0.5947   |
|            **knee**            |  0.6875   |  0.7941   |     1     |  0.6623   |
|           **ankle**            |  0.4838   |  0.5947   |  0.6623   |     1     |
|           **biceps**           |  0.6883   |  0.7662   |  0.6667   |  0.5101   |
|          **forearm**           |   0.476   |  0.5393   |  0.5013   |  0.4396   |
|           **wrist**            |  0.5975   |  0.5608   |  0.6605   |  0.6013   |
|            **rs1**             | -0.01009  |  0.03311  | 0.008891  | -0.01349  |
|            **rs2**             | -0.08265  | -0.04418  | -0.03461  |  -0.0314  |
|            **rs3**             | -0.07015  |  -0.1288  | -0.08027  |  -0.1296  |
|            **rs4**             | -0.03774  | -0.004073 | -0.005299 | -0.06772  |
|            **rs5**             |  -0.1226  |  -0.0847  | -0.08809  |  -0.1073  |
|            **rs6**             |  0.02775  | -0.03063  | -0.02925  | -0.04941  |
|        **num_children**        | 0.0008927 | -0.05564  |  0.04627  |   0.122   |
|    **family_mealtime_freq**    | -0.02505  | -0.06108  | 0.003448  | -0.05825  |
| **physical_environment_score** | -0.02332  | -0.04462  |  -0.1356  | -0.05918  |
|        **stress_level**        | -0.04379  | -0.06745  |  -0.1395  | -0.07937  |
|     **physical_activity**      |  0.05438  |  0.0898   |  0.03015  |  0.0824   |
|       **sleep_duration**       |  0.07987  |  0.01817  | -0.01762  |  -0.0371  |
|    **healthy_eating_index**    | -0.04941  | -0.01847  |  0.01922  |  0.01243  |
|         **gxe_effect**         |  0.01293  |  0.03944  |  0.01225  | 0.0009521 |

Table continues below

|                                |  biceps  |  forearm  |   wrist   |   rs1    |
|:------------------------------:|:--------:|:---------:|:---------:|:--------:|
|            **age**             | -0.01943 |  -0.1262  |  0.2168   | -0.01434 |
|         **weight_kg**          |  0.8078  |  0.6209   |  0.7238   | 0.007234 |
|         **height_cm**          |  0.2773  |  0.2644   |  0.3819   | 0.09656  |
|            **neck**            |  0.7315  |  0.6174   |  0.7346   | 0.02961  |
|          **abdomen**           |  0.6883  |   0.476   |  0.5975   | -0.01009 |
|           **thigh**            |  0.7662  |  0.5393   |  0.5608   | 0.03311  |
|            **knee**            |  0.6667  |  0.5013   |  0.6605   | 0.008891 |
|           **ankle**            |  0.5101  |  0.4396   |  0.6013   | -0.01349 |
|           **biceps**           |    1     |  0.6434   |  0.6338   | 0.01249  |
|          **forearm**           |  0.6434  |     1     |   0.525   | -0.03855 |
|           **wrist**            |  0.6338  |   0.525   |     1     | -0.05724 |
|            **rs1**             | 0.01249  | -0.03855  | -0.05724  |    1     |
|            **rs2**             | -0.03105 |  -0.1358  | -0.002128 | 0.08062  |
|            **rs3**             | -0.1113  |  -0.1008  | -0.03533  | 0.003719 |
|            **rs4**             | 0.05571  |  0.02124  |  0.03285  |  0.1287  |
|            **rs5**             | -0.07457 |  -0.1288  | -0.05387  | 0.08475  |
|            **rs6**             | -0.01449 | 0.006098  | -0.01727  | 0.03585  |
|        **num_children**        | -0.03495 | -0.008623 |  0.04174  | -0.1206  |
|    **family_mealtime_freq**    | -0.00216 | -0.02213  |  0.01035  | 0.03253  |
| **physical_environment_score** | 0.01088  | -0.01056  | -0.07075  | 0.02125  |
|        **stress_level**        | -0.06364 |  0.02742  |  -0.1021  | -0.1386  |
|     **physical_activity**      |  0.1093  |  0.1041   |  0.08991  | -0.1103  |
|       **sleep_duration**       | 0.06898  |  0.1191   | -0.009479 | 0.08817  |
|    **healthy_eating_index**    | 0.02064  |  0.01588  |  0.08796  | -0.08155 |
|         **gxe_effect**         | 0.04156  | 0.0005114 | -0.03647  |  0.9521  |

Table continues below

|                                |    rs2    |    rs3    |    rs4    |   rs5    |
|:------------------------------:|:---------:|:---------:|:---------:|:--------:|
|            **age**             |  -0.0263  |  0.1448   | -0.009383 | -0.01368 |
|         **weight_kg**          |  -0.041   |  -0.1237  | -0.02824  | -0.1076  |
|         **height_cm**          |  0.04244  |  -0.1345  |  0.04312  | 0.003148 |
|            **neck**            | -0.009341 | -0.03345  | 0.004635  | -0.0907  |
|          **abdomen**           | -0.08265  | -0.07015  | -0.03774  | -0.1226  |
|           **thigh**            | -0.04418  |  -0.1288  | -0.004073 | -0.0847  |
|            **knee**            | -0.03461  | -0.08027  | -0.005299 | -0.08809 |
|           **ankle**            |  -0.0314  |  -0.1296  | -0.06772  | -0.1073  |
|           **biceps**           | -0.03105  |  -0.1113  |  0.05571  | -0.07457 |
|          **forearm**           |  -0.1358  |  -0.1008  |  0.02124  | -0.1288  |
|           **wrist**            | -0.002128 | -0.03533  |  0.03285  | -0.05387 |
|            **rs1**             |  0.08062  | 0.003719  |  0.1287   | 0.08475  |
|            **rs2**             |     1     |  0.1411   |  0.02408  | 0.006181 |
|            **rs3**             |  0.1411   |     1     | 0.001738  | 0.04025  |
|            **rs4**             |  0.02408  | 0.001738  |     1     | 0.07163  |
|            **rs5**             | 0.006181  |  0.04025  |  0.07163  |    1     |
|            **rs6**             |  0.02182  |  -0.0139  | -0.05527  | -0.01235 |
|        **num_children**        | -0.05755  | -0.08919  |  0.01325  | -0.03335 |
|    **family_mealtime_freq**    |  -0.0412  | -0.04675  |  -0.1007  | -0.0598  |
| **physical_environment_score** |  0.06366  | -0.04117  | -0.07579  | -0.05139 |
|        **stress_level**        |  0.1337   |  0.02588  |  0.1737   | -0.08746 |
|     **physical_activity**      | -0.03517  |  -0.1039  |  0.0457   | -0.1049  |
|       **sleep_duration**       |  -0.142   |  0.01076  |   0.106   | -0.1257  |
|    **healthy_eating_index**    |  0.04107  |  0.1098   |  0.0307   | -0.09141 |
|         **gxe_effect**         |  0.03138  | -0.008795 |  0.1739   | 0.04818  |

Table continues below

|                                |    rs6    | num_children |
|:------------------------------:|:---------:|:------------:|
|            **age**             |  0.06632  |   0.01667    |
|         **weight_kg**          | -0.002413 |   0.008047   |
|         **height_cm**          | -0.07218  |   0.05035    |
|            **neck**            |  0.04036  |   -0.07331   |
|          **abdomen**           |  0.02775  |  0.0008927   |
|           **thigh**            | -0.03063  |   -0.05564   |
|            **knee**            | -0.02925  |   0.04627    |
|           **ankle**            | -0.04941  |    0.122     |
|           **biceps**           | -0.01449  |   -0.03495   |
|          **forearm**           | 0.006098  |  -0.008623   |
|           **wrist**            | -0.01727  |   0.04174    |
|            **rs1**             |  0.03585  |   -0.1206    |
|            **rs2**             |  0.02182  |   -0.05755   |
|            **rs3**             |  -0.0139  |   -0.08919   |
|            **rs4**             | -0.05527  |   0.01325    |
|            **rs5**             | -0.01235  |   -0.03335   |
|            **rs6**             |     1     |   0.06119    |
|        **num_children**        |  0.06119  |      1       |
|    **family_mealtime_freq**    |  0.1483   |   -0.0524    |
| **physical_environment_score** |  0.05591  |    -0.13     |
|        **stress_level**        |  0.02225  |   -0.07869   |
|     **physical_activity**      |  0.01894  |   0.02586    |
|       **sleep_duration**       | -0.003312 |   0.07246    |
|    **healthy_eating_index**    |  0.0138   |   0.01976    |
|         **gxe_effect**         |  0.02395  |   -0.1018    |

Table continues below

|                                | family_mealtime_freq |
|:------------------------------:|:--------------------:|
|            **age**             |        0.0355        |
|         **weight_kg**          |       -0.02631       |
|         **height_cm**          |       -0.03879       |
|            **neck**            |       0.05389        |
|          **abdomen**           |       -0.02505       |
|           **thigh**            |       -0.06108       |
|            **knee**            |       0.003448       |
|           **ankle**            |       -0.05825       |
|           **biceps**           |       -0.00216       |
|          **forearm**           |       -0.02213       |
|           **wrist**            |       0.01035        |
|            **rs1**             |       0.03253        |
|            **rs2**             |       -0.0412        |
|            **rs3**             |       -0.04675       |
|            **rs4**             |       -0.1007        |
|            **rs5**             |       -0.0598        |
|            **rs6**             |        0.1483        |
|        **num_children**        |       -0.0524        |
|    **family_mealtime_freq**    |          1           |
| **physical_environment_score** |       0.008395       |
|        **stress_level**        |       -0.05178       |
|     **physical_activity**      |      -0.002868       |
|       **sleep_duration**       |       -0.07952       |
|    **healthy_eating_index**    |       0.01505        |
|         **gxe_effect**         |       -0.02233       |

Table continues below

|                                | physical_environment_score | stress_level |
|:------------------------------:|:--------------------------:|:------------:|
|            **age**             |          -0.02363          |   -0.06532   |
|         **weight_kg**          |          -0.04044          |   -0.08934   |
|         **height_cm**          |          -0.1022           |    -0.147    |
|            **neck**            |         -0.001796          |   -0.08747   |
|          **abdomen**           |          -0.02332          |   -0.04379   |
|           **thigh**            |          -0.04462          |   -0.06745   |
|            **knee**            |          -0.1356           |   -0.1395    |
|           **ankle**            |          -0.05918          |   -0.07937   |
|           **biceps**           |          0.01088           |   -0.06364   |
|          **forearm**           |          -0.01056          |   0.02742    |
|           **wrist**            |          -0.07075          |   -0.1021    |
|            **rs1**             |          0.02125           |   -0.1386    |
|            **rs2**             |          0.06366           |    0.1337    |
|            **rs3**             |          -0.04117          |   0.02588    |
|            **rs4**             |          -0.07579          |    0.1737    |
|            **rs5**             |          -0.05139          |   -0.08746   |
|            **rs6**             |          0.05591           |   0.02225    |
|        **num_children**        |           -0.13            |   -0.07869   |
|    **family_mealtime_freq**    |          0.008395          |   -0.05178   |
| **physical_environment_score** |             1              |   0.09289    |
|        **stress_level**        |          0.09289           |      1       |
|     **physical_activity**      |          0.07132           |   -0.02681   |
|       **sleep_duration**       |          -0.09233          |   0.02224    |
|    **healthy_eating_index**    |          -0.08028          |   0.06358    |
|         **gxe_effect**         |          -0.02522          |   -0.1467    |

Table continues below

|                                | physical_activity | sleep_duration |
|:------------------------------:|:-----------------:|:--------------:|
|            **age**             |      -0.1253      |    -0.04537    |
|         **weight_kg**          |      0.1012       |    0.05496     |
|         **height_cm**          |      0.06056      |   -0.009391    |
|            **neck**            |      0.06895      |     0.1458     |
|          **abdomen**           |      0.05438      |    0.07987     |
|           **thigh**            |      0.0898       |    0.01817     |
|            **knee**            |      0.03015      |    -0.01762    |
|           **ankle**            |      0.0824       |    -0.0371     |
|           **biceps**           |      0.1093       |    0.06898     |
|          **forearm**           |      0.1041       |     0.1191     |
|           **wrist**            |      0.08991      |   -0.009479    |
|            **rs1**             |      -0.1103      |    0.08817     |
|            **rs2**             |     -0.03517      |     -0.142     |
|            **rs3**             |      -0.1039      |    0.01076     |
|            **rs4**             |      0.0457       |     0.106      |
|            **rs5**             |      -0.1049      |    -0.1257     |
|            **rs6**             |      0.01894      |   -0.003312    |
|        **num_children**        |      0.02586      |    0.07246     |
|    **family_mealtime_freq**    |     -0.002868     |    -0.07952    |
| **physical_environment_score** |      0.07132      |    -0.09233    |
|        **stress_level**        |     -0.02681      |    0.02224     |
|     **physical_activity**      |         1         |     0.1047     |
|       **sleep_duration**       |      0.1047       |       1        |
|    **healthy_eating_index**    |      0.05909      |    -0.0411     |
|         **gxe_effect**         |     -0.07217      |     0.2956     |

Table continues below

|                                | healthy_eating_index | gxe_effect |
|:------------------------------:|:--------------------:|:----------:|
|            **age**             |       0.01366        |  -0.03343  |
|         **weight_kg**          |      -0.001971       |  0.03397   |
|         **height_cm**          |       0.009323       |   0.1271   |
|            **neck**            |      -0.008544       |   0.0835   |
|          **abdomen**           |       -0.04941       |  0.01293   |
|           **thigh**            |       -0.01847       |  0.03944   |
|            **knee**            |       0.01922        |  0.01225   |
|           **ankle**            |       0.01243        | 0.0009521  |
|           **biceps**           |       0.02064        |  0.04156   |
|          **forearm**           |       0.01588        | 0.0005114  |
|           **wrist**            |       0.08796        |  -0.03647  |
|            **rs1**             |       -0.08155       |   0.9521   |
|            **rs2**             |       0.04107        |  0.03138   |
|            **rs3**             |        0.1098        | -0.008795  |
|            **rs4**             |        0.0307        |   0.1739   |
|            **rs5**             |       -0.09141       |  0.04818   |
|            **rs6**             |        0.0138        |  0.02395   |
|        **num_children**        |       0.01976        |  -0.1018   |
|    **family_mealtime_freq**    |       0.01505        |  -0.02233  |
| **physical_environment_score** |       -0.08028       |  -0.02522  |
|        **stress_level**        |       0.06358        |  -0.1467   |
|     **physical_activity**      |       0.05909        |  -0.07217  |
|       **sleep_duration**       |       -0.0411        |   0.2956   |
|    **healthy_eating_index**    |          1           |  -0.07845  |
|         **gxe_effect**         |       -0.07845       |     1      |

``` r
cor_matrix <- cor(case1.bodyfat[, pred], use = "complete.obs")
upper_tri <- cor_matrix[upper.tri(cor_matrix)]
num_strong_correlations <- sum(upper_tri > 0.5)
corrplot(round(cor(case1.bodyfat[, pred]),2))
```

![](CaseStudy_files/figure-gfm/correlation-1.png)<!-- -->

An interesting feature of this data set is that many of the
anthropometric measures are intrinsically correlated. For example,
weight and abdomen circumference have a Pearson correlation coefficient
of 0.88. Further, a group of 34 pairs have correlation coefficients
greater than 0.5. As described above, these high correlations impose
some challenges in model development and interpretation.

## Regression

Based on our scientific expertise, let’s say we believe that abdominal
circumference and height are two central independent variables for
estimating body fat proportion. As such, we will not subject these to
variable selection (i.e., we will “force” them into all models). We
further believe that all other independent variables may be strongly
interrelated and exchangeable when used for body fat estimation.
Therefore, we will subject them to backward elimination with AIC as
stopping criterion.

### Full model

First estimate the global (i.e., full) model which includes all
candidate independent variables.

``` r
# Estimate full model ----------------------------------------------
formula <- paste("siri_simulated~", paste(pred, collapse = "+"))
full_mod <- lm(formula, data = case1.bodyfat, x = T, y = T)
full_est <- coef(full_mod)
full_se <- coef(summary(full_mod))[, "Std. Error"]
pander(summary(full_mod))
```

|                                |  Estimate  | Std. Error | t value  | Pr(\>\|t\|) |
|:------------------------------:|:----------:|:----------:|:--------:|:-----------:|
|        **(Intercept)**         |   -17.44   |   15.22    |  -1.146  |   0.2524    |
|            **age**             |  0.05218   |  0.02608   |  2.001   |   0.04601   |
|         **weight_kg**          |  -0.2042   |   0.0974   |  -2.096  |   0.0366    |
|         **height_cm**          |  -0.01277  |  0.05681   | -0.2248  |   0.8222    |
|            **neck**            |  -0.3987   |   0.1964   |  -2.03   |   0.04289   |
|          **abdomen**           |   0.9577   |  0.07291   |  13.14   |  7.553e-34  |
|           **thigh**            |   0.1858   |   0.1073   |  1.732   |   0.08392   |
|            **knee**            |  -0.1661   |   0.2006   | -0.8281  |    0.408    |
|           **ankle**            |   -0.124   |   0.2163   | -0.5731  |   0.5669    |
|           **biceps**           |   0.1464   |    0.14    |  1.046   |   0.2961    |
|          **forearm**           |   0.4145   |   0.1495   |  2.773   |  0.005776   |
|           **wrist**            |   -1.827   |   0.4413   |  -4.141  |  4.087e-05  |
|            **rs1**             |   3.874    |   1.556    |   2.49   |   0.01312   |
|            **rs2**             |   -1.274   |   0.3303   |  -3.857  |  0.0001308  |
|            **rs3**             |   0.7106   |   0.3953   |  1.798   |   0.07288   |
|            **rs4**             |  -0.3164   |   0.3285   | -0.9631  |    0.336    |
|            **rs5**             |   0.8264   |   0.355    |  2.328   |   0.02035   |
|            **rs6**             |   0.0807   |   0.3259   |  0.2477  |   0.8045    |
|        **num_children**        |   0.1501   |   0.1549   |  0.9693  |   0.3329    |
|    **family_mealtime_freq**    |  0.02933   |   0.0531   |  0.5524  |    0.581    |
| **physical_environment_score** | -0.0004278 |  0.02195   | -0.01949 |   0.9845    |
|        **stress_level**        |  -0.2942   |   0.1136   |  -2.589  |  0.009911   |
|     **physical_activity**      |   0.3706   |   0.4843   |  0.7652  |   0.4445    |
|       **sleep_duration**       |   -1.171   |   0.1414   |  -8.281  |  1.25e-15   |
|    **healthy_eating_index**    |  0.03348   |  0.02266   |  1.478   |   0.1402    |
|         **gxe_effect**         |   0.8576   |   0.1676   |  5.118   |  4.478e-07  |

| Observations | Residual Std. Error | $R^2$  | Adjusted $R^2$ |
|:------------:|:-------------------:|:------:|:--------------:|
|     502      |        4.717        | 0.8469 |     0.8389     |

Fitting linear model: formula

``` r
summary(full_mod)
```

    ## 
    ## Call:
    ## lm(formula = formula, data = case1.bodyfat, x = T, y = T)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -18.6516  -2.8468  -0.1146   2.7869  10.7738 
    ## 
    ## Coefficients:
    ##                              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                -1.744e+01  1.522e+01  -1.146 0.252428    
    ## age                         5.218e-02  2.608e-02   2.001 0.046010 *  
    ## weight_kg                  -2.042e-01  9.740e-02  -2.096 0.036605 *  
    ## height_cm                  -1.277e-02  5.681e-02  -0.225 0.822196    
    ## neck                       -3.987e-01  1.964e-01  -2.030 0.042890 *  
    ## abdomen                     9.577e-01  7.291e-02  13.135  < 2e-16 ***
    ## thigh                       1.858e-01  1.073e-01   1.732 0.083920 .  
    ## knee                       -1.661e-01  2.006e-01  -0.828 0.408029    
    ## ankle                      -1.240e-01  2.163e-01  -0.573 0.566853    
    ## biceps                      1.464e-01  1.400e-01   1.046 0.296145    
    ## forearm                     4.145e-01  1.495e-01   2.773 0.005776 ** 
    ## wrist                      -1.827e+00  4.413e-01  -4.141 4.09e-05 ***
    ## rs1                         3.874e+00  1.556e+00   2.490 0.013122 *  
    ## rs2                        -1.274e+00  3.303e-01  -3.857 0.000131 ***
    ## rs3                         7.106e-01  3.953e-01   1.798 0.072877 .  
    ## rs4                        -3.164e-01  3.285e-01  -0.963 0.335982    
    ## rs5                         8.264e-01  3.550e-01   2.328 0.020347 *  
    ## rs6                         8.070e-02  3.259e-01   0.248 0.804497    
    ## num_children                1.501e-01  1.549e-01   0.969 0.332913    
    ## family_mealtime_freq        2.933e-02  5.310e-02   0.552 0.580965    
    ## physical_environment_score -4.278e-04  2.195e-02  -0.019 0.984455    
    ## stress_level               -2.942e-01  1.136e-01  -2.589 0.009911 ** 
    ## physical_activity           3.706e-01  4.843e-01   0.765 0.444543    
    ## sleep_duration             -1.171e+00  1.414e-01  -8.281 1.25e-15 ***
    ## healthy_eating_index        3.348e-02  2.266e-02   1.478 0.140170    
    ## gxe_effect                  8.576e-01  1.676e-01   5.118 4.48e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 4.717 on 476 degrees of freedom
    ## Multiple R-squared:  0.8469, Adjusted R-squared:  0.8389 
    ## F-statistic: 105.3 on 25 and 476 DF,  p-value: < 2.2e-16

### Selected model

Next, use backwards elimination to identify the “selected” model.

``` r
# Selected model ---------------------------------------------------
sel_mod <- step(lm(formula, data = case1.bodyfat,  x=T,y=T), 
                direction = "backward",
                scope = list(upper = formula, 
                             # Force height and abdomen circumference into the model
                             lower = formula(siri_simulated~abdomen+height_cm)),
                trace = 0)
pander(summary(sel_mod))
```

|                          | Estimate | Std. Error | t value | Pr(\>\|t\|) |
|:------------------------:|:--------:|:----------:|:-------:|:-----------:|
|     **(Intercept)**      |  -15.01  |   14.15    |  -1.06  |   0.2895    |
|         **age**          | 0.04878  |  0.02479   |  1.967  |   0.04973   |
|      **weight_kg**       | -0.1879  |  0.08853   | -2.123  |   0.03427   |
|      **height_cm**       | -0.04205 |  0.05298   | -0.7936 |   0.4278    |
|         **neck**         | -0.3581  |   0.1869   | -1.916  |   0.05594   |
|       **abdomen**        |  0.949   |  0.06988   |  13.58  |  7.884e-36  |
|        **thigh**         |  0.1473  |  0.09353   |  1.575  |    0.116    |
|       **forearm**        |  0.4439  |   0.1445   |  3.071  |  0.002251   |
|        **wrist**         |  -1.92   |   0.4151   | -4.625  |  4.805e-06  |
|         **rs1**          |  4.094   |   1.486    |  2.754  |  0.006109   |
|         **rs2**          |  -1.261  |   0.3271   | -3.855  |  0.0001315  |
|         **rs3**          |  0.6013  |   0.3879   |  1.55   |   0.1218    |
|         **rs5**          |  0.8003  |   0.3469   |  2.307  |   0.02149   |
|     **stress_level**     | -0.3297  |   0.1086   | -3.034  |   0.00254   |
|    **sleep_duration**    |  -1.138  |   0.1387   | -8.208  |  2.044e-15  |
| **healthy_eating_index** | 0.03545  |   0.0224   |  1.582  |   0.1142    |
|      **gxe_effect**      |  0.824   |   0.1593   |  5.173  |  3.378e-07  |

| Observations | Residual Std. Error | $R^2$  | Adjusted $R^2$ |
|:------------:|:-------------------:|:------:|:--------------:|
|     502      |        4.698        | 0.8453 |     0.8402     |

Fitting linear model: siri_simulated ~ age + weight_kg + height_cm +
neck + abdomen + thigh + forearm + wrist + rs1 + rs2 + rs3 + rs5 +
stress_level + sleep_duration + healthy_eating_index + gxe_effect

``` r
sel_est <- coef(sel_mod)[c("(Intercept)", pred)]
sel_est[is.na(sel_est)] <- 0
names(sel_est) <- c("(Intercept)", pred)
sel_se <- coef(summary(sel_mod))[, "Std. Error"][c("(Intercept)", pred)]
sel_se[is.na(sel_se)] <- 0

summary(sel_mod)
```

    ## 
    ## Call:
    ## lm(formula = siri_simulated ~ age + weight_kg + height_cm + neck + 
    ##     abdomen + thigh + forearm + wrist + rs1 + rs2 + rs3 + rs5 + 
    ##     stress_level + sleep_duration + healthy_eating_index + gxe_effect, 
    ##     data = case1.bodyfat, x = T, y = T)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -18.859  -2.954  -0.070   2.825  11.145 
    ## 
    ## Coefficients:
    ##                       Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)          -15.00633   14.15225  -1.060 0.289513    
    ## age                    0.04878    0.02479   1.967 0.049726 *  
    ## weight_kg             -0.18793    0.08853  -2.123 0.034274 *  
    ## height_cm             -0.04205    0.05298  -0.794 0.427794    
    ## neck                  -0.35806    0.18687  -1.916 0.055942 .  
    ## abdomen                0.94901    0.06988  13.580  < 2e-16 ***
    ## thigh                  0.14729    0.09353   1.575 0.115953    
    ## forearm                0.44394    0.14454   3.071 0.002251 ** 
    ## wrist                 -1.91995    0.41510  -4.625 4.80e-06 ***
    ## rs1                    4.09367    1.48647   2.754 0.006109 ** 
    ## rs2                   -1.26072    0.32707  -3.855 0.000132 ***
    ## rs3                    0.60125    0.38794   1.550 0.121826    
    ## rs5                    0.80026    0.34694   2.307 0.021494 *  
    ## stress_level          -0.32966    0.10864  -3.034 0.002540 ** 
    ## sleep_duration        -1.13807    0.13865  -8.208 2.04e-15 ***
    ## healthy_eating_index   0.03545    0.02240   1.582 0.114213    
    ## gxe_effect             0.82397    0.15929   5.173 3.38e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 4.698 on 485 degrees of freedom
    ## Multiple R-squared:  0.8453, Adjusted R-squared:  0.8402 
    ## F-statistic: 165.6 on 16 and 485 DF,  p-value: < 2.2e-16

Observation: We identify a model where several variables were dropped
(reducing from 25 variables to 16). The adjusted R^2 only slightly
increases from from the global model to the selected model.

(Note that the coefficients and p-values for this model are identical to
those if we a priori selected the variables for glm)

### Bootstrap model

Finally, repeat backwards elimination using bootstrapping for a
stability investigation and post-selection inference.

``` r
# Bootstrap ----------------------------------------------------------
# Set number of bootstraps 
bootnum <- 1000
# Set up matrix for results 
boot_est <-  boot_se <- matrix(0, ncol = length(pred) + 1, nrow = bootnum,
                               dimnames = list(NULL, c("(Intercept)", pred)))
# Set seed for reproducible results 
set.seed(5437854)
# Perform bootstrapping (sampling with replacement)
for (i in 1:bootnum) {
  data_id <- sample(1:dim(case1.bodyfat)[1], replace = T)
  boot_mod <- step(lm(formula, data = case1.bodyfat[data_id, ], 
                             x=T, y=T), 
                         scope = list(upper = formula, 
                                      lower = 
                                        # Again forcing height and abdominal circumference
                                        # into all models 
                                        formula(siri_simulated ~ abdomen + height_cm)),
                         direction = "backward", trace = 0)
  boot_est[i, names(coef(boot_mod))] <- coef(boot_mod)
  boot_se[i, names(coef(boot_mod))] <- coef(summary(boot_mod))[, "Std. Error"]
}

# Compute final bootstrap estimates as the median 
boot_median <- apply(boot_est, 2, median)
boot_025per <- apply(boot_est, 2, function(x) quantile(x, 0.025))
boot_975per <- apply(boot_est, 2, function(x) quantile(x, 0.975))
boot.temp <- cbind(boot_median, boot_025per, boot_975per)
pander(boot.temp)
```

|                                | boot_median | boot_025per | boot_975per |
|:------------------------------:|:-----------:|:-----------:|:-----------:|
|        **(Intercept)**         |    -14.5    |   -46.37    |    20.22    |
|            **age**             |   0.05123   |      0      |   0.1044    |
|         **weight_kg**          |   -0.1929   |   -0.3852   |      0      |
|         **height_cm**          |  -0.02863   |   -0.1471   |   0.08253   |
|            **neck**            |   -0.4262   |   -0.8142   |      0      |
|          **abdomen**           |    0.953    |    0.779    |    1.102    |
|           **thigh**            |   0.1643    |      0      |   0.3446    |
|            **knee**            |      0      |   -0.5581   |      0      |
|           **ankle**            |      0      |   -0.5872   |   0.3582    |
|           **biceps**           |      0      |      0      |   0.4185    |
|          **forearm**           |   0.4022    |      0      |   0.6483    |
|           **wrist**            |   -1.834    |   -2.726    |   -0.9196   |
|            **rs1**             |    3.951    |      0      |    7.696    |
|            **rs2**             |   -1.312    |   -1.929    |   -0.6549   |
|            **rs3**             |   0.7165    |      0      |    1.435    |
|            **rs4**             |      0      |   -0.9303   |      0      |
|            **rs5**             |   0.8288    |      0      |    1.538    |
|            **rs6**             |      0      |   -0.5637   |   0.7257    |
|        **num_children**        |      0      |   -0.216    |   0.4981    |
|    **family_mealtime_freq**    |      0      |  -0.08284   |   0.1368    |
| **physical_environment_score** |      0      |  -0.04272   |   0.0427    |
|        **stress_level**        |    -0.3     |   -0.5542   |      0      |
|     **physical_activity**      |      0      |      0      |    1.288    |
|       **sleep_duration**       |   -1.179    |    -1.64    |   -0.7553   |
|    **healthy_eating_index**    |   0.0352    |      0      |   0.08102   |
|         **gxe_effect**         |   0.8443    |   0.4451    |    1.321    |

``` r
boot.temp
```

    ##                             boot_median  boot_025per boot_975per
    ## (Intercept)                -14.49748948 -46.36720000 20.22024775
    ## age                          0.05122974   0.00000000  0.10443557
    ## weight_kg                   -0.19286373  -0.38524101  0.00000000
    ## height_cm                   -0.02863082  -0.14714882  0.08253272
    ## neck                        -0.42619473  -0.81423700  0.00000000
    ## abdomen                      0.95304475   0.77897646  1.10211156
    ## thigh                        0.16428842   0.00000000  0.34458273
    ## knee                         0.00000000  -0.55812681  0.00000000
    ## ankle                        0.00000000  -0.58718315  0.35815279
    ## biceps                       0.00000000   0.00000000  0.41853519
    ## forearm                      0.40217735   0.00000000  0.64834995
    ## wrist                       -1.83441979  -2.72611623 -0.91962243
    ## rs1                          3.95078789   0.00000000  7.69587481
    ## rs2                         -1.31220687  -1.92881271 -0.65491007
    ## rs3                          0.71648638   0.00000000  1.43504567
    ## rs4                          0.00000000  -0.93029182  0.00000000
    ## rs5                          0.82877450   0.00000000  1.53765218
    ## rs6                          0.00000000  -0.56371002  0.72569382
    ## num_children                 0.00000000  -0.21595614  0.49808378
    ## family_mealtime_freq         0.00000000  -0.08284351  0.13676319
    ## physical_environment_score   0.00000000  -0.04272427  0.04270156
    ## stress_level                -0.29999599  -0.55420879  0.00000000
    ## physical_activity            0.00000000   0.00000000  1.28820797
    ## sleep_duration              -1.17856867  -1.63990320 -0.75528398
    ## healthy_eating_index         0.03520326   0.00000000  0.08102228
    ## gxe_effect                   0.84431282   0.44505203  1.32100354

Interpretation: Bootstrap estimates for the 2.5th (boot_025per) and
97.5th (boot_975per) percentiles can be interpreted as limits of 95%
confidence intervals obtained by resampling-based multi-model inference
estimated via bootstrap medians.

## Calculate bias estimators

### Bootstrap inclusion frequency

This number represents the percentage of time a variable was selected
across all bootstrapped realizations.

``` r
# Calculate bootstrap inclusion frequency 
boot_01 <- (boot_est != 0) * 1
boot_inclusion <- apply(boot_01, 2, function(x) sum(x) / length(x) * 100)
pander(boot_inclusion)
```

| (Intercept) | age | weight_kg | height_cm | neck | abdomen | thigh | knee |
|:-----------:|:---:|:---------:|:---------:|:----:|:-------:|:-----:|:----:|
|     100     | 76  |   75.8    |    100    | 77.8 |   100   | 62.6  | 29.3 |

Table continues below

| ankle | biceps | forearm | wrist | rs1  | rs2  | rs3  | rs4  | rs5  | rs6  |
|:-----:|:------:|:-------:|:-----:|:----:|:----:|:----:|:----:|:----:|:----:|
| 30.1  |  37.5  |  95.3   | 99.9  | 84.6 | 99.2 | 67.5 | 33.7 | 81.3 | 21.4 |

Table continues below

| num_children | family_mealtime_freq | physical_environment_score |
|:------------:|:--------------------:|:--------------------------:|
|     36.9     |         24.5         |            18.1            |

Table continues below

| stress_level | physical_activity | sleep_duration | healthy_eating_index |
|:------------:|:-----------------:|:--------------:|:--------------------:|
|     87.3     |       29.6        |      100       |         57.4         |

Table continues below

| gxe_effect |
|:----------:|
|    100     |

Interpretation: These numbers quantify how likely an independent
variable is to be selected (as %).

### Root mean squared difference (RMSD) ratio

The root mean squared difference (RMSD) ratio is computed as the
variable-specific RMSD of the bootstrap estimates divided by the
standard error in the global model (assumed unbiased), representing the
variance inflation/deflation consequent to variable selection.

``` r
# Calculate RMSD ratio 
sqe <- (t(boot_est) - full_est) ^ 2 # squared vector of the difference between boot and full 
rmsd <- apply(sqe, 1, function(x) sqrt(mean(x))) # find mean of that 
rmsdratio <- rmsd / full_se # evaluate relative to the SE (interpretation: if the selected model didn't mess anything up relative to full model, we end up with 1)
pander(rmsdratio)
```

| (Intercept) |  age  | weight_kg | height_cm | neck  | abdomen | thigh | knee  |
|:-----------:|:-----:|:---------:|:---------:|:-----:|:-------:|:-----:|:-----:|
|    1.224    | 1.235 |   1.266   |   1.124   | 1.249 |  1.176  | 1.207 | 1.035 |

Table continues below

| ankle | biceps | forearm | wrist |  rs1  |  rs2  |  rs3  |  rs4  |  rs5  |
|:-----:|:------:|:-------:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|
| 1.053 | 1.087  | 0.9516  | 1.018 | 1.364 | 1.033 | 1.221 | 1.055 | 1.278 |

Table continues below

|  rs6   | num_children | family_mealtime_freq | physical_environment_score |
|:------:|:------------:|:--------------------:|:--------------------------:|
| 0.8669 |    1.187     |        0.9654        |           0.7789           |

Table continues below

| stress_level | physical_activity | sleep_duration | healthy_eating_index |
|:------------:|:-----------------:|:--------------:|:--------------------:|
|    1.309     |      0.9899       |     1.556      |        1.261         |

Table continues below

| gxe_effect |
|:----------:|
|   1.393    |

Interpretation: Variable selection adds to uncertainty about the
regression coefficients, which is evidenced by RMSD ratios above 1
(except for ankle, wrist, rs4, rs6, family mealtime frequency, physical
environment score, and physical activity - which have values \<1 and
indicate variance deflation). So our process made estimates less certain
for all coefficients except for those with values \<1 (which were
included in fewer models … and when they were, they were included with a
set of variables that were similar to the full model).

RMSD ratios above 1: indicate increased uncertainty (variance inflation)
for those variables, meaning variable selection has made the estimates
less certain for coefficients with ratios above 1. RMSD ratios below 1:
indicate decreased uncertainty (variance deflation), suggesting that
variable selection made the estimates more stable (more certain) for
coefficients with ratios below 1. This often happens because those
variables were included less frequently or under similar conditions in
bootstrap models.

### Relative conditional bias

The relative conditional bias is calculated as the differences of the
mean of resampled regression coefficients and the global model
regression coefficient, divided by the global model regression
coefficient, representing the bias present if a variable was selected
because its regression coefficient appeared extreme in a specific
resample.

``` r
# Compute relative conditional bias (%)
boot_mean <- apply(boot_est, 2, mean)
boot_meanratio <- boot_mean / full_est
boot_relbias <- (boot_meanratio / (boot_inclusion / 100) - 1) * 100
pander(boot_relbias)
```

| (Intercept) |  age  | weight_kg | height_cm | neck  | abdomen | thigh | knee  |
|:-----------:|:-----:|:---------:|:---------:|:-----:|:-------:|:-----:|:-----:|
|   -24.79    | 20.37 |   15.23   |   157.3   | 25.16 | -1.257  | 19.66 | 109.5 |

Table continues below

| ankle | biceps | forearm | wrist  |  rs1  |  rs2  |  rs3  |  rs4  |  rs5  |
|:-----:|:------:|:-------:|:------:|:-----:|:-----:|:-----:|:-----:|:-----:|
| 145.1 | 95.16  | 0.7588  | 0.3199 | 15.14 | 2.917 | 27.47 | 103.3 | 14.61 |

Table continues below

|  rs6  | num_children | family_mealtime_freq | physical_environment_score |
|:-----:|:------------:|:--------------------:|:--------------------------:|
| 259.4 |    87.26     |        143.4         |           -504.1           |

Table continues below

| stress_level | physical_activity | sleep_duration | healthy_eating_index |
|:------------:|:-----------------:|:--------------:|:--------------------:|
|    13.24     |       132.2       |     0.9042     |        53.52         |

Table continues below

| gxe_effect |
|:----------:|
|   1.225    |

Interpretation: Relative conditional bias is negligible for abdomen,
wrist, rs2, sleep_duration, gxe_effect, forearm, rs1, stress_level, age,
neck, rs5, weight_kg, rs3, thigh, and healthy_eating_index (which we see
above all have bootstrap inclusion frequencies greater than 50%), but
becomes more relevant in variables for which selection is less sure. In
other words, these variables were selected less often, but when they
were they had larger effects.

# Results

## Overview

Next we will put together an overview table so we can see and discuss
the results side-by-side. This code replicates the style of that shown
in Table 5 of [Heinze, Wallisch, and Dunkler
(2018)](https://onlinelibrary.wiley.com/doi/full/10.1002/bimj.201700067).
This would also be the main results table that you could present in a
manuscript.

``` r
# Overview of estimates and measures --------------------------------
overview <- round(cbind(full_est, full_se, boot_inclusion, sel_est, sel_se, 
                        rmsdratio, boot_relbias, boot_median, boot_025per, 
                        boot_975per), 4)
overview <- overview[order(overview[,"boot_inclusion"], decreasing=T),] # Sort 

# Similar to Table 5 from Heinze, Wallisch, and Dunkler (2018) --------
pander(overview)
```

|                                | full_est | full_se | boot_inclusion | sel_est |
|:------------------------------:|:--------:|:-------:|:--------------:|:-------:|
|        **(Intercept)**         |  -17.44  |  15.22  |      100       | -15.01  |
|         **height_cm**          | -0.0128  | 0.0568  |      100       | -0.0421 |
|          **abdomen**           |  0.9577  | 0.0729  |      100       |  0.949  |
|       **sleep_duration**       |  -1.171  | 0.1414  |      100       | -1.138  |
|         **gxe_effect**         |  0.8576  | 0.1676  |      100       |  0.824  |
|           **wrist**            |  -1.827  | 0.4413  |      99.9      |  -1.92  |
|            **rs2**             |  -1.274  | 0.3303  |      99.2      | -1.261  |
|          **forearm**           |  0.4145  | 0.1495  |      95.3      | 0.4439  |
|        **stress_level**        | -0.2942  | 0.1136  |      87.3      | -0.3297 |
|            **rs1**             |  3.874   |  1.556  |      84.6      |  4.094  |
|            **rs5**             |  0.8264  |  0.355  |      81.3      | 0.8003  |
|            **neck**            | -0.3987  | 0.1964  |      77.8      | -0.3581 |
|            **age**             |  0.0522  | 0.0261  |       76       | 0.0488  |
|         **weight_kg**          | -0.2042  | 0.0974  |      75.8      | -0.1879 |
|            **rs3**             |  0.7106  | 0.3953  |      67.5      | 0.6013  |
|           **thigh**            |  0.1858  | 0.1073  |      62.6      | 0.1473  |
|    **healthy_eating_index**    |  0.0335  | 0.0227  |      57.4      | 0.0354  |
|           **biceps**           |  0.1464  |  0.14   |      37.5      |    0    |
|        **num_children**        |  0.1501  | 0.1549  |      36.9      |    0    |
|            **rs4**             | -0.3164  | 0.3285  |      33.7      |    0    |
|           **ankle**            |  -0.124  | 0.2163  |      30.1      |    0    |
|     **physical_activity**      |  0.3706  | 0.4843  |      29.6      |    0    |
|            **knee**            | -0.1661  | 0.2006  |      29.3      |    0    |
|    **family_mealtime_freq**    |  0.0293  | 0.0531  |      24.5      |    0    |
|            **rs6**             |  0.0807  | 0.3259  |      21.4      |    0    |
| **physical_environment_score** |  -4e-04  | 0.0219  |      18.1      |    0    |

Table continues below

|                                | sel_se | rmsdratio | boot_relbias |
|:------------------------------:|:------:|:---------:|:------------:|
|        **(Intercept)**         | 14.15  |   1.224   |    -24.79    |
|         **height_cm**          | 0.053  |   1.124   |    157.3     |
|          **abdomen**           | 0.0699 |   1.176   |    -1.257    |
|       **sleep_duration**       | 0.1387 |   1.556   |    0.9042    |
|         **gxe_effect**         | 0.1593 |   1.393   |    1.225     |
|           **wrist**            | 0.4151 |   1.018   |    0.3199    |
|            **rs2**             | 0.3271 |   1.032   |    2.917     |
|          **forearm**           | 0.1445 |  0.9516   |    0.7588    |
|        **stress_level**        | 0.1086 |   1.309   |    13.24     |
|            **rs1**             | 1.486  |   1.364   |    15.14     |
|            **rs5**             | 0.3469 |   1.278   |    14.61     |
|            **neck**            | 0.1869 |   1.249   |    25.16     |
|            **age**             | 0.0248 |   1.235   |    20.37     |
|         **weight_kg**          | 0.0885 |   1.266   |    15.23     |
|            **rs3**             | 0.3879 |   1.222   |    27.47     |
|           **thigh**            | 0.0935 |   1.207   |    19.66     |
|    **healthy_eating_index**    | 0.0224 |   1.261   |    53.52     |
|           **biceps**           |   0    |   1.087   |    95.16     |
|        **num_children**        |   0    |   1.187   |    87.26     |
|            **rs4**             |   0    |   1.055   |    103.3     |
|           **ankle**            |   0    |   1.053   |    145.1     |
|     **physical_activity**      |   0    |  0.9899   |    132.2     |
|            **knee**            |   0    |   1.035   |    109.5     |
|    **family_mealtime_freq**    |   0    |  0.9654   |    143.4     |
|            **rs6**             |   0    |  0.8669   |    259.4     |
| **physical_environment_score** |   0    |  0.7789   |    -504.1    |

Table continues below

|                                | boot_median | boot_025per | boot_975per |
|:------------------------------:|:-----------:|:-----------:|:-----------:|
|        **(Intercept)**         |    -14.5    |   -46.37    |    20.22    |
|         **height_cm**          |   -0.0286   |   -0.1471   |   0.0825    |
|          **abdomen**           |    0.953    |    0.779    |    1.102    |
|       **sleep_duration**       |   -1.179    |    -1.64    |   -0.7553   |
|         **gxe_effect**         |   0.8443    |   0.4451    |    1.321    |
|           **wrist**            |   -1.834    |   -2.726    |   -0.9196   |
|            **rs2**             |   -1.312    |   -1.929    |   -0.6549   |
|          **forearm**           |   0.4022    |      0      |   0.6483    |
|        **stress_level**        |    -0.3     |   -0.5542   |      0      |
|            **rs1**             |    3.951    |      0      |    7.696    |
|            **rs5**             |   0.8288    |      0      |    1.538    |
|            **neck**            |   -0.4262   |   -0.8142   |      0      |
|            **age**             |   0.0512    |      0      |   0.1044    |
|         **weight_kg**          |   -0.1929   |   -0.3852   |      0      |
|            **rs3**             |   0.7165    |      0      |    1.435    |
|           **thigh**            |   0.1643    |      0      |   0.3446    |
|    **healthy_eating_index**    |   0.0352    |      0      |    0.081    |
|           **biceps**           |      0      |      0      |   0.4185    |
|        **num_children**        |      0      |   -0.216    |   0.4981    |
|            **rs4**             |      0      |   -0.9303   |      0      |
|           **ankle**            |      0      |   -0.5872   |   0.3582    |
|     **physical_activity**      |      0      |      0      |    1.288    |
|            **knee**            |      0      |   -0.5581   |      0      |
|    **family_mealtime_freq**    |      0      |   -0.0828   |   0.1368    |
|            **rs6**             |      0      |   -0.5637   |   0.7257    |
| **physical_environment_score** |      0      |   -0.0427   |   0.0427    |

``` r
overview
```

    ##                            full_est full_se boot_inclusion  sel_est  sel_se
    ## (Intercept)                -17.4415 15.2213          100.0 -15.0063 14.1522
    ## height_cm                   -0.0128  0.0568          100.0  -0.0421  0.0530
    ## abdomen                      0.9577  0.0729          100.0   0.9490  0.0699
    ## sleep_duration              -1.1706  0.1414          100.0  -1.1381  0.1387
    ## gxe_effect                   0.8576  0.1676          100.0   0.8240  0.1593
    ## wrist                       -1.8274  0.4413           99.9  -1.9200  0.4151
    ## rs2                         -1.2739  0.3303           99.2  -1.2607  0.3271
    ## forearm                      0.4145  0.1495           95.3   0.4439  0.1445
    ## stress_level                -0.2942  0.1136           87.3  -0.3297  0.1086
    ## rs1                          3.8744  1.5561           84.6   4.0937  1.4865
    ## rs5                          0.8264  0.3550           81.3   0.8003  0.3469
    ## neck                        -0.3987  0.1964           77.8  -0.3581  0.1869
    ## age                          0.0522  0.0261           76.0   0.0488  0.0248
    ## weight_kg                   -0.2042  0.0974           75.8  -0.1879  0.0885
    ## rs3                          0.7106  0.3953           67.5   0.6013  0.3879
    ## thigh                        0.1858  0.1073           62.6   0.1473  0.0935
    ## healthy_eating_index         0.0335  0.0227           57.4   0.0354  0.0224
    ## biceps                       0.1464  0.1400           37.5   0.0000  0.0000
    ## num_children                 0.1501  0.1549           36.9   0.0000  0.0000
    ## rs4                         -0.3164  0.3285           33.7   0.0000  0.0000
    ## ankle                       -0.1240  0.2163           30.1   0.0000  0.0000
    ## physical_activity            0.3706  0.4843           29.6   0.0000  0.0000
    ## knee                        -0.1661  0.2006           29.3   0.0000  0.0000
    ## family_mealtime_freq         0.0293  0.0531           24.5   0.0000  0.0000
    ## rs6                          0.0807  0.3259           21.4   0.0000  0.0000
    ## physical_environment_score  -0.0004  0.0219           18.1   0.0000  0.0000
    ##                            rmsdratio boot_relbias boot_median boot_025per
    ## (Intercept)                   1.2239     -24.7901    -14.4975    -46.3672
    ## height_cm                     1.1243     157.3476     -0.0286     -0.1471
    ## abdomen                       1.1762      -1.2567      0.9530      0.7790
    ## sleep_duration                1.5563       0.9042     -1.1786     -1.6399
    ## gxe_effect                    1.3926       1.2249      0.8443      0.4451
    ## wrist                         1.0180       0.3199     -1.8344     -2.7261
    ## rs2                           1.0325       2.9169     -1.3122     -1.9288
    ## forearm                       0.9516       0.7588      0.4022      0.0000
    ## stress_level                  1.3090      13.2390     -0.3000     -0.5542
    ## rs1                           1.3636      15.1413      3.9508      0.0000
    ## rs5                           1.2783      14.6072      0.8288      0.0000
    ## neck                          1.2493      25.1567     -0.4262     -0.8142
    ## age                           1.2351      20.3693      0.0512      0.0000
    ## weight_kg                     1.2659      15.2317     -0.1929     -0.3852
    ## rs3                           1.2215      27.4729      0.7165      0.0000
    ## thigh                         1.2068      19.6571      0.1643      0.0000
    ## healthy_eating_index          1.2606      53.5160      0.0352      0.0000
    ## biceps                        1.0869      95.1585      0.0000      0.0000
    ## num_children                  1.1870      87.2621      0.0000     -0.2160
    ## rs4                           1.0552     103.2616      0.0000     -0.9303
    ## ankle                         1.0528     145.0509      0.0000     -0.5872
    ## physical_activity             0.9899     132.2329      0.0000      0.0000
    ## knee                          1.0349     109.4968      0.0000     -0.5581
    ## family_mealtime_freq          0.9654     143.3856      0.0000     -0.0828
    ## rs6                           0.8669     259.4143      0.0000     -0.5637
    ## physical_environment_score    0.7789    -504.1432      0.0000     -0.0427
    ##                            boot_975per
    ## (Intercept)                    20.2202
    ## height_cm                       0.0825
    ## abdomen                         1.1021
    ## sleep_duration                 -0.7553
    ## gxe_effect                      1.3210
    ## wrist                          -0.9196
    ## rs2                            -0.6549
    ## forearm                         0.6483
    ## stress_level                    0.0000
    ## rs1                             7.6959
    ## rs5                             1.5377
    ## neck                            0.0000
    ## age                             0.1044
    ## weight_kg                       0.0000
    ## rs3                             1.4350
    ## thigh                           0.3446
    ## healthy_eating_index            0.0810
    ## biceps                          0.4185
    ## num_children                    0.4981
    ## rs4                             0.0000
    ## ankle                           0.3582
    ## physical_activity               1.2882
    ## knee                            0.0000
    ## family_mealtime_freq            0.1368
    ## rs6                             0.7257
    ## physical_environment_score      0.0427

Looking at everything together in a sorted/organized table, we note that
the bootstrapped results resemble the global model —- but interestingly
we see that some of the effect estimates have opposite signs in the
selected model!!!!!! This further supports the importance of adding
stability investigations to your models when using variable
selection!!!! (I didn’t even plan that through data simulation!!!!!).

## Model selection frequency

Data related to the number of times a group of variables were selected
together in backwards elimination can also be obtained. The following
code will replicate the style of Table 6 of [Heinze, Wallisch, and
Dunkler
(2018)](https://onlinelibrary.wiley.com/doi/full/10.1002/bimj.201700067).

``` r
# Model frequency ---------------------------------------------------
# Similar to Table 6 from Heinze, Wallisch, and Dunkler (2018) ------
pred_ord <- names(boot_inclusion)[order(boot_inclusion, decreasing = T)]
boot_01 <- boot_01[, pred_ord]
boot_01 <- cbind(boot_01, count = rep(1, times = bootnum))
boot_modfreq <- aggregate(count ~ ., data = boot_01, sum)
boot_modfreq[, "percent"] <- boot_modfreq$count / bootnum * 100
boot_modfreq <- boot_modfreq[order(boot_modfreq[, "percent"], decreasing = T), ]
boot_modfreq[, "cum_percent"] <- cumsum(boot_modfreq$percent)
boot_modfreq <- boot_modfreq[boot_modfreq[, "cum_percent"] <= 80, ]
if (dim(boot_modfreq)[1] > 20) boot_modfreq <- boot_modfreq[1:20, ]


pander(cbind("Predictors"= apply(boot_modfreq[,c(2:14)], 1, 
                          function(x) paste(names(x[x==1]), collapse=" ")),
      boot_modfreq[,c("count", "percent", "cum_percent")]))
```

|         |                                              Predictors                                               | count | percent | cum_percent |
|:-------:|:-----------------------------------------------------------------------------------------------------:|:-----:|:-------:|:-----------:|
| **20**  |   height_cm abdomen sleep_duration gxe_effect wrist rs2 forearm stress_level rs1 rs5 neck weight_kg   |   4   |   0.4   |     0.4     |
| **13**  | height_cm abdomen sleep_duration gxe_effect wrist rs2 forearm stress_level rs1 rs5 neck age weight_kg |   3   |   0.3   |     0.7     |
| **49**  | height_cm abdomen sleep_duration gxe_effect wrist rs2 forearm stress_level rs1 rs5 neck age weight_kg |   3   |   0.3   |      1      |
| **64**  |   height_cm abdomen sleep_duration gxe_effect wrist rs2 forearm stress_level rs1 rs5 age weight_kg    |   3   |   0.3   |     1.3     |
| **66**  | height_cm abdomen sleep_duration gxe_effect wrist rs2 forearm stress_level rs1 rs5 neck age weight_kg |   3   |   0.3   |     1.6     |
| **93**  | height_cm abdomen sleep_duration gxe_effect wrist rs2 forearm stress_level rs1 rs5 neck age weight_kg |   3   |   0.3   |     1.9     |
| **153** |   height_cm abdomen sleep_duration gxe_effect wrist rs2 forearm stress_level rs1 neck age weight_kg   |   3   |   0.3   |     2.2     |
| **342** |   height_cm abdomen sleep_duration gxe_effect wrist rs2 forearm stress_level rs1 rs5 age weight_kg    |   3   |   0.3   |     2.5     |
| **465** | height_cm abdomen sleep_duration gxe_effect wrist rs2 forearm stress_level rs1 rs5 neck age weight_kg |   3   |   0.3   |     2.8     |
| **15**  |     height_cm abdomen sleep_duration gxe_effect wrist rs2 forearm stress_level rs1 rs5 weight_kg      |   2   |   0.2   |      3      |
| **16**  |   height_cm abdomen sleep_duration gxe_effect wrist rs2 forearm stress_level rs1 rs5 neck weight_kg   |   2   |   0.2   |     3.2     |
| **25**  |   height_cm abdomen sleep_duration gxe_effect wrist rs2 forearm stress_level rs1 rs5 age weight_kg    |   2   |   0.2   |     3.4     |
| **30**  |   height_cm abdomen sleep_duration gxe_effect wrist rs2 forearm stress_level rs1 rs5 age weight_kg    |   2   |   0.2   |     3.6     |
| **31**  |   height_cm abdomen sleep_duration gxe_effect wrist rs2 forearm stress_level rs1 neck age weight_kg   |   2   |   0.2   |     3.8     |
| **33**  | height_cm abdomen sleep_duration gxe_effect wrist rs2 forearm stress_level rs1 rs5 neck age weight_kg |   2   |   0.2   |      4      |
| **44**  |     height_cm abdomen sleep_duration gxe_effect wrist rs2 forearm stress_level rs1 rs5 weight_kg      |   2   |   0.2   |     4.2     |
| **46**  | height_cm abdomen sleep_duration gxe_effect wrist rs2 forearm stress_level rs1 rs5 neck age weight_kg |   2   |   0.2   |     4.4     |
| **55**  | height_cm abdomen sleep_duration gxe_effect wrist rs2 forearm stress_level rs1 rs5 neck age weight_kg |   2   |   0.2   |     4.6     |
| **60**  |   height_cm abdomen sleep_duration gxe_effect wrist rs2 forearm stress_level rs5 neck age weight_kg   |   2   |   0.2   |     4.8     |
| **85**  | height_cm abdomen sleep_duration gxe_effect wrist rs2 forearm stress_level rs1 rs5 neck age weight_kg |   2   |   0.2   |      5      |

``` r
cbind("Predictors"= apply(boot_modfreq[,c(2:14)], 1, 
                          function(x) paste(names(x[x==1]), collapse=" ")),
      boot_modfreq[,c("count", "percent", "cum_percent")])
```

    ##                                                                                                Predictors
    ## 20      height_cm abdomen sleep_duration gxe_effect wrist rs2 forearm stress_level rs1 rs5 neck weight_kg
    ## 13  height_cm abdomen sleep_duration gxe_effect wrist rs2 forearm stress_level rs1 rs5 neck age weight_kg
    ## 49  height_cm abdomen sleep_duration gxe_effect wrist rs2 forearm stress_level rs1 rs5 neck age weight_kg
    ## 64       height_cm abdomen sleep_duration gxe_effect wrist rs2 forearm stress_level rs1 rs5 age weight_kg
    ## 66  height_cm abdomen sleep_duration gxe_effect wrist rs2 forearm stress_level rs1 rs5 neck age weight_kg
    ## 93  height_cm abdomen sleep_duration gxe_effect wrist rs2 forearm stress_level rs1 rs5 neck age weight_kg
    ## 153     height_cm abdomen sleep_duration gxe_effect wrist rs2 forearm stress_level rs1 neck age weight_kg
    ## 342      height_cm abdomen sleep_duration gxe_effect wrist rs2 forearm stress_level rs1 rs5 age weight_kg
    ## 465 height_cm abdomen sleep_duration gxe_effect wrist rs2 forearm stress_level rs1 rs5 neck age weight_kg
    ## 15           height_cm abdomen sleep_duration gxe_effect wrist rs2 forearm stress_level rs1 rs5 weight_kg
    ## 16      height_cm abdomen sleep_duration gxe_effect wrist rs2 forearm stress_level rs1 rs5 neck weight_kg
    ## 25       height_cm abdomen sleep_duration gxe_effect wrist rs2 forearm stress_level rs1 rs5 age weight_kg
    ## 30       height_cm abdomen sleep_duration gxe_effect wrist rs2 forearm stress_level rs1 rs5 age weight_kg
    ## 31      height_cm abdomen sleep_duration gxe_effect wrist rs2 forearm stress_level rs1 neck age weight_kg
    ## 33  height_cm abdomen sleep_duration gxe_effect wrist rs2 forearm stress_level rs1 rs5 neck age weight_kg
    ## 44           height_cm abdomen sleep_duration gxe_effect wrist rs2 forearm stress_level rs1 rs5 weight_kg
    ## 46  height_cm abdomen sleep_duration gxe_effect wrist rs2 forearm stress_level rs1 rs5 neck age weight_kg
    ## 55  height_cm abdomen sleep_duration gxe_effect wrist rs2 forearm stress_level rs1 rs5 neck age weight_kg
    ## 60      height_cm abdomen sleep_duration gxe_effect wrist rs2 forearm stress_level rs5 neck age weight_kg
    ## 85  height_cm abdomen sleep_duration gxe_effect wrist rs2 forearm stress_level rs1 rs5 neck age weight_kg
    ##     count percent cum_percent
    ## 20      4     0.4         0.4
    ## 13      3     0.3         0.7
    ## 49      3     0.3         1.0
    ## 64      3     0.3         1.3
    ## 66      3     0.3         1.6
    ## 93      3     0.3         1.9
    ## 153     3     0.3         2.2
    ## 342     3     0.3         2.5
    ## 465     3     0.3         2.8
    ## 15      2     0.2         3.0
    ## 16      2     0.2         3.2
    ## 25      2     0.2         3.4
    ## 30      2     0.2         3.6
    ## 31      2     0.2         3.8
    ## 33      2     0.2         4.0
    ## 44      2     0.2         4.2
    ## 46      2     0.2         4.4
    ## 55      2     0.2         4.6
    ## 60      2     0.2         4.8
    ## 85      2     0.2         5.0

The above table shows the combinations of predictors and the number of
times/percentage of time those predictors were selected together. The
highest selection frequency is only 2%.

``` r
# Model frequency in % of selected model ----------------------------
sel_modfreq <- sum(apply(boot_01[, -dim(boot_01)[2]], 1, function(x)
    identical(((sel_est != 0) * 1)[pred_ord], x))) / bootnum * 100
sel_modfreq
```

    ## [1] 0.2

Note also that our selected model (via Backwards elimination) was never
exactly reproduced in the bootstrapping approach.

## Pairwise inclusion frequency

Pairwise inclusion frequencies can also be calculated with this
approach. The following code will produce a table similar to Table S2 of
[Heinze, Wallisch, and Dunkler
(2018)](https://onlinelibrary.wiley.com/doi/full/10.1002/bimj.201700067).

``` r
# Pairwise inclusion frequency in % ----------------------------------
# Similar to Supplementary Table 2 in supporting information from Heinze, Wallisch, and Dunkler (2018) --------------------
pval <- 0.01
boot_pairfreq <- matrix(100, ncol = length(pred), nrow = length(pred),
                        dimnames = list(pred_ord[-1], 
                                        pred_ord[-1]))

expect_pairfreq <- NULL
combis <- combn(pred_ord[-1], 2)

for (i in 1:dim(combis)[2]) {
  boot_pairfreq[combis[1, i], combis[2, i]] <-
    sum(apply(boot_01[, combis[, i]], 1, sum) == 2) / bootnum * 100
  expect_pairfreq[i] <-
    boot_inclusion[grepl(combis[1, i], pred_ord)][1] *
    boot_inclusion[grepl(combis[2, i], pred_ord)][1] / 100
  boot_pairfreq[combis[2, i], combis[1, i]] <-
    ifelse(is(suppressWarnings(try(chisq.test(boot_01[, combis[1, i]], 
                                              boot_01[, combis[2, i]]), 
                                   silent = T)), "try-error"),
      NA, ifelse(suppressWarnings(
                   chisq.test(boot_01[, combis[1, i]],
                              boot_01[, combis[2, i]])$p.value) > pval,
        "", ifelse(as.numeric(boot_pairfreq[combis[1, i], combis[2, i]]) < 
                     expect_pairfreq[i], "-", "+")))
}
diag(boot_pairfreq) <- boot_inclusion[pred_ord][-1]
pander(boot_pairfreq, quote = F)
```

|                                | height_cm | abdomen | sleep_duration |
|:------------------------------:|:---------:|:-------:|:--------------:|
|         **height_cm**          |    100    |   100   |      100       |
|          **abdomen**           |    NA     |   100   |      100       |
|       **sleep_duration**       |    NA     |   NA    |      100       |
|         **gxe_effect**         |    NA     |   NA    |       NA       |
|           **wrist**            |    NA     |   NA    |       NA       |
|            **rs2**             |    NA     |   NA    |       NA       |
|          **forearm**           |    NA     |   NA    |       NA       |
|        **stress_level**        |    NA     |   NA    |       NA       |
|            **rs1**             |    NA     |   NA    |       NA       |
|            **rs5**             |    NA     |   NA    |       NA       |
|            **neck**            |    NA     |   NA    |       NA       |
|            **age**             |    NA     |   NA    |       NA       |
|         **weight_kg**          |    NA     |   NA    |       NA       |
|            **rs3**             |    NA     |   NA    |       NA       |
|           **thigh**            |    NA     |   NA    |       NA       |
|    **healthy_eating_index**    |    NA     |   NA    |       NA       |
|           **biceps**           |    NA     |   NA    |       NA       |
|        **num_children**        |    NA     |   NA    |       NA       |
|            **rs4**             |    NA     |   NA    |       NA       |
|           **ankle**            |    NA     |   NA    |       NA       |
|     **physical_activity**      |    NA     |   NA    |       NA       |
|            **knee**            |    NA     |   NA    |       NA       |
|    **family_mealtime_freq**    |    NA     |   NA    |       NA       |
|            **rs6**             |    NA     |   NA    |       NA       |
| **physical_environment_score** |    NA     |   NA    |       NA       |

Table continues below

|                                | gxe_effect | wrist | rs2  | forearm |
|:------------------------------:|:----------:|:-----:|:----:|:-------:|
|         **height_cm**          |    100     | 99.9  | 99.2 |  95.3   |
|          **abdomen**           |    100     | 99.9  | 99.2 |  95.3   |
|       **sleep_duration**       |    100     | 99.9  | 99.2 |  95.3   |
|         **gxe_effect**         |    100     | 99.9  | 99.2 |  95.3   |
|           **wrist**            |     NA     | 99.9  | 99.1 |  95.2   |
|            **rs2**             |     NA     |       | 99.2 |  94.5   |
|          **forearm**           |     NA     |       |      |  95.3   |
|        **stress_level**        |     NA     |       |      |         |
|            **rs1**             |     NA     |       |      |         |
|            **rs5**             |     NA     |       |      |         |
|            **neck**            |     NA     |       |      |         |
|            **age**             |     NA     |       |      |         |
|         **weight_kg**          |     NA     |       |      |         |
|            **rs3**             |     NA     |       |      |         |
|           **thigh**            |     NA     |       |      |         |
|    **healthy_eating_index**    |     NA     |       |      |         |
|           **biceps**           |     NA     |       |      |   \+    |
|        **num_children**        |     NA     |       |      |         |
|            **rs4**             |     NA     |       |      |         |
|           **ankle**            |     NA     |       |      |         |
|     **physical_activity**      |     NA     |       |      |   \+    |
|            **knee**            |     NA     |       |      |         |
|    **family_mealtime_freq**    |     NA     |       |      |         |
|            **rs6**             |     NA     |       |      |         |
| **physical_environment_score** |     NA     |       |      |         |

Table continues below

|                                | stress_level | rs1  | rs5  | neck | age  |
|:------------------------------:|:------------:|:----:|:----:|:----:|:----:|
|         **height_cm**          |     87.3     | 84.6 | 81.3 | 77.8 |  76  |
|          **abdomen**           |     87.3     | 84.6 | 81.3 | 77.8 |  76  |
|       **sleep_duration**       |     87.3     | 84.6 | 81.3 | 77.8 |  76  |
|         **gxe_effect**         |     87.3     | 84.6 | 81.3 | 77.8 |  76  |
|           **wrist**            |     87.2     | 84.5 | 81.3 | 77.7 |  76  |
|            **rs2**             |     86.6     | 83.9 | 80.7 |  77  | 75.2 |
|          **forearm**           |     83.2     | 80.5 | 77.6 | 74.7 |  73  |
|        **stress_level**        |     87.3     | 76.2 | 71.9 | 69.4 | 64.7 |
|            **rs1**             |      \+      | 84.6 | 68.7 | 65.1 | 62.8 |
|            **rs5**             |              |      | 81.3 | 62.3 | 61.5 |
|            **neck**            |      \+      |      |      | 77.8 | 59.8 |
|            **age**             |      \+      |  \+  |      |      |  76  |
|         **weight_kg**          |              |      |      |  \-  |  \-  |
|            **rs3**             |              |      |      |  \-  |      |
|           **thigh**            |      \+      |  \+  |      |  \+  |  \+  |
|    **healthy_eating_index**    |              |      |      |      |      |
|           **biceps**           |              |      |  \+  |      |  \+  |
|        **num_children**        |      \+      |      |      |  \-  |  \-  |
|            **rs4**             |      \+      |      |      |      |      |
|           **ankle**            |              |      |  \+  |      |      |
|     **physical_activity**      |              |  \-  |  \-  |  \-  |      |
|            **knee**            |              |      |      |      |  \+  |
|    **family_mealtime_freq**    |              |  \-  |      |      |      |
|            **rs6**             |              |      |      |      |      |
| **physical_environment_score** |              |      |      |      |      |

Table continues below

|                                | weight_kg | rs3  | thigh |
|:------------------------------:|:---------:|:----:|:-----:|
|         **height_cm**          |   75.8    | 67.5 | 62.6  |
|          **abdomen**           |   75.8    | 67.5 | 62.6  |
|       **sleep_duration**       |   75.8    | 67.5 | 62.6  |
|         **gxe_effect**         |   75.8    | 67.5 | 62.6  |
|           **wrist**            |   75.7    | 67.4 | 62.5  |
|            **rs2**             |   75.2    | 66.9 | 62.1  |
|          **forearm**           |   72.9    | 63.7 | 60.1  |
|        **stress_level**        |    66     | 58.1 | 52.5  |
|            **rs1**             |   63.7    | 55.7 | 50.6  |
|            **rs5**             |   61.3    | 54.1 | 50.7  |
|            **neck**            |   54.5    | 54.4 | 45.2  |
|            **age**             |   53.7    | 51.9 | 54.2  |
|         **weight_kg**          |   75.8    | 50.6 | 50.6  |
|            **rs3**             |           | 67.5 | 43.8  |
|           **thigh**            |    \+     |      | 62.6  |
|    **healthy_eating_index**    |           |  \-  |       |
|           **biceps**           |    \+     |      |  \+   |
|        **num_children**        |           |      |  \+   |
|            **rs4**             |    \+     |      |  \+   |
|           **ankle**            |    \+     |      |       |
|     **physical_activity**      |           |      |       |
|            **knee**            |    \-     |      |  \+   |
|    **family_mealtime_freq**    |           |      |  \-   |
|            **rs6**             |           |      |       |
| **physical_environment_score** |           |      |       |

Table continues below

|                                | healthy_eating_index | biceps | num_children |
|:------------------------------:|:--------------------:|:------:|:------------:|
|         **height_cm**          |         57.4         |  37.5  |     36.9     |
|          **abdomen**           |         57.4         |  37.5  |     36.9     |
|       **sleep_duration**       |         57.4         |  37.5  |     36.9     |
|         **gxe_effect**         |         57.4         |  37.5  |     36.9     |
|           **wrist**            |         57.3         |  37.5  |     36.9     |
|            **rs2**             |         56.7         |  37.2  |     36.8     |
|          **forearm**           |         55.4         |  34.7  |     35.3     |
|        **stress_level**        |         50.1         |  31.9  |     29.8     |
|            **rs1**             |         48.6         |  32.7  |     31.1     |
|            **rs5**             |         47.1         |  28.2  |     28.9     |
|            **neck**            |         43.1         |  30.4  |     26.5     |
|            **age**             |         42.9         |   26   |     29.9     |
|         **weight_kg**          |         44.6         |  32.5  |     28.2     |
|            **rs3**             |         34.9         |  25.5  |     26.6     |
|           **thigh**            |         35.4         |  17.6  |     27.7     |
|    **healthy_eating_index**    |         57.4         |  19.5  |     20.9     |
|           **biceps**           |          \+          |  37.5  |     12.9     |
|        **num_children**        |                      |        |     36.9     |
|            **rs4**             |                      |   \+   |              |
|           **ankle**            |                      |        |              |
|     **physical_activity**      |                      |        |              |
|            **knee**            |                      |        |      \+      |
|    **family_mealtime_freq**    |                      |   \-   |              |
|            **rs6**             |                      |   \-   |              |
| **physical_environment_score** |                      |        |              |

Table continues below

|                                | rs4  | ankle | physical_activity | knee |
|:------------------------------:|:----:|:-----:|:-----------------:|:----:|
|         **height_cm**          | 33.7 | 30.1  |       29.6        | 29.3 |
|          **abdomen**           | 33.7 | 30.1  |       29.6        | 29.3 |
|       **sleep_duration**       | 33.7 | 30.1  |       29.6        | 29.3 |
|         **gxe_effect**         | 33.7 | 30.1  |       29.6        | 29.3 |
|           **wrist**            | 33.7 |  30   |       29.6        | 29.2 |
|            **rs2**             | 33.4 | 29.7  |       29.3        | 29.2 |
|          **forearm**           | 31.8 | 29.2  |       27.3        | 28.7 |
|        **stress_level**        | 25.5 | 25.9  |       25.8        | 25.4 |
|            **rs1**             | 27.5 | 26.3  |       23.5        |  25  |
|            **rs5**             | 28.3 | 22.7  |       25.6        | 23.4 |
|            **neck**            | 25.1 | 24.3  |       20.9        | 22.9 |
|            **age**             | 25.2 | 22.2  |       22.4        | 25.1 |
|         **weight_kg**          | 28.2 | 18.5  |       23.6        | 19.9 |
|            **rs3**             | 22.2 | 20.4  |       20.5        |  21  |
|           **thigh**            | 23.8 | 18.3  |        19         | 24.1 |
|    **healthy_eating_index**    | 20.1 | 17.7  |       17.2        | 16.5 |
|           **biceps**           | 15.7 |  9.5  |       11.2        | 10.2 |
|        **num_children**        | 13.2 | 12.2  |        9.7        | 13.2 |
|            **rs4**             | 33.7 | 10.4  |       12.4        | 9.2  |
|           **ankle**            |      | 30.1  |        8.8        | 7.1  |
|     **physical_activity**      |  \-  |       |       29.6        | 7.3  |
|            **knee**            |      |       |                   | 29.3 |
|    **family_mealtime_freq**    |      |       |                   |      |
|            **rs6**             |      |       |                   |      |
| **physical_environment_score** |      |       |                   |  \-  |

Table continues below

|                                | family_mealtime_freq | rs6  |
|:------------------------------:|:--------------------:|:----:|
|         **height_cm**          |         24.5         | 21.4 |
|          **abdomen**           |         24.5         | 21.4 |
|       **sleep_duration**       |         24.5         | 21.4 |
|         **gxe_effect**         |         24.5         | 21.4 |
|           **wrist**            |         24.5         | 21.4 |
|            **rs2**             |         24.3         | 21.1 |
|          **forearm**           |         23.9         | 20.5 |
|        **stress_level**        |         21.5         | 18.9 |
|            **rs1**             |         18.9         | 18.9 |
|            **rs5**             |         21.3         | 17.5 |
|            **neck**            |         20.5         | 16.5 |
|            **age**             |         19.9         | 15.4 |
|         **weight_kg**          |         17.7         | 17.4 |
|            **rs3**             |         17.7         | 14.6 |
|           **thigh**            |         17.7         | 12.6 |
|    **healthy_eating_index**    |         13.3         |  12  |
|           **biceps**           |         7.1          | 9.8  |
|        **num_children**        |         9.4          | 8.1  |
|            **rs4**             |         8.5          | 6.7  |
|           **ankle**            |         6.1          | 6.4  |
|     **physical_activity**      |         7.1          | 5.2  |
|            **knee**            |          8           |  6   |
|    **family_mealtime_freq**    |         24.5         | 6.4  |
|            **rs6**             |                      | 21.4 |
| **physical_environment_score** |                      |      |

Table continues below

|                                | physical_environment_score |
|:------------------------------:|:--------------------------:|
|         **height_cm**          |            18.1            |
|          **abdomen**           |            18.1            |
|       **sleep_duration**       |            18.1            |
|         **gxe_effect**         |            18.1            |
|           **wrist**            |            18.1            |
|            **rs2**             |            17.7            |
|          **forearm**           |            17.1            |
|        **stress_level**        |             16             |
|            **rs1**             |            14.8            |
|            **rs5**             |            15.2            |
|            **neck**            |            13.4            |
|            **age**             |            14.7            |
|         **weight_kg**          |            13.1            |
|            **rs3**             |            12.4            |
|           **thigh**            |            10.9            |
|    **healthy_eating_index**    |            10.1            |
|           **biceps**           |            6.8             |
|        **num_children**        |            6.3             |
|            **rs4**             |            5.8             |
|           **ankle**            |             5              |
|     **physical_activity**      |            5.5             |
|            **knee**            |             7              |
|    **family_mealtime_freq**    |            4.7             |
|            **rs6**             |            3.4             |
| **physical_environment_score** |            18.1            |

``` r
print(boot_pairfreq, quote = F)
```

    ##                            height_cm abdomen sleep_duration gxe_effect wrist
    ## height_cm                  100       100     100            100        99.9 
    ## abdomen                    <NA>      100     100            100        99.9 
    ## sleep_duration             <NA>      <NA>    100            100        99.9 
    ## gxe_effect                 <NA>      <NA>    <NA>           100        99.9 
    ## wrist                      <NA>      <NA>    <NA>           <NA>       99.9 
    ## rs2                        <NA>      <NA>    <NA>           <NA>            
    ## forearm                    <NA>      <NA>    <NA>           <NA>            
    ## stress_level               <NA>      <NA>    <NA>           <NA>            
    ## rs1                        <NA>      <NA>    <NA>           <NA>            
    ## rs5                        <NA>      <NA>    <NA>           <NA>            
    ## neck                       <NA>      <NA>    <NA>           <NA>            
    ## age                        <NA>      <NA>    <NA>           <NA>            
    ## weight_kg                  <NA>      <NA>    <NA>           <NA>            
    ## rs3                        <NA>      <NA>    <NA>           <NA>            
    ## thigh                      <NA>      <NA>    <NA>           <NA>            
    ## healthy_eating_index       <NA>      <NA>    <NA>           <NA>            
    ## biceps                     <NA>      <NA>    <NA>           <NA>            
    ## num_children               <NA>      <NA>    <NA>           <NA>            
    ## rs4                        <NA>      <NA>    <NA>           <NA>            
    ## ankle                      <NA>      <NA>    <NA>           <NA>            
    ## physical_activity          <NA>      <NA>    <NA>           <NA>            
    ## knee                       <NA>      <NA>    <NA>           <NA>            
    ## family_mealtime_freq       <NA>      <NA>    <NA>           <NA>            
    ## rs6                        <NA>      <NA>    <NA>           <NA>            
    ## physical_environment_score <NA>      <NA>    <NA>           <NA>            
    ##                            rs2  forearm stress_level rs1  rs5  neck age 
    ## height_cm                  99.2 95.3    87.3         84.6 81.3 77.8 76  
    ## abdomen                    99.2 95.3    87.3         84.6 81.3 77.8 76  
    ## sleep_duration             99.2 95.3    87.3         84.6 81.3 77.8 76  
    ## gxe_effect                 99.2 95.3    87.3         84.6 81.3 77.8 76  
    ## wrist                      99.1 95.2    87.2         84.5 81.3 77.7 76  
    ## rs2                        99.2 94.5    86.6         83.9 80.7 77   75.2
    ## forearm                         95.3    83.2         80.5 77.6 74.7 73  
    ## stress_level                            87.3         76.2 71.9 69.4 64.7
    ## rs1                                     +            84.6 68.7 65.1 62.8
    ## rs5                                                       81.3 62.3 61.5
    ## neck                                    +                      77.8 59.8
    ## age                                     +            +              76  
    ## weight_kg                                                      -    -   
    ## rs3                                                            -        
    ## thigh                                   +            +         +    +   
    ## healthy_eating_index                                                    
    ## biceps                          +                         +         +   
    ## num_children                            +                      -    -   
    ## rs4                                     +                               
    ## ankle                                                     +             
    ## physical_activity               +                    -    -    -        
    ## knee                                                                +   
    ## family_mealtime_freq                                 -                  
    ## rs6                                                                     
    ## physical_environment_score                                              
    ##                            weight_kg rs3  thigh healthy_eating_index biceps
    ## height_cm                  75.8      67.5 62.6  57.4                 37.5  
    ## abdomen                    75.8      67.5 62.6  57.4                 37.5  
    ## sleep_duration             75.8      67.5 62.6  57.4                 37.5  
    ## gxe_effect                 75.8      67.5 62.6  57.4                 37.5  
    ## wrist                      75.7      67.4 62.5  57.3                 37.5  
    ## rs2                        75.2      66.9 62.1  56.7                 37.2  
    ## forearm                    72.9      63.7 60.1  55.4                 34.7  
    ## stress_level               66        58.1 52.5  50.1                 31.9  
    ## rs1                        63.7      55.7 50.6  48.6                 32.7  
    ## rs5                        61.3      54.1 50.7  47.1                 28.2  
    ## neck                       54.5      54.4 45.2  43.1                 30.4  
    ## age                        53.7      51.9 54.2  42.9                 26    
    ## weight_kg                  75.8      50.6 50.6  44.6                 32.5  
    ## rs3                                  67.5 43.8  34.9                 25.5  
    ## thigh                      +              62.6  35.4                 17.6  
    ## healthy_eating_index                 -          57.4                 19.5  
    ## biceps                     +              +     +                    37.5  
    ## num_children                              +                                
    ## rs4                        +              +                          +     
    ## ankle                      +                                               
    ## physical_activity                                                          
    ## knee                       -              +                                
    ## family_mealtime_freq                      -                          -     
    ## rs6                                                                  -     
    ## physical_environment_score                                                 
    ##                            num_children rs4  ankle physical_activity knee
    ## height_cm                  36.9         33.7 30.1  29.6              29.3
    ## abdomen                    36.9         33.7 30.1  29.6              29.3
    ## sleep_duration             36.9         33.7 30.1  29.6              29.3
    ## gxe_effect                 36.9         33.7 30.1  29.6              29.3
    ## wrist                      36.9         33.7 30    29.6              29.2
    ## rs2                        36.8         33.4 29.7  29.3              29.2
    ## forearm                    35.3         31.8 29.2  27.3              28.7
    ## stress_level               29.8         25.5 25.9  25.8              25.4
    ## rs1                        31.1         27.5 26.3  23.5              25  
    ## rs5                        28.9         28.3 22.7  25.6              23.4
    ## neck                       26.5         25.1 24.3  20.9              22.9
    ## age                        29.9         25.2 22.2  22.4              25.1
    ## weight_kg                  28.2         28.2 18.5  23.6              19.9
    ## rs3                        26.6         22.2 20.4  20.5              21  
    ## thigh                      27.7         23.8 18.3  19                24.1
    ## healthy_eating_index       20.9         20.1 17.7  17.2              16.5
    ## biceps                     12.9         15.7 9.5   11.2              10.2
    ## num_children               36.9         13.2 12.2  9.7               13.2
    ## rs4                                     33.7 10.4  12.4              9.2 
    ## ankle                                        30.1  8.8               7.1 
    ## physical_activity                       -          29.6              7.3 
    ## knee                       +                                         29.3
    ## family_mealtime_freq                                                     
    ## rs6                                                                      
    ## physical_environment_score                                           -   
    ##                            family_mealtime_freq rs6  physical_environment_score
    ## height_cm                  24.5                 21.4 18.1                      
    ## abdomen                    24.5                 21.4 18.1                      
    ## sleep_duration             24.5                 21.4 18.1                      
    ## gxe_effect                 24.5                 21.4 18.1                      
    ## wrist                      24.5                 21.4 18.1                      
    ## rs2                        24.3                 21.1 17.7                      
    ## forearm                    23.9                 20.5 17.1                      
    ## stress_level               21.5                 18.9 16                        
    ## rs1                        18.9                 18.9 14.8                      
    ## rs5                        21.3                 17.5 15.2                      
    ## neck                       20.5                 16.5 13.4                      
    ## age                        19.9                 15.4 14.7                      
    ## weight_kg                  17.7                 17.4 13.1                      
    ## rs3                        17.7                 14.6 12.4                      
    ## thigh                      17.7                 12.6 10.9                      
    ## healthy_eating_index       13.3                 12   10.1                      
    ## biceps                     7.1                  9.8  6.8                       
    ## num_children               9.4                  8.1  6.3                       
    ## rs4                        8.5                  6.7  5.8                       
    ## ankle                      6.1                  6.4  5                         
    ## physical_activity          7.1                  5.2  5.5                       
    ## knee                       8                    6    7                         
    ## family_mealtime_freq       24.5                 6.4  4.7                       
    ## rs6                                             21.4 3.4                       
    ## physical_environment_score                           18.1

LACEY START HERE

Pairwise inclusion frequencies inform about “rope teams” and
“competitors” among the independent variables. For example, weight and
neck circumference were both selected in only 54.5% of the resamples,
while one would expect a frequency of 59.0% ( = 75.8% × 77.8%) given
independent selection. Therefore, the pair is flagged with “-” in the
lower triangle of this table.

Alternatively, thigh and neck are flagged with “+” because they are
simultaneously selected in 48.7% of the resamples, while the expectation
under independence is only 45.2%.

Interestingly, age forms a “rope team” with thigh, biceps, and knee, but
age is a competitor to weight and number of children.

In this table, a significance of a chi-squared test at the 0.01 level is
the formal criterion for the flags.

## Shrinkage factors

See Section 4.5 (p. 75) of Frank E. Harrell, Jr. Regression Modeling
Strategies:
<https://antivirus.uclv.edu.cu/update/libros/Mathematics%20and%20Statistics/Regression%20Modeling%20Strategies%20-%20Frank%20E.%20Harrell%20%2C%20Jr.%2C%202nd%20ed.%202015%20-%20978-3-319-19425-7.pdf>

Shrinkage factors are a tool to improve the stability and
generalizability of regression models, particularly in situations where
there are concerns about overfitting or multicollinearity. They help
strike a balance between bias and variance in model estimation,
ultimately leading to more reliable and interpretable results.

Post-estimation shrinkage factors (global and parameterwise) can be
calculated for this approach. The shrinkage factors are obtained by
leave‐one‐out resampling.

Global Shrinkage Factor: This modifies all regression coefficients by
the same factor. It’s like a uniform adjustment applied to all
coefficients in the model. The idea is to shrink coefficients towards
zero or some other value to reduce their variance, thereby improving the
model’s generalization ability.

Parameterwise Shrinkage Factor: This adjusts each regression coefficient
individually. Regression coefficients for which selection is less stable
are shrunken more strongly (further away from 1) than coefficients for
which selection is stable (closer to 1).

NOTE: The authors of this paper DO NOT recommend parameterwise shrinkage
factors for data sets with high correlations between the variables.
Nevertheless, the following code computes these data similar to that
shown in Table S3 of [Heinze, Wallisch, and Dunkler
(2018)](https://onlinelibrary.wiley.com/doi/full/10.1002/bimj.201700067).

``` r
# Shrinkage factors --------------------------------------------------
# Similar to Supplementary Table 3 in supporting information from Heinze, Wallisch, and Dunkler (2018) -------------
sel_mod_shrinkg <- shrink(sel_mod, type="global")
sel_mod_shrinkp <- shrink(sel_mod, type="parameterwise")

pander(sel_mod_shrinkp$ShrinkageFactors)
```

|  age   | weight_kg | height_cm |  neck  | abdomen | thigh  | forearm | wrist  |
|:------:|:---------:|:---------:|:------:|:-------:|:------:|:-------:|:------:|
| 0.6546 |   1.297   |  0.1723   | 0.3685 |  1.037  | 0.8249 | 0.7874  | 0.9557 |

Table continues below

|  rs1  |  rs2   |  rs3   |  rs5   | stress_level | sleep_duration |
|:-----:|:------:|:------:|:------:|:------------:|:--------------:|
| 1.138 | 0.9504 | 0.7255 | 0.7039 |    0.9186    |     0.9851     |

Table continues below

| healthy_eating_index | gxe_effect |
|:--------------------:|:----------:|
|        0.7163        |   0.9147   |

``` r
# Global Shrinkage Factor
sel_mod_shrinkg$ShrinkageFactors # Very close to 1, so can be neglected
```

    ## [1] 0.9860103

``` r
# Parameterise shrinkage
sel_mod_shrinkp_vcov<-vcov(sel_mod_shrinkp)
```

    ##                      (Intercept)          age    weight_kg   height_cm
    ## (Intercept)           62.1646000  0.556150000 -1.894310000  1.60834000
    ## age                    0.5561500  0.262251000 -0.063194100 -0.02540940
    ## weight_kg             -1.8943100 -0.063194100  0.116084000 -0.03477290
    ## height_cm              1.6083400 -0.025409400 -0.034772900  0.21928400
    ## neck                   0.6090950  0.006219820 -0.013755500 -0.00331756
    ## abdomen               -0.1847900 -0.016668900  0.015077200 -0.00486033
    ## thigh                 -1.7260500  0.077310900  0.080121100 -0.04510570
    ## forearm               -0.5466410  0.019235500  0.023690700  0.00375648
    ## wrist                  0.8560780  0.047303300 -0.029484500 -0.00952228
    ## rs1                   -0.3726960 -0.013788200 -0.004123790  0.00134136
    ## rs2                    0.0516987 -0.004610280 -0.005151580 -0.00438268
    ## rs3                    0.1493690 -0.015775200 -0.013022400  0.00323501
    ## rs5                   -0.1602260  0.004681840 -0.002885720  0.01495380
    ## stress_level           0.2557460 -0.014666600  0.009420880  0.00826791
    ## sleep_duration         0.1787890  0.000547481  0.000960486  0.00133905
    ## healthy_eating_index  -0.6584210 -0.009355370  0.010256700 -0.00890055
    ## gxe_effect             0.1512000  0.006050840  0.003822690 -0.00276127
    ##                             neck      abdomen        thigh     forearm
    ## (Intercept)           0.60909500 -0.184790000 -1.726050000 -0.54664100
    ## age                   0.00621982 -0.016668900  0.077310900  0.01923550
    ## weight_kg            -0.01375550  0.015077200  0.080121100  0.02369070
    ## height_cm            -0.00331756 -0.004860330 -0.045105700  0.00375648
    ## neck                  0.15138600  0.003997920 -0.009912500  0.03431510
    ## abdomen               0.00399792  0.003642770 -0.002279840  0.00277608
    ## thigh                -0.00991250 -0.002279840  0.312984000 -0.00795969
    ## forearm               0.03431510  0.002776080 -0.007959690  0.11990500
    ## wrist                -0.01358840 -0.003735480  0.003516160  0.00902730
    ## rs1                  -0.01541110 -0.000585117 -0.018320100 -0.00931379
    ## rs2                  -0.00809331 -0.001462920 -0.000359204 -0.01687410
    ## rs3                   0.01144830 -0.001891040  0.016137100  0.00457923
    ## rs5                   0.01794160  0.001014640  0.000531560  0.01393900
    ## stress_level          0.01423970  0.002530820 -0.004145720  0.01470250
    ## sleep_duration       -0.00335334  0.000223223  0.000930981  0.00302982
    ## healthy_eating_index -0.00289833  0.004173060 -0.000744044  0.00371087
    ## gxe_effect            0.00773809  0.000569885  0.009190630  0.00528891
    ##                             wrist          rs1          rs2         rs3
    ## (Intercept)           8.56078e-01 -0.372696000  0.051698700  0.14936900
    ## age                   4.73033e-02 -0.013788200 -0.004610280 -0.01577520
    ## weight_kg            -2.94845e-02 -0.004123790 -0.005151580 -0.01302240
    ## height_cm            -9.52228e-03  0.001341360 -0.004382680  0.00323501
    ## neck                 -1.35884e-02 -0.015411100 -0.008093310  0.01144830
    ## abdomen              -3.73548e-03 -0.000585117 -0.001462920 -0.00189104
    ## thigh                 3.51616e-03 -0.018320100 -0.000359204  0.01613710
    ## forearm               9.02730e-03 -0.009313790 -0.016874100  0.00457923
    ## wrist                 4.21142e-02 -0.004785400 -0.003811140  0.00223431
    ## rs1                  -4.78540e-03  0.122974000  0.008144740 -0.02044980
    ## rs2                  -3.81114e-03  0.008144740  0.070886900  0.02298840
    ## rs3                   2.23431e-03 -0.020449800  0.022988400  0.43996800
    ## rs5                  -7.33157e-05 -0.009496880 -0.005291240 -0.01035170
    ## stress_level          2.17044e-03  0.004914510 -0.015298500 -0.00138155
    ## sleep_duration        2.94876e-03 -0.028176500  0.002213840  0.00757581
    ## healthy_eating_index  1.26799e-02  0.015203400 -0.000670400 -0.05823500
    ## gxe_effect            1.76133e-03 -0.063557800 -0.002994470  0.01100710
    ##                               rs5 stress_level sleep_duration
    ## (Intercept)          -1.60226e-01   0.25574600    0.178789000
    ## age                   4.68184e-03  -0.01466660    0.000547481
    ## weight_kg            -2.88572e-03   0.00942088    0.000960486
    ## height_cm             1.49538e-02   0.00826791    0.001339050
    ## neck                  1.79416e-02   0.01423970   -0.003353340
    ## abdomen               1.01464e-03   0.00253082    0.000223223
    ## thigh                 5.31560e-04  -0.00414572    0.000930981
    ## forearm               1.39390e-02   0.01470250    0.003029820
    ## wrist                -7.33157e-05   0.00217044    0.002948760
    ## rs1                  -9.49688e-03   0.00491451   -0.028176500
    ## rs2                  -5.29124e-03  -0.01529850    0.002213840
    ## rs3                  -1.03517e-02  -0.00138155    0.007575810
    ## rs5                   2.02033e-01  -0.00808484   -0.004136090
    ## stress_level         -8.08484e-03   0.11529900   -0.003850510
    ## sleep_duration       -4.13609e-03  -0.00385051    0.016051600
    ## healthy_eating_index  2.67561e-02   0.01525840   -0.003920890
    ## gxe_effect            3.36597e-03  -0.00501830    0.015958200
    ##                      healthy_eating_index   gxe_effect
    ## (Intercept)                  -0.658421000  0.151200000
    ## age                          -0.009355370  0.006050840
    ## weight_kg                     0.010256700  0.003822690
    ## height_cm                    -0.008900550 -0.002761270
    ## neck                         -0.002898330  0.007738090
    ## abdomen                       0.004173060  0.000569885
    ## thigh                        -0.000744044  0.009190630
    ## forearm                       0.003710870  0.005288910
    ## wrist                         0.012679900  0.001761330
    ## rs1                           0.015203400 -0.063557800
    ## rs2                          -0.000670400 -0.002994470
    ## rs3                          -0.058235000  0.011007100
    ## rs5                           0.026756100  0.003365970
    ## stress_level                  0.015258400 -0.005018300
    ## sleep_duration               -0.003920890  0.015958200
    ## healthy_eating_index          0.419332000 -0.006807270
    ## gxe_effect                   -0.006807270  0.034941000

``` r
pander(round(cbind("Shrinkage factors" = sel_mod_shrinkp$ShrinkageFactors,
            "SE of shrinkage factors" = sqrt(diag(sel_mod_shrinkp_vcov))[-1],
            "Correlation matrix of shrinkage factors" = 
              cov2cor(sel_mod_shrinkp_vcov)[-1,-1]),4))
```

|                          | Shrinkage factors | SE of shrinkage factors |
|:------------------------:|:-----------------:|:-----------------------:|
|         **age**          |      0.6546       |         0.5121          |
|      **weight_kg**       |       1.297       |         0.3407          |
|      **height_cm**       |      0.1723       |         0.4683          |
|         **neck**         |      0.3685       |         0.3891          |
|       **abdomen**        |       1.037       |         0.0604          |
|        **thigh**         |      0.8249       |         0.5594          |
|       **forearm**        |      0.7874       |         0.3463          |
|        **wrist**         |      0.9557       |         0.2052          |
|         **rs1**          |       1.138       |         0.3507          |
|         **rs2**          |      0.9504       |         0.2662          |
|         **rs3**          |      0.7255       |         0.6633          |
|         **rs5**          |      0.7039       |         0.4495          |
|     **stress_level**     |      0.9186       |         0.3396          |
|    **sleep_duration**    |      0.9851       |         0.1267          |
| **healthy_eating_index** |      0.7163       |         0.6476          |
|      **gxe_effect**      |      0.9147       |         0.1869          |

Table continues below

|                          |   age   | weight_kg | height_cm |  neck   | abdomen |
|:------------------------:|:-------:|:---------:|:---------:|:-------:|:-------:|
|         **age**          |    1    |  -0.3622  |  -0.106   | 0.0312  | -0.5393 |
|      **weight_kg**       | -0.3622 |     1     |  -0.2179  | -0.1038 | 0.7332  |
|      **height_cm**       | -0.106  |  -0.2179  |     1     | -0.0182 | -0.172  |
|         **neck**         | 0.0312  |  -0.1038  |  -0.0182  |    1    | 0.1702  |
|       **abdomen**        | -0.5393 |  0.7332   |  -0.172   | 0.1702  |    1    |
|        **thigh**         | 0.2698  |  0.4203   |  -0.1722  | -0.0455 | -0.0675 |
|       **forearm**        | 0.1085  |  0.2008   |  0.0232   | 0.2547  | 0.1328  |
|        **wrist**         | 0.4501  |  -0.4217  |  -0.0991  | -0.1702 | -0.3016 |
|         **rs1**          | -0.0768 |  -0.0345  |  0.0082   | -0.1129 | -0.0276 |
|         **rs2**          | -0.0338 |  -0.0568  |  -0.0352  | -0.0781 | -0.091  |
|         **rs3**          | -0.0464 |  -0.0576  |  0.0104   | 0.0444  | -0.0472 |
|         **rs5**          | 0.0203  |  -0.0188  |   0.071   | 0.1026  | 0.0374  |
|     **stress_level**     | -0.0843 |  0.0814   |   0.052   | 0.1078  | 0.1235  |
|    **sleep_duration**    | 0.0084  |  0.0223   |  0.0226   | -0.068  | 0.0292  |
| **healthy_eating_index** | -0.0282 |  0.0465   |  -0.0294  | -0.0115 | 0.1068  |
|      **gxe_effect**      | 0.0632  |   0.06    |  -0.0315  | 0.1064  | 0.0505  |

Table continues below

|                          |  thigh  | forearm |  wrist  |   rs1   |   rs2   |
|:------------------------:|:-------:|:-------:|:-------:|:-------:|:-------:|
|         **age**          | 0.2698  | 0.1085  | 0.4501  | -0.0768 | -0.0338 |
|      **weight_kg**       | 0.4203  | 0.2008  | -0.4217 | -0.0345 | -0.0568 |
|      **height_cm**       | -0.1722 | 0.0232  | -0.0991 | 0.0082  | -0.0352 |
|         **neck**         | -0.0455 | 0.2547  | -0.1702 | -0.1129 | -0.0781 |
|       **abdomen**        | -0.0675 | 0.1328  | -0.3016 | -0.0276 | -0.091  |
|        **thigh**         |    1    | -0.0411 | 0.0306  | -0.0934 | -0.0024 |
|       **forearm**        | -0.0411 |    1    |  0.127  | -0.0767 | -0.183  |
|        **wrist**         | 0.0306  |  0.127  |    1    | -0.0665 | -0.0698 |
|         **rs1**          | -0.0934 | -0.0767 | -0.0665 |    1    | 0.0872  |
|         **rs2**          | -0.0024 | -0.183  | -0.0698 | 0.0872  |    1    |
|         **rs3**          | 0.0435  | 0.0199  | 0.0164  | -0.0879 | 0.1302  |
|         **rs5**          | 0.0021  | 0.0896  | -8e-04  | -0.0603 | -0.0442 |
|     **stress_level**     | -0.0218 |  0.125  | 0.0311  | 0.0413  | -0.1692 |
|    **sleep_duration**    | 0.0131  | 0.0691  | 0.1134  | -0.6342 | 0.0656  |
| **healthy_eating_index** | -0.0021 | 0.0165  | 0.0954  |  0.067  | -0.0039 |
|      **gxe_effect**      | 0.0879  | 0.0817  | 0.0459  | -0.9696 | -0.0602 |

Table continues below

|                          |   rs3   |   rs5   | stress_level | sleep_duration |
|:------------------------:|:-------:|:-------:|:------------:|:--------------:|
|         **age**          | -0.0464 | 0.0203  |   -0.0843    |     0.0084     |
|      **weight_kg**       | -0.0576 | -0.0188 |    0.0814    |     0.0223     |
|      **height_cm**       | 0.0104  |  0.071  |    0.052     |     0.0226     |
|         **neck**         | 0.0444  | 0.1026  |    0.1078    |     -0.068     |
|       **abdomen**        | -0.0472 | 0.0374  |    0.1235    |     0.0292     |
|        **thigh**         | 0.0435  | 0.0021  |   -0.0218    |     0.0131     |
|       **forearm**        | 0.0199  | 0.0896  |    0.125     |     0.0691     |
|        **wrist**         | 0.0164  | -8e-04  |    0.0311    |     0.1134     |
|         **rs1**          | -0.0879 | -0.0603 |    0.0413    |    -0.6342     |
|         **rs2**          | 0.1302  | -0.0442 |   -0.1692    |     0.0656     |
|         **rs3**          |    1    | -0.0347 |   -0.0061    |     0.0901     |
|         **rs5**          | -0.0347 |    1    |    -0.053    |    -0.0726     |
|     **stress_level**     | -0.0061 | -0.053  |      1       |    -0.0895     |
|    **sleep_duration**    | 0.0901  | -0.0726 |   -0.0895    |       1        |
| **healthy_eating_index** | -0.1356 | 0.0919  |    0.0694    |    -0.0478     |
|      **gxe_effect**      | 0.0888  | 0.0401  |   -0.0791    |     0.6738     |

Table continues below

|                          | healthy_eating_index | gxe_effect |
|:------------------------:|:--------------------:|:----------:|
|         **age**          |       -0.0282        |   0.0632   |
|      **weight_kg**       |        0.0465        |    0.06    |
|      **height_cm**       |       -0.0294        |  -0.0315   |
|         **neck**         |       -0.0115        |   0.1064   |
|       **abdomen**        |        0.1068        |   0.0505   |
|        **thigh**         |       -0.0021        |   0.0879   |
|       **forearm**        |        0.0165        |   0.0817   |
|        **wrist**         |        0.0954        |   0.0459   |
|         **rs1**          |        0.067         |  -0.9696   |
|         **rs2**          |       -0.0039        |  -0.0602   |
|         **rs3**          |       -0.1356        |   0.0888   |
|         **rs5**          |        0.0919        |   0.0401   |
|     **stress_level**     |        0.0694        |  -0.0791   |
|    **sleep_duration**    |       -0.0478        |   0.6738   |
| **healthy_eating_index** |          1           |  -0.0562   |
|      **gxe_effect**      |       -0.0562        |     1      |

``` r
round(cbind("Shrinkage factors" = sel_mod_shrinkp$ShrinkageFactors,
            "SE of shrinkage factors" = sqrt(diag(sel_mod_shrinkp_vcov))[-1],
            "Correlation matrix of shrinkage factors" = 
              cov2cor(sel_mod_shrinkp_vcov)[-1,-1]),4)
```

    ##                      Shrinkage factors SE of shrinkage factors     age
    ## age                             0.6546                  0.5121  1.0000
    ## weight_kg                       1.2973                  0.3407 -0.3622
    ## height_cm                       0.1723                  0.4683 -0.1060
    ## neck                            0.3685                  0.3891  0.0312
    ## abdomen                         1.0366                  0.0604 -0.5393
    ## thigh                           0.8249                  0.5594  0.2698
    ## forearm                         0.7874                  0.3463  0.1085
    ## wrist                           0.9557                  0.2052  0.4501
    ## rs1                             1.1384                  0.3507 -0.0768
    ## rs2                             0.9504                  0.2662 -0.0338
    ## rs3                             0.7255                  0.6633 -0.0464
    ## rs5                             0.7039                  0.4495  0.0203
    ## stress_level                    0.9186                  0.3396 -0.0843
    ## sleep_duration                  0.9851                  0.1267  0.0084
    ## healthy_eating_index            0.7163                  0.6476 -0.0282
    ## gxe_effect                      0.9147                  0.1869  0.0632
    ##                      weight_kg height_cm    neck abdomen   thigh forearm
    ## age                    -0.3622   -0.1060  0.0312 -0.5393  0.2698  0.1085
    ## weight_kg               1.0000   -0.2179 -0.1038  0.7332  0.4203  0.2008
    ## height_cm              -0.2179    1.0000 -0.0182 -0.1720 -0.1722  0.0232
    ## neck                   -0.1038   -0.0182  1.0000  0.1702 -0.0455  0.2547
    ## abdomen                 0.7332   -0.1720  0.1702  1.0000 -0.0675  0.1328
    ## thigh                   0.4203   -0.1722 -0.0455 -0.0675  1.0000 -0.0411
    ## forearm                 0.2008    0.0232  0.2547  0.1328 -0.0411  1.0000
    ## wrist                  -0.4217   -0.0991 -0.1702 -0.3016  0.0306  0.1270
    ## rs1                    -0.0345    0.0082 -0.1129 -0.0276 -0.0934 -0.0767
    ## rs2                    -0.0568   -0.0352 -0.0781 -0.0910 -0.0024 -0.1830
    ## rs3                    -0.0576    0.0104  0.0444 -0.0472  0.0435  0.0199
    ## rs5                    -0.0188    0.0710  0.1026  0.0374  0.0021  0.0896
    ## stress_level            0.0814    0.0520  0.1078  0.1235 -0.0218  0.1250
    ## sleep_duration          0.0223    0.0226 -0.0680  0.0292  0.0131  0.0691
    ## healthy_eating_index    0.0465   -0.0294 -0.0115  0.1068 -0.0021  0.0165
    ## gxe_effect              0.0600   -0.0315  0.1064  0.0505  0.0879  0.0817
    ##                        wrist     rs1     rs2     rs3     rs5 stress_level
    ## age                   0.4501 -0.0768 -0.0338 -0.0464  0.0203      -0.0843
    ## weight_kg            -0.4217 -0.0345 -0.0568 -0.0576 -0.0188       0.0814
    ## height_cm            -0.0991  0.0082 -0.0352  0.0104  0.0710       0.0520
    ## neck                 -0.1702 -0.1129 -0.0781  0.0444  0.1026       0.1078
    ## abdomen              -0.3016 -0.0276 -0.0910 -0.0472  0.0374       0.1235
    ## thigh                 0.0306 -0.0934 -0.0024  0.0435  0.0021      -0.0218
    ## forearm               0.1270 -0.0767 -0.1830  0.0199  0.0896       0.1250
    ## wrist                 1.0000 -0.0665 -0.0698  0.0164 -0.0008       0.0311
    ## rs1                  -0.0665  1.0000  0.0872 -0.0879 -0.0603       0.0413
    ## rs2                  -0.0698  0.0872  1.0000  0.1302 -0.0442      -0.1692
    ## rs3                   0.0164 -0.0879  0.1302  1.0000 -0.0347      -0.0061
    ## rs5                  -0.0008 -0.0603 -0.0442 -0.0347  1.0000      -0.0530
    ## stress_level          0.0311  0.0413 -0.1692 -0.0061 -0.0530       1.0000
    ## sleep_duration        0.1134 -0.6342  0.0656  0.0901 -0.0726      -0.0895
    ## healthy_eating_index  0.0954  0.0670 -0.0039 -0.1356  0.0919       0.0694
    ## gxe_effect            0.0459 -0.9696 -0.0602  0.0888  0.0401      -0.0791
    ##                      sleep_duration healthy_eating_index gxe_effect
    ## age                          0.0084              -0.0282     0.0632
    ## weight_kg                    0.0223               0.0465     0.0600
    ## height_cm                    0.0226              -0.0294    -0.0315
    ## neck                        -0.0680              -0.0115     0.1064
    ## abdomen                      0.0292               0.1068     0.0505
    ## thigh                        0.0131              -0.0021     0.0879
    ## forearm                      0.0691               0.0165     0.0817
    ## wrist                        0.1134               0.0954     0.0459
    ## rs1                         -0.6342               0.0670    -0.9696
    ## rs2                          0.0656              -0.0039    -0.0602
    ## rs3                          0.0901              -0.1356     0.0888
    ## rs5                         -0.0726               0.0919     0.0401
    ## stress_level                -0.0895               0.0694    -0.0791
    ## sleep_duration               1.0000              -0.0478     0.6738
    ## healthy_eating_index        -0.0478               1.0000    -0.0562
    ## gxe_effect                   0.6738              -0.0562     1.0000

Interpretation: In this example, model shrinkage is 0.9860103 which is
very close to 1 and hence can be neglected.

These quantities were obtained using the R package
[`shrink`](https://cran.r-project.org/web/packages/shrink/index.html)
[(Dunkler, Sauerbrei and Heinze,
2016)](https://www.jstatsoft.org/article/view/v069i08).

# References

Courvoisier, D. S., Combescure, C., Agoritsas, T., Gayet-Ageron, A., &
Perneger, T. V. (2011). Performance of logistic regression modeling:
Beyond the number of events per variable, the role of data structure.
Journal of Clinical Epidemiology, 64, 993–1000.

Dunkler, D., Sauerbrei, W., and Heinze, G. (2016). Global, parameterwise
and joint shrinkage factor estimation. Journal of Statistical Software
69, 1-19.

Heinze G, Wallisch C, Dunkler D. Variable selection - A review and
recommendations for the practicing statistician. Biom J. 2018
May;60(3):431-449.

Johnson, R. W. (1996). Fitting percentage of body fat to simple body
measurements. Journal of Statistics Education, 4(1), 265–266.

Siri, W. E. (1956). The gross composition of the body. Advances in
Biological and Medical Physics, 4, 239–280.

# Session Information

``` r
sessionInfo()
```

    ## R version 4.4.2 (2024-10-31)
    ## Platform: aarch64-apple-darwin20
    ## Running under: macOS Sonoma 14.3.1
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: America/New_York
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] dplyr_1.1.4   corrplot_0.95 pander_0.6.5  shrink_1.2.3 
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] gtable_0.3.6        xfun_0.49           ggplot2_3.5.1      
    ##  [4] htmlwidgets_1.6.4   lattice_0.22-6      numDeriv_2016.8-1.1
    ##  [7] vctrs_0.6.5         tools_4.4.2         generics_0.1.3     
    ## [10] mfp_1.5.4.1         sandwich_3.1-1      tibble_3.2.1       
    ## [13] fansi_1.0.6         highr_0.11          cluster_2.1.6      
    ## [16] pkgconfig_2.0.3     Matrix_1.7-1        data.table_1.16.2  
    ## [19] checkmate_2.3.2     lifecycle_1.0.4     compiler_4.4.2     
    ## [22] stringr_1.5.1       MatrixModels_0.5-3  munsell_0.5.1      
    ## [25] codetools_0.2-20    SparseM_1.84-2      quantreg_5.99      
    ## [28] htmltools_0.5.8.1   yaml_2.3.10         htmlTable_2.4.3    
    ## [31] Formula_1.2-5       pillar_1.9.0        MASS_7.3-61        
    ## [34] rms_6.8-2           Hmisc_5.2-0         rpart_4.1.23       
    ## [37] multcomp_1.4-26     nlme_3.1-166        tidyselect_1.2.1   
    ## [40] digest_0.6.37       polspline_1.1.25    mvtnorm_1.3-2      
    ## [43] stringi_1.8.4       splines_4.4.2       fastmap_1.2.0      
    ## [46] grid_4.4.2          colorspace_2.1-1    cli_3.6.3          
    ## [49] magrittr_2.0.3      base64enc_0.1-3     survival_3.7-0     
    ## [52] utf8_1.2.4          TH.data_1.1-2       foreign_0.8-87     
    ## [55] scales_1.3.0        backports_1.5.0     rmarkdown_2.29     
    ## [58] nnet_7.3-19         gridExtra_2.3       zoo_1.8-12         
    ## [61] evaluate_1.0.1      knitr_1.48          rlang_1.1.4        
    ## [64] Rcpp_1.0.13-1       glue_1.8.0          rstudioapi_0.17.1  
    ## [67] R6_2.5.1
