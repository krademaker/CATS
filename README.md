# Classification Assessment of Tumor Subtypes (CATS)

## Background
This was part of the course '[Bioinformatics for Translational Medicine](https://studiegids.vu.nl/en/2019-2020/courses/X_405092)', taught at the Vrije Universiteit Amsterdam during the MSc Bioinformatics and Systems Biology.

Breast cancer is a heterogeneous disease and classification of breast cancer tumors in their molecular
subtypes has important implications for the treatment and prognosis. Three receptors play a pivotal role in
these subtypes: the Estrogen Receptor (ER), Progesterone Receptor (PR) and Human Epidermal growth
factor Receptor 2 (HER2). After removal of the tumor, the pathology department of the hospital tests these
samples for presence of ER, PR and HER2. The three main subtypes in breast cancer on which the
treatment decision will be based are:
- HER2 positive: HER2+
- Hormone receptor positive (HR+): ER+ and/or PR+, and HER2-
- Triple negative (TN): ER-, PR- and HER2-
Each of the three subtypes reacts differently to different types of treatment.

In this assignment, we trained a classifier that is able to predict subtypes of breast cancer
tumors based on array CGH (aCGH) data.

## About our work
You can read about our work in the paper that we wrote [over here](https://github.com/krademaker/CATS/blob/master/CATS-Classifying%20breast%20cancer%20subgroups%20with%20chromosomal%20aberration%20patterns.pdf).

In summary, our _best-performing model was a single hidden layer feed-forward neural network (SLFN) with Boruta_ as feature selection method (38 features selected) that had an _accuracy score of 0.895_.

In the competetion that was organized for this course, **we came in first with 47 correct predictions out of 57 samples** (82%, validation in `/data` folder, results unavailable to us.)

## Members
- Will Harley
- Aamir Hasan
- Saarika Prathivadi Bhayankaram
- Koen Rademaker
- Joris Visser
