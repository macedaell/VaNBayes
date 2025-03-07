NOTE: 

This is a new addition to the VaNBayes_v1 paper, added as a response to the reviewers asking for comparisons with other more established methods. 






This will be the VaNBayes implementation of a non-spatial SIR model, which will be compared to BayesFlow's estimate of the non-spatial SIR model. 

This is a work-in-progress, and currently on hold until BayesFlow gets good diagnostics. 

When Bayesflow has good diagnostics, do the following: 

- compare the prior distributions of each of the parameters to ensure those are the same

- compare the likelihood distributions of each of the datapoints to ensure those are the same (compare I1, I2, etc.)

- find a way to get simulation study datasets into our python


What we need to compare: 
- MAD and COV of each of the methods (can be done independently
