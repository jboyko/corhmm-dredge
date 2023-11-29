# corhmm-dredge
Automatically search discrete character evolution models and regularize the parameter values

## Simulation study and associated code
Each simulation setting has an identifier following a '-' and has a corresponding dataset in /data/ and model fit in /fits/. I vary several parameters for each script in order to test a variety of empirically relavent scenarios. Rates of evolution and tree size are two of the varried parameters. We are also interested in detecting particular model types when the data structure allows. For example, we may be interested in testing for correlation when there are two or more characters in the dataset. Finally we need to vary the number of rate classes. 
**code/simulate-01:** nChar=1, nStates=c(2), qRates=list(c(0.01,0.01), c(0.1, 0.1), c(1,1))