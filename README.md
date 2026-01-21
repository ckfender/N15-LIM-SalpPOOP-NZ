# N15-LIM-SalpPOOP-NZ
This code is largely a modified use case of Stukel et al 2022's N15 Bloofinz Gulf of Mexico script (https://github.com/mstukel/N15-LIM-BLOOFINZ-GoM). Like that project, the code here was written to solve a linear inverse ecosystem model that includes non-linear approximate equality constraints (the δ15N of each model compartment). It uses a Markov Chain Monte Carlo approach adapted from the MCMC+15N approach developed in Stukel et al. (2018a). The MCMC+15N, approach was itself an adaptation of the MCMC approach of Van Den Meersche et al. (2009). The MCMC+15N approach added a second simultaneous random walk through δ15N space. The δ15N vector in turn affected the approximate equalities (Ax≈b) in Van Den Meersche et al. (Van den Meersche et al., 2009). For a thorough discussion of the MCMC+15N approach, I refer readers to Stukel et al. (2018a; 2018b). This version of the code applies a similar topology to that of the Bloofinz project at the lower trophic levels, diverging at the macrozooplankton to better resolve the food web of New Zealand's Chatham Rise. Particular emphasis is placed on salps and the flows connected to them, especially their grazing, export, and consumption by higher trophic levels.

The goal of this code was to develop an approach sufficient to constrain an underdetermined ecosystem model for the open-ocean Chatham Rise with a specific focus on trophic pathways of salps within the euphotic zone (upper 70 meters). The model includes a total of 140 unknown ecosystem nitrogen flows to be solved for as well as 24 mass balance (exact equality) equations, 64 approximate equality constraints (including 20 field rate measurements, 25 nitrogen istope mass balance equations, and 19 fish ingestion estimates), and 138 inequality constraints. Full details of the model are available in ---. If you do not have access to that manuscript, please contact me for a copy.

Files

N15InverseModelNZ.xlsx and N15InverseModelNZC1: These files contain the model structure in an easily readable excel format. Note that because of slight differences in measurements Cycle 1 and the rest of the cruise (no evidence of migratory mesozooplankton in C1), structure differs slightly between the two configurations as there is one less compartment in the C1 version. This is also why there are often 2 versions of the following codes as well.

InputsNZ.xlsx and InputsNZC1: This excel file includes all the field data used to constrain the inverse model.

SetMatricesN15SalpPoop.m and SetMatricesN15SalpPoopC1.m: These Matlab files can be used to read the aforementioned excel files and create the matrices needed to run the LIM-N15-SalpPOOP-NZ model. This will creat initialization files for each of the 5 cycles following a naming convention of, in the case of Cycle 1, N15NZInverseCycle1.mat.

N15NZInverseCycleX.mat: These files (created by SetMatricesN15SalpPoop.m and SetMatricesN15SalpPoopC1.m) contain the matrices for Cycle X that must be solved by the LIM+15N procedure, as well as the δ15N values that add additional constraint to the random walk.

RunN15InverseNZ.R: This code runs the actual inverse model in R. It relies on the presence of the above input files (the .mat files mentioned previously) as well as the xsampleN15outputsNZ.R, ExternalFunctionsNZ.R, and ExternalFunctionsNZC1.R files. It will create files named N15NZInverseCycleXRoutputs.mat that contain the model solutions. These solutions are stored in the variables MCMCmat and del15N. Please note that the model is configured as it was run for the paper, including simulations >100,000,000 runs long which will likely take days to weeks to run in parallel for each cycle depending on your computational power.

xsampleN15outputsNZ.R: This is the file that actually runs the MCMC+15N inverse. It is called by RunN15InverseNZ.R.

ExternalFunctionsNZ.R and ExternalFunctionsNZC1: These files are called by xsampleN15outputsNZ.R and serve to re-set the approximate equality matrix (Aa) based on changes in the estimated δ15N of each model compartment as stored in the d15N vector.

Additional scripts related to post-processing and the creation of the figures used in the manuscripts are also provided in the Extras folder.

References:
Stukel, M.R., Décima, M. and Kelly, T.B., 2018a. A new approach for incorporating 15N isotopic data into linear inverse ecosystem models with Markov Chain Monte Carlo sampling. PloS one, 13(6): e0199123. 

Stukel, M.R., Décima, M., Landry, M.R. and Selph, K.E., 2018b. Nitrogen and isotope flows through the Costa Rica Dome upwelling ecosystem: The crucial mesozooplankton role in export flux. Global Biogeochemical Cycles, 32: 1815-1832. 

Stukel, M.R., Gerard, T., Kelly, T.B., Knapp, A. N., Laiz-Carrión, R., Lamkin, J. T., ... & Swalethorp, R. (2022). Plankton food webs in the oligotrophic Gulf of Mexico spawning grounds of Atlantic bluefin tuna. Journal of Plankton Research, 44(5), 763-781.

Van den Meersche, K., Soetaert, K. and Van Oevelen, D., 2009. xSample(): An R function for sampling linear inverse problems. Journal of Statistal Software, Code Snippets, 30(1): 1-15.
