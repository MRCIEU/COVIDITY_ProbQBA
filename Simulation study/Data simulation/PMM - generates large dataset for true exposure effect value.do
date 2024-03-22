/*************************************************************************************************************************************
                               COMPUTING THE EXPOSURE EFFECT ESTIMATE FOR A LARGE DATASET SIMULATED 
							             USING THE PMM DATA GENERATING MODEL

Logistic regression                                 Number of obs = 50,000,000
                                                    LR chi2(5)    = 2855501.39
                                                    Prob > chi2   =     0.0000
Log likelihood = -8867606.3                         Pseudo R2     =     0.1387

------------------------------------------------------------------------------
    infected | Coefficient  Std. err.      z    P>|z|     [95% conf. interval]
-------------+----------------------------------------------------------------
       sdbmi |   1.094689   .0007344  1490.53   0.000      1.09325    1.096129
         men |   .3109957    .001334   233.13   0.000     .3083812    .3136103
       sdage |   .0425092   .0006693    63.51   0.000     .0411974     .043821
      degree |  -.3137361   .0015862  -197.79   0.000     -.316845   -.3106272
      smoker |   .2089741   .0020903    99.97   0.000     .2048772     .213071
       _cons |  -3.495863   .0012446 -2808.90   0.000    -3.498302   -3.493424
------------------------------------------------------------------------------
Note: _cons estimates baseline odds.


Logistic regression                                 Number of obs = 50,000,000
                                                    LR chi2(5)    =  130711.37
                                                    Prob > chi2   =     0.0000
Log likelihood = -10057423                          Pseudo R2     =     0.0065

------------------------------------------------------------------------------
    infected | Coefficient  Std. err.      z    P>|z|     [95% conf. interval]
-------------+----------------------------------------------------------------
       sdbmi |   -.000208   .0006502    -0.32   0.749    -.0014823    .0010663
         men |   .3120372   .0012934   241.25   0.000     .3095022    .3145723
       sdage |   .0432621   .0006498    66.58   0.000     .0419886    .0445356
      degree |  -.3161138   .0014763  -214.13   0.000    -.3190073   -.3132204
      smoker |   .2109778   .0019497   108.21   0.000     .2071564    .2147991
       _cons |  -3.001566   .0010321 -2908.32   0.000    -3.003589   -2.999544
------------------------------------------------------------------------------
Note: _cons estimates baseline odds.
*************************************************************************************************************************************/
* SET WORKING DIRECTORY
	* DO FILES STORED IN "$mydirectory\Do files\"
	* SIMULATED DATASETS STORED IN "$mydirectory\Simulated data\"
cd "$mydirectory"

/********************
  NOT NULL SCENARIO
********************/
 set seed 546

local causal 1		
local samplesize 50000000
run "Do files\PMM - sets true values of data generation models.do" `causal'

* p(infected│sdbmi,smoker,degree,sdage,men)
local obs_asthma_infected = 1.058
local obs_diabetes_infected = .668
local obs_hyperten_infected = .451
foreach name in hyperten_infected diabetes_infected asthma_infected {
		global true_`name' = `obs_`name''
}

run "Do files\PMM - simulates a dataset.do" `samplesize'

* FIT SUBSTANTIVE ANALYSIS OF INTEREST
logit infected sdbmi men sdage degree smoker 

/********************
  NULL SCENARIO
********************/
set seed 546

local causal 0		
local samplesize 50000000
run "Do files\PMM - sets true values of data generation models.do" `causal'

* p(infected│sdbmi,smoker,degree,sdage,men)
local obs_asthma_infected = 1.058
local obs_diabetes_infected = .668
local obs_hyperten_infected = .451
foreach name in hyperten_infected diabetes_infected asthma_infected {
		global true_`name' = `obs_`name''
}

run "Do files\PMM - simulates a dataset.do" `samplesize'

* FIT SUBSTANTIVE ANALYSIS OF INTEREST
logit infected sdbmi men sdage degree smoker 

* END OF DO-FILE