/*********************************************************************************************************************************
                         SET THE TRUE VALUES OF THE DATA GENERATION MODEL UNDER THE PATTERN MIXTURE MODEL FRAMEWORK
							

(1) THE JOINT DISTRIBUTION OF THE COMPLETE DATA AND itest IS FACTORISED AS FOLLOWS: 								  
	p(men,sdage,degree,smoker,BMI,itest, infected, asthma,diabetes,hyperten)=
		p(men)p(age│men)p(degree│age,men)p(smoker│degree,sdage,men)
		p(sdbmi│smoker,degree,sdage,men)p(itest|sdbmi,smoker,degree,sdage,men)
		p(infected│itest,sdbmi,smoker,degree,sdage,men)
		p(asthma│infected,itest,sdbmi,smoker,degree,sdage,men)  
		p(diabetes│asthma,infected,itest,sdbmi,smoker,degree,sdage,men)
		p(hyperten│diabetes,asthma,infected,itest,sdbmi,smoker,degree,sdage,men)

(3) MISSING DATA MECHANISM FOR hyperten, sdbmi, AND smoker
	p(R_hyperten|sdage,men,degree,asthma,diabetes)
	p(R_smoker|sdage,men,degree,asthma,diabetes)
	p(R_smoker|sdage,men,degree,asthma,diabetes)

ARGUMENTS OF THE DO-FILE
	causal
	        !=1 SETS true_infected_sdbmi=0; 
			==1 SETS true_infected_sdbmi=log(3)
***********************************************************************************************************************************/
args causal

/******************************************************************************
  (1) SET TRUE VALUES OF THE JOINT DISTRIBUTION OF THE COMPLETE DATA AND itest
*******************************************************************************/
* p(men)
global true_pr_male = 0.457
global true_men_cons = logit($true_pr_male)		

* p(sdage│men)
global true_sdage_cons = -.021
global true_sdage_men = .046
global true_sdage_rmse = 1

* p(degree│sdage,men)
global true_degree_cons = -.795
global true_degree_men = .14
global true_degree_sdage = -.272

* p(smoker│degree,sdage,men)
global true_smoker_cons = -2.221
global true_smoker_men = .406
global true_smoker_sdage = -.303
global true_smoker_degree = -.6

* p(sdbmi│smoker,degree,sdage,men)
global true_sdbmi_cons = .018
global true_sdbmi_men = .168
global true_sdbmi_sdage = .024
global true_sdbmi_degree = -.26
global true_sdbmi_smoker = -.123
global true_sdbmi_rmse = .987

* PARAMETER VALUES THAT DIFFER ACCORDING TO NOT NULL AND NULL SCENARIOS
if `causal'==1 {
	* p(T|X,C)
	global true_itest_cons = -3.4586145
	global true_itest_men = .30176466
	global true_itest_sdage = .05649763
	global true_itest_degree = -.2985004
	global true_itest_smoker = .21488435
	global true_itest_sdbmi = 1.0381909

	* p(Y|X,C,T)
	global true_infected_cons = -5.9713074
	global true_infected_men = .22224755
	global true_infected_sdage = -.06950772
	global true_infected_degree = -.26313165
	global true_infected_smoker = .08290057
	global true_infected_sdbmi = .94515627
	global true_infected_itest = 8.1205697

	* p(Aa|X,C,T,Y)	
	global true_asthma_cons = -1.1900994
	global true_asthma_men = .06807318
	global true_asthma_sdage = .13082294
	global true_asthma_degree = -.09155731
	global true_asthma_smoker = .41412888
	global true_asthma_sdbmi = .1875341
	global true_asthma_itest = 1.063066
	global true_asthma_infected = .64013071

	* p(Ad|X,C,T,Y,Aa)
	global true_diabetes_cons = -3.0968974
	global true_diabetes_men = .61559818
	global true_diabetes_sdage = .42560668
	global true_diabetes_degree = -.25613796
	global true_diabetes_smoker = .33248871
	global true_diabetes_sdbmi = .75632042
	global true_diabetes_itest = .60952188
	global true_diabetes_infected = .65927736
	global true_diabetes_asthma = .48836744

	* p(Ah|X,C,T,Y,Aa,Ad)
	global true_hyperten_cons = -.89582443
	global true_hyperten_men = .34913013
	global true_hyperten_sdage = .64914671
	global true_hyperten_degree = -.19932663
	global true_hyperten_smoker = .09417527
	global true_hyperten_sdbmi = .52706477
	global true_hyperten_itest = .43888419
	global true_hyperten_infected = .61158019
	global true_hyperten_asthma = .37620764
	global true_hyperten_diabetes = 1.3718944
}
else {	
	* p(T|X,C)
	global true_itest_sdbmi = .02435542
	global true_itest_men = .29608111
	global true_itest_sdage = .05733592
	global true_itest_degree = -.2956811
	global true_itest_smoker = .21216036
	global true_itest_cons = -3.0138464
	
	* p(Y|X,C,T)
	global true_infected_cons = -5.284835
	global true_infected_men = .24745199
	global true_infected_sdage = -.04659723
	global true_infected_degree = -.27876376
	global true_infected_smoker = .11609946
	global true_infected_sdbmi = -.1182128
	global true_infected_itest = 7.5952799

	* p(Aa|X,C,T,Y)	
	global true_asthma_cons = -1.1918436
	global true_asthma_men = .06834257
	global true_asthma_sdage = .13095646
	global true_asthma_degree = -.09107356
	global true_asthma_smoker = .41405647
	global true_asthma_sdbmi = .18758709
	global true_asthma_itest = 1.0511195
	global true_asthma_infected = -.96697196

	* p(Ad|X,C,T,Y,Aa)
	global true_diabetes_cons = -3.0989876
	global true_diabetes_men = .61560474
	global true_diabetes_sdage = .42576045
	global true_diabetes_degree = -.25614707
	global true_diabetes_smoker = .33244211
	global true_diabetes_sdbmi = .75627339
	global true_diabetes_itest = .63080886
	global true_diabetes_infected = -.583318
	global true_diabetes_asthma = .4906264

	* p(Ah|X,C,T,Y,Aa,Ad)
	global true_hyperten_cons = -.89607342
	global true_hyperten_men = .34883197
	global true_hyperten_sdage = .64909335
	global true_hyperten_degree = -.19945855
	global true_hyperten_smoker = .09398336
	global true_hyperten_sdbmi = .52711023
	global true_hyperten_itest = .44045975
	global true_hyperten_infected = -.39637508
	global true_hyperten_asthma = .37577773
	global true_hyperten_diabetes = 1.3714082
}

/**************************************************************************************
  (2) SET TRUE VALUES OF THE MISSING DATA MECHANISM FOR BMI, HYPERTENSION AND SMOKER
**************************************************************************************/
* p(Rsdbmi|men, sdage, degree, asthma, diabetes)
global true_Rsdbmi_men = .233
global true_Rsdbmi_sdage = -.057
global true_Rsdbmi_degree = -.133
global true_Rsdbmi_asthma = .363
global true_Rsdbmi_diabetes = .763
global true_Rsdbmi_cons -3.2
 
* MISSING DATA MECHANISM FOR hyperten - DEPENDS ON FULLY OBSERVED VARIABLES
global true_Rhyperten_men = -.059
global true_Rhyperten_sdage = -.029
global true_Rhyperten_degree = -.19
global true_Rhyperten_asthma = .13
global true_Rhyperten_diabetes = .192
global true_Rhyperten_cons = -2.953
 
* MISSING DATA MECHANISM FOR smoker - DEPENDS ON FULLY OBSERVED VARIABLES
global true_Rsmoker_men = .072
global true_Rsmoker_sdage = .232
global true_Rsmoker_degree = -.774
global true_Rsmoker_asthma = .169
global true_Rsmoker_diabetes = .417
global true_Rsmoker_cons = -2.9

* END OF DO FILE