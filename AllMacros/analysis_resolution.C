/* 
   Winston DeGraw (wdegraw@lbl.gov)
   Rey Cruz-Torres (reynier@lbl.gov)

   Heavily modified by Karthik Suresh (ksuresh@wm.edu)
*/


// Forward-declaring functions
void prettyTH1F( TH1F * h1 , int color , int marker , float min , float max );
//

Double_t DoubleGauss(Double_t *x, Double_t *par)
{
	return par[0]*exp(-0.5*pow((x[0]-par[1])/par[2], 2)) + par[3]*exp(-0.5*pow((x[0]-par[1])/par[4], 2));
}

//============================================================================================================================================
void analysis_resolution(string FileName, bool isDefault = 0, float Bfield = 1.4){
	// -------------------------------------------------------------
	// Some important parameters
	const TString partic = "pi-";	// particle to be studied
	//const float Bfield = 3.0;	// [T] Magnetic field
	bool use_widths = true;
	bool update_tab = true;
	// -------------------------
	// Binning
	//float eta_bin[] = {0.,0.5,1.0,1.5,2.0,2.5,3.0,3.5};		const int size_eta_bin = sizeof(eta_bin)/sizeof(*eta_bin);
	//float eta_bin[] = {0., 1.};
	//float eta_bin[] = {-3.5, -3.0, -2.5, 2.0, -1.5, -1.0, -0.5, 0., 0.5, 1., 1.5, 2.0, 2.5, 3.0, 3.5};               const int size_eta_bin = sizeof(eta_bin)/sizeof(*eta_bin);
	float eta_bin[] = {-3.4, -2., -1.0, 1.0, 2., 3.4};               const int size_eta_bin = sizeof(eta_bin)/sizeof(*eta_bin);
	float mom_bin[] = {1., 2., 4., 6., 8., 10., 12., 14., 16., 18., 20.};	const int size_mom_bin = sizeof(mom_bin)/sizeof(*mom_bin);
	float pT_bin[] = {0.,2.,4.,6.,8.,10.,12.,14.,16.,18.,20.}; const int size_pT_bin = sizeof(pT_bin)/sizeof(*pT_bin);
	// ------------------------- Bfield is only used to name the files
	string Bfield_str = Form("%.1f",Bfield);
	Bfield_str.replace(Bfield_str.find("."), sizeof(".") - 1, "_");
	TString tab_name = Form("sigma_eta_%i_p_%i_",size_eta_bin-1,size_mom_bin-1) + partic + Form("_B%.1fT.txt",Bfield);
	// -------------------------------------------------------------
	// Some settings
	TH1::SetDefaultSumw2();
	TH2::SetDefaultSumw2();
	gStyle -> SetOptStat(0);	
	// -------------------------------------------------------------
	// Loading all the needed info from the root file
	TFile * F = new TFile(FileName.data());
	TTree * T = (TTree*) F -> Get("tracks");
	int trackID;
	float gpx, gpy, gpz, px, py, pz;
	float mRICHProj_proj_px, mRICHProj_proj_py, mRICHProj_proj_pz, mRICHProj_px, mRICHProj_py, mRICHProj_pz;
	float DIRCProj_proj_px, DIRCProj_proj_py, DIRCProj_proj_pz, DIRCProj_px, DIRCProj_py, DIRCProj_pz;
	float dRICHProj_proj_px, dRICHProj_proj_py, dRICHProj_proj_pz, dRICHProj_px, dRICHProj_py, dRICHProj_pz;
	float dca2d;
	T -> SetBranchAddress("trackID",&trackID);
	T -> SetBranchAddress("gpx",&gpx);
	T -> SetBranchAddress("gpy",&gpy);
	T -> SetBranchAddress("gpz",&gpz);
	T -> SetBranchAddress("px" ,&px );
	T -> SetBranchAddress("py" ,&py );
	T -> SetBranchAddress("pz" ,&pz );
	T -> SetBranchAddress("dca2d",&dca2d); // Added on June 11 2021 to study the dca2d

	int nEntries = T -> GetEntries();
	// -------------------------------------------------------------
	fstream tab;
	float approx_sig_dpp[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float approx_sig_dth[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float approx_sig_dph[size_eta_bin-1][size_mom_bin-1] = {{0}};
	TString temp_str;
	if(use_widths){
		tab.open("tables/"+tab_name);
		if(!tab){cout << "Could not find file '" << tab_name << "'" << endl; use_widths = false; update_tab = true;}
		else{
			cout << "Loading parameters from file '" << tab_name << "'" << endl;
			for(int et = 0 ; et < size_eta_bin-1 ; et++){ for(int p = 0 ; p < size_mom_bin-1 ; p++){tab >> approx_sig_dpp[et][p];}}	//tab >> temp_str;
			for(int et = 0 ; et < size_eta_bin-1 ; et++){ for(int p = 0 ; p < size_mom_bin-1 ; p++){tab >> approx_sig_dth[et][p];}}	//tab >> temp_str;
			for(int et = 0 ; et < size_eta_bin-1 ; et++){ for(int p = 0 ; p < size_mom_bin-1 ; p++){tab >> approx_sig_dph[et][p];}}
		}
		tab.close();
	}

	float approx_sig_dpp_3_0[size_eta_bin-1][size_mom_bin-1] = {0}; float approx_sig_dpp_1_2[size_eta_bin-1][size_mom_bin-1] = {0};
	float approx_sig_dth_3_0[size_eta_bin-1][size_mom_bin-1] = {0}; float approx_sig_dth_1_2[size_eta_bin-1][size_mom_bin-1] = {0};
	float approx_sig_dph_3_0[size_eta_bin-1][size_mom_bin-1] = {0}; float approx_sig_dph_1_2[size_eta_bin-1][size_mom_bin-1] = {0};

	for(int et = 0 ; et < size_eta_bin-1 ; et++){
		for(int p = 0 ; p < size_mom_bin-1 ; p++){
			approx_sig_dpp_3_0[et][p] = 3.0*approx_sig_dpp[et][p]; approx_sig_dpp_1_2[et][p] = 1.2*approx_sig_dpp[et][p];
			approx_sig_dth_3_0[et][p] = 3.0*approx_sig_dth[et][p]; approx_sig_dth_1_2[et][p] = 1.2*approx_sig_dth[et][p];
			approx_sig_dph_3_0[et][p] = 3.0*approx_sig_dph[et][p]; approx_sig_dph_1_2[et][p] = 1.2*approx_sig_dph[et][p];
		}
	}
	// -------------------------------------------------------------
	// Defining histograms
	TH1F *** h1_dpp_p_et_bins = new TH1F**[size_eta_bin-1];	// delta p / p vs. p in eta bins
	TH1F *** h1_dth_p_et_bins = new TH1F**[size_eta_bin-1];	// delta theta vs. p in eta bins
	TH1F *** h1_dph_p_et_bins = new TH1F**[size_eta_bin-1];	// delta phi   vs. p in eta bins
	TH1F *** h1_dca2d_p_et_bins = new TH1F**[size_eta_bin-1]; // dca2d vs. p in eta bins
	TH1F *** h1_dca2d_pT_et_bins = new TH1F**[size_eta_bin-1]; // dca2d vs. pT in eta bins
	for(int et = 0 ; et < size_eta_bin-1 ; et++){
		h1_dpp_p_et_bins[et] = new TH1F*[size_mom_bin-1];
		h1_dth_p_et_bins[et] = new TH1F*[size_mom_bin-1];
		h1_dph_p_et_bins[et] = new TH1F*[size_mom_bin-1];
		h1_dca2d_p_et_bins[et] = new TH1F*[size_mom_bin-1];
		h1_dca2d_pT_et_bins[et] = new TH1F*[size_mom_bin-1];
		for(int p = 0 ; p < size_mom_bin-1 ; p++){
			if(use_widths){
				h1_dpp_p_et_bins[et][p] = new TH1F(Form("h1_dpp_p_et_bins_%i_%i",et,p),";dp/p;Counts"         ,100,-approx_sig_dpp_3_0[et][p],approx_sig_dpp_3_0[et][p]);
				h1_dth_p_et_bins[et][p] = new TH1F(Form("h1_dth_p_et_bins_%i_%i",et,p),";d#theta [rad];Counts",100,-approx_sig_dth_3_0[et][p],approx_sig_dth_3_0[et][p]);
				h1_dph_p_et_bins[et][p] = new TH1F(Form("h1_dph_p_et_bins_%i_%i",et,p),";d#phi [rad];Counts"  ,100,-approx_sig_dph_3_0[et][p],approx_sig_dph_3_0[et][p]);
			}
			else{

				double pbins = 500;
				if(fabs(eta_bin[et]) >=1. && fabs(eta_bin[et] < 2.))
                                {
                                        pbins = 400;
                                }
				if(fabs(eta_bin[et]) >=2.5)
				{
					pbins = 300;
				}
                                if(fabs(eta_bin[et]) >=3.0)
                                {
                                        pbins = 250;
                                }




				h1_dpp_p_et_bins[et][p] = new TH1F(Form("h1_dpp_p_et_bins_%i_%i",et,p),";dp/p;Counts" ,pbins, -1, 1);
				
				

		// The dph, dth, and dca2d are not defined as roboust as p. It is customised though after seeing into a number of different distributions corresponding to different design points

				double ph_range = 0.005; // This if or projection phi at PID detectors
				/*
				if(fabs(eta_bin[et]) >= 2.5)
				{
					if(mom_bin[p] <= 4) ph_range = 0.05;
					else ph_range = 0.025;
				}
				if(fabs(eta_bin[et]) >= 1 && fabs(eta_bin[et]) < 2.5)
				{
					if(mom_bin[p] <= 8) ph_range = 0.01;
					else ph_range = 0.005;
				}*/
				if(mom_bin[p] <= 4) ph_range = 0.025;
				h1_dph_p_et_bins[et][p] = new TH1F(Form("h1_dph_p_et_bins_%i_%i",et,p),";d#phi [rad];Counts"  ,200,-1*ph_range  ,ph_range  );


				double th_range = 0.0014;
				int th_bins = 200;
				if(mom_bin[p] < 4)
				{
					th_range = 0.005;
					th_bins = 100;
				}
				h1_dth_p_et_bins[et][p] = new TH1F(Form("h1_dth_p_et_bins_%i_%i",et,p),";d#theta [rad];Counts", th_bins,-1.*th_range,1.*th_range);

				double dca2d_range = 0.25;
				int bins = 500;
				if(fabs(eta_bin[et]) < 1.0) 
					{
						dca2d_range = 0.01;
						if(mom_bin[p] < 4.) 
						{
							dca2d_range = 0.05;
							bins = 100;
						}
					}
				else if(fabs(eta_bin[et]) >=1.0 && fabs(eta_bin[et]) < 2.5) 
					{
						if(mom_bin[p] < 4.) 
						{
							dca2d_range = 0.06;
							bins = 100;						
						}
						else dca2d_range = 0.02;
					}
				else if(fabs(eta_bin[et]) >= 2.5 && fabs(eta_bin[et]) < 3.5) 
					{
						if(mom_bin[p] < 3)
						{
							dca2d_range = 0.5;
							bins = 100;
						}
						else if(mom_bin[p] >=3 && mom_bin[p] < 6.) 
						{
							dca2d_range = 0.1;
							bins = 100;
						}
						
						else dca2d_range = 0.06;
					}
				else dca2d_range = 0.1;
				

				h1_dca2d_p_et_bins[et][p] = new TH1F(Form("h1_dca2d_p_et_bins_%i_%i",et,p),";dca2d [cms];Counts"  ,bins,-dca2d_range  ,dca2d_range  );
				h1_dca2d_pT_et_bins[et][p] = new TH1F(Form("h1_dca2d_pT_et_bins_%i_%i",et,p),";dca2d [cms];Counts"  ,bins,-dca2d_range  ,dca2d_range  );
			}

			h1_dpp_p_et_bins[et][p] -> SetTitle(Form("%.1f < |#eta| < %.1f, %.1f < p < %.1f GeV/c",eta_bin[et],eta_bin[et+1],mom_bin[p],mom_bin[p+1]));
			h1_dth_p_et_bins[et][p] -> SetTitle(Form("%.1f < |#eta| < %.1f, %.1f < p < %.1f GeV/c",eta_bin[et],eta_bin[et+1],mom_bin[p],mom_bin[p+1]));
			h1_dph_p_et_bins[et][p] -> SetTitle(Form("%.1f < |#eta| < %.1f, %.1f < p < %.1f GeV/c",eta_bin[et],eta_bin[et+1],mom_bin[p],mom_bin[p+1]));
			h1_dca2d_p_et_bins[et][p] -> SetTitle(Form("%.1f < |#eta| < %.1f, %.1f < p < %.1f GeV/c",eta_bin[et],eta_bin[et+1],mom_bin[p],mom_bin[p+1]));
			h1_dca2d_pT_et_bins[et][p] -> SetTitle(Form("%.1f < |#eta| < %.1f, %.1f < pT < %.1f GeV/c",eta_bin[et],eta_bin[et+1],mom_bin[p],mom_bin[p+1]));
		}
	}	
	// -------------------------------------------------------------	
	int color[] = {1,2,62,8,95,52,6,28,209,92,15,1,2,62,8,95,52};
	int marker[] = {20,21,23,24,25,26,27,28,29,30,20,21,23,24,25,26};

	TH1F ** h1_dpp_v_p_et_bins = new TH1F*[size_eta_bin-1];
	TH1F ** h1_dth_v_p_et_bins = new TH1F*[size_eta_bin-1];
	TH1F ** h1_dph_v_p_et_bins = new TH1F*[size_eta_bin-1];
	TH1F ** h1_dca2d_v_p_et_bins = new TH1F*[size_eta_bin-1];
	TH1F ** h1_dca2d_v_pT_et_bins = new TH1F*[size_eta_bin-1];

	for(int et = 0 ; et < size_eta_bin-1 ; et++){
		h1_dpp_v_p_et_bins[et] = new TH1F(Form("h1_dpp_v_p_et_bins_%i",et),";p [GeV/c];dp/p [%]"      ,size_mom_bin-1,mom_bin);	prettyTH1F( h1_dpp_v_p_et_bins[et] , color[et] , marker[et] , 0. , 10. );
		h1_dth_v_p_et_bins[et] = new TH1F(Form("h1_dth_v_p_et_bins_%i",et),";p [GeV/c];d#theta [mrad]",size_mom_bin-1,mom_bin);	prettyTH1F( h1_dth_v_p_et_bins[et] , color[et] , marker[et] , 0. , 1.  );
		h1_dph_v_p_et_bins[et] = new TH1F(Form("h1_dph_v_p_et_bins_%i",et),";p [GeV/c];d#phi [mrad]"  ,size_mom_bin-1,mom_bin);	prettyTH1F( h1_dph_v_p_et_bins[et] , color[et] , marker[et] , 0. , 25. );
		h1_dca2d_v_p_et_bins[et] = new TH1F(Form("h1_dca2d_v_p_et_bins_%i",et),Form("%f < |#eta| < %f;p [GeV/c];dca2d [ums]", eta_bin[et], eta_bin[et+1]) ,size_mom_bin-1,mom_bin);	prettyTH1F( h1_dca2d_v_p_et_bins[et] , color[et] , marker[et] , 0. , 5000. );
		h1_dca2d_v_pT_et_bins[et] = new TH1F(Form("h1_dca2d_v_pT_et_bins_%i",et),Form("%f < |#eta| < %f;pT [GeV/c];dca2d [ums]", eta_bin[et], eta_bin[et+1]) ,size_mom_bin-1,mom_bin);	prettyTH1F( h1_dca2d_v_pT_et_bins[et] , color[et] , marker[et] , 0. , 5000. );
	}

	TH1F ** h1_dpp_v_et_p_bins = new TH1F*[size_mom_bin-1];
	TH1F ** h1_dth_v_et_p_bins = new TH1F*[size_mom_bin-1];
	TH1F ** h1_dph_v_et_p_bins = new TH1F*[size_mom_bin-1];
	TH1F ** h1_dca2d_v_et_p_bins = new TH1F*[size_mom_bin-1];
	TH1F ** h1_dca2d_v_et_pT_bins = new TH1F*[size_mom_bin-1];

	for(int p = 0 ; p < size_mom_bin-1 ; p++){
		h1_dpp_v_et_p_bins[p] = new TH1F(Form("h1_dpp_v_et_p_bins_%i",p),";#eta;dp/p [%]"      ,size_eta_bin-1,eta_bin);	prettyTH1F( h1_dpp_v_et_p_bins[p] , color[p] , marker[p] , 0. , 10. );
		h1_dth_v_et_p_bins[p] = new TH1F(Form("h1_dth_v_et_p_bins_%i",p),";#eta;d#theta [mrad]",size_eta_bin-1,eta_bin);	prettyTH1F( h1_dth_v_et_p_bins[p] , color[p] , marker[p] , 0. , 1.  );
		h1_dph_v_et_p_bins[p] = new TH1F(Form("h1_dph_v_et_p_bins_%i",p),";#eta;d#phi [mrad]"  ,size_eta_bin-1,eta_bin);	prettyTH1F( h1_dph_v_et_p_bins[p] , color[p] , marker[p] , 0. , 25. );
		h1_dca2d_v_et_p_bins[p] = new TH1F(Form("h1_dca2d_v_et_p_bins_%i",p),";#eta;dca2d [ums]"  ,size_eta_bin-1,eta_bin);	prettyTH1F( h1_dca2d_v_et_p_bins[p] , color[p] , marker[p] , 0. , 5000. );
		h1_dca2d_v_et_pT_bins[p] = new TH1F(Form("h1_dca2d_v_et_pT_bins_%i",p),";#eta;dca2d [ums]"  ,size_eta_bin-1,eta_bin);	prettyTH1F( h1_dca2d_v_et_pT_bins[p] , color[p] , marker[p] , 0. , 5000. );
	}

	// -------------------------------------------------------------
	// Declaring other useful variables and functions
	float width_dpp[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float error_dpp[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float width_dth[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float error_dth[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float width_dph[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float error_dph[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float width_dca2d[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float error_dca2d[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float width_dca2d_pT[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float error_dca2d_pT[size_eta_bin-1][size_mom_bin-1] = {{0}};


	// Declaring the chi2 and NDF for fits
	
	float Chi2_dpp[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float NDF_dpp[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float Chi2_dth[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float NDF_dth[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float Chi2_dph[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float NDF_dph[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float Chi2_dca2d[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float NDF_dca2d[size_eta_bin-1][size_mom_bin-1] = {{0}};

	float sigma1_dpp[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float sigma2_dpp[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float sigma1_dth[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float sigma2_dth[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float sigma1_dph[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float sigma2_dph[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float sigma1_dca2d[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float sigma2_dca2d[size_eta_bin-1][size_mom_bin-1] = {{0}};

	float A1_dpp[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float A2_dpp[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float A1_dth[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float A2_dth[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float A1_dph[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float A2_dph[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float A1_dca2d[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float A2_dca2d[size_eta_bin-1][size_mom_bin-1] = {{0}};

	float A1_inf[size_eta_bin-1][size_mom_bin-1] = {{0}}; // These correspond to the one multiplying in the PseudoKFInEfficiency
	float A2_inf[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float Entry_inf[size_eta_bin-1][size_mom_bin-1] = {{0}};

	bool OverFlowUnderFlowFlag[size_eta_bin-1][size_mom_bin-1] = {{0}};

	float PseudoGlobalKFInEff[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float error_PseudoGlobalKFInEff[size_eta_bin-1][size_mom_bin-1] = {{0}};


	float TotalCounts[size_eta_bin-1][size_mom_bin-1] = {{0}}; // this will be used to calculate the reconstruction efficiency https://indico.bnl.gov/event/11590/contributions/49428/attachments/34393/55804/ECCETrkWG_pingLANL_05102021.pdf for details
	float mean_dpp[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float Yield[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float Eff[size_eta_bin-1][size_mom_bin-1] = {{0}};
	
	float CrisEff[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float CrisTotalCounts[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float error_CrisEff[size_eta_bin-1][size_mom_bin-1] = {{0}};

	float KF_Failed_Counts[size_eta_bin-1][size_mom_bin-1] = {{0}};
	float error_KF_Eff[size_eta_bin-1][size_mom_bin-1] = {{0}};

	TF1 *** f_gaus_dpp = new TF1**[size_eta_bin-1];
	TF1 *** f_gaus1_dpp = new TF1**[size_eta_bin-1];
	TF1 *** f_gaus2_dpp = new TF1**[size_eta_bin-1];
	TF1 *** f_gaus_dth = new TF1**[size_eta_bin-1];
	TF1 *** f_gaus1_dth = new TF1**[size_eta_bin-1];
	TF1 *** f_gaus2_dth = new TF1**[size_eta_bin-1];
	TF1 *** f_gaus_dph = new TF1**[size_eta_bin-1];
	TF1 *** f_gaus1_dph = new TF1**[size_eta_bin-1];
	TF1 *** f_gaus2_dph = new TF1**[size_eta_bin-1];
	TF1 *** f_gaus_dca2d = new TF1**[size_eta_bin-1];
	TF1 *** f_gaus1_dca2d = new TF1**[size_eta_bin-1];
	TF1 *** f_gaus2_dca2d = new TF1**[size_eta_bin-1];
	TF1 *** f_gaus_dca2d_pT = new TF1**[size_eta_bin-1];

	for(int et = 0 ; et < size_eta_bin-1 ; et++){
		f_gaus_dpp[et] = new TF1*[size_mom_bin-1];
		f_gaus1_dpp[et] = new TF1*[size_mom_bin-1];
		f_gaus2_dpp[et] = new TF1*[size_mom_bin-1];

		f_gaus_dth[et] = new TF1*[size_mom_bin-1];
		f_gaus1_dth[et] = new TF1*[size_mom_bin-1];
		f_gaus2_dth[et] = new TF1*[size_mom_bin-1];

		f_gaus_dph[et] = new TF1*[size_mom_bin-1];
		f_gaus1_dph[et] = new TF1*[size_mom_bin-1];
		f_gaus2_dph[et] = new TF1*[size_mom_bin-1];

		f_gaus_dca2d[et] = new TF1*[size_mom_bin-1];
		f_gaus1_dca2d[et] = new TF1*[size_mom_bin-1];
		f_gaus2_dca2d[et] = new TF1*[size_mom_bin-1];

		f_gaus_dca2d_pT[et] = new TF1*[size_mom_bin-1];

		for(int p = 0 ; p < size_mom_bin-1 ; p++){
			if(use_widths){
				f_gaus_dpp[et][p] = new TF1(Form("f_gaus_dpp_%i_%i",et,p),"gaus",-approx_sig_dpp_1_2[et][p],approx_sig_dpp_1_2[et][p]);
				f_gaus_dth[et][p] = new TF1(Form("f_gaus_dth_%i_%i",et,p),"gaus",-approx_sig_dth_1_2[et][p],approx_sig_dth_1_2[et][p]);
				f_gaus_dph[et][p] = new TF1(Form("f_gaus_dph_%i_%i",et,p),"gaus",-approx_sig_dph_1_2[et][p],approx_sig_dph_1_2[et][p]);
			}
			else{
				//f_gaus_dpp[et][p] = new TF1(Form("f_gaus_dpp_%i_%i",et,p),"gaus",-10.*0.007 ,10.*0.007 );
				f_gaus_dpp[et][p] = new TF1(Form("f_gaus_dpp_%i_%i",et,p),DoubleGauss,-1 ,1, 5);
				f_gaus1_dpp[et][p] = new TF1(Form("f_gaus1_dpp_%i_%i",et,p),"gaus",-1 ,1);
				f_gaus2_dpp[et][p] = new TF1(Form("f_gaus2_dpp_%i_%i",et,p),"gaus",-1 ,1);
				f_gaus1_dpp[et][p]->SetLineColor(kGreen);
				f_gaus2_dpp[et][p]->SetLineColor(kYellow);


				//f_gaus_dth[et][p] = new TF1(Form("f_gaus_dth_%i_%i",et,p),"gaus",-2.*0.0007,2.*0.0007);
				f_gaus_dth[et][p] = new TF1(Form("f_gaus_dth_%i_%i",et,p),DoubleGauss,-0.01,0.01, 5);
				f_gaus1_dth[et][p] = new TF1(Form("f_gaus1_dth_%i_%i",et,p),"gaus",-0.01,0.01);
				f_gaus2_dth[et][p] = new TF1(Form("f_gaus2_dth_%i_%i",et,p),"gaus",-0.01,0.01);
				f_gaus1_dth[et][p]->SetLineColor(kGreen);
				f_gaus2_dth[et][p]->SetLineColor(kYellow);

				f_gaus_dph[et][p] = new TF1(Form("f_gaus_dph_%i_%i",et,p),DoubleGauss,-0.05  ,0.05, 5 );
				f_gaus1_dph[et][p] = new TF1(Form("f_gaus1_dph_%i_%i",et,p),"gaus",-0.05  ,0.05 );
				f_gaus2_dph[et][p] = new TF1(Form("f_gaus2_dph_%i_%i",et,p),"gaus",-0.05  ,0.05 );
				f_gaus1_dph[et][p]->SetLineColor(kGreen);
				f_gaus2_dph[et][p]->SetLineColor(kYellow);

				f_gaus_dca2d[et][p] = new TF1(Form("f_gaus_dca2d_%i_%i",et,p),DoubleGauss,-0.1  ,0.1, 5  );
				f_gaus1_dca2d[et][p] = new TF1(Form("f_gaus1_dca2d_%i_%i",et,p),"gaus",-0.1  ,0.1  );
				f_gaus2_dca2d[et][p] = new TF1(Form("f_gaus2_dca2d_%i_%i",et,p),"gaus",-0.1  ,0.1  );
				f_gaus1_dca2d[et][p]->SetLineColor(kGreen);
				f_gaus2_dca2d[et][p]->SetLineColor(kYellow);


				f_gaus_dca2d_pT[et][p] = new TF1(Form("f_gaus_dca2d_pT_%i_%i",et,p),"gaus",-0.1  ,0.1  );
			}
		}
	}
	// -------------------------------------------------------------
	// Loop over entries of the tree
	for(int ev = 0 ; ev < nEntries ; ev++){
		T -> GetEntry(ev);
		if(ev%1000000==0) cout << "Looping over entry " << ev << " out of " << nEntries << endl;
		//if(trackID<0) continue; // this is to make sure kalman filter had worked	
		// Calculating some variables
		float gtheta = TMath::ACos(gpz/sqrt(gpx*gpx+gpy*gpy+gpz*gpz));	
		float theta = TMath::ACos(pz/sqrt(px*px+py*py+pz*pz));
		float dth = theta - gtheta;

		float geta = -TMath::Log(TMath::Tan(gtheta/2.));

		float p_reco = sqrt(px*px+py*py+pz*pz);
		float p_truth = sqrt(gpx*gpx+gpy*gpy+gpz*gpz);
		float pT_truth = sqrt(gpx*gpx + gpy*gpy);
		float dp_p = (p_reco-p_truth)/p_truth;

		float gphi = TMath::ATan(gpy/gpx);
		float phi = TMath::ATan(py/px);
		float dph = phi - gphi;
		
		// Filling histograms
		// Here the looping is needed since we have uneven binning in eta It can still be improved. 
		for(int et = 0 ; et < size_eta_bin-1 ; et++){
			if( geta > eta_bin[et] &&  geta <= eta_bin[et+1] ){
				for(int p = 0 ; p < size_mom_bin-1 ; p++){
					if( p_truth > mom_bin[p] && p_truth <= mom_bin[p+1] ){
						if(trackID>=0){ // Make sure we only fill in passed KF Tracks. 
						h1_dpp_p_et_bins[et][p] -> Fill( dp_p );
						h1_dth_p_et_bins[et][p] -> Fill( dth  );
						h1_dph_p_et_bins[et][p] -> Fill( dph  );
						h1_dca2d_p_et_bins[et][p] -> Fill( dca2d  );
						}
						if(trackID<0) 
						{
							KF_Failed_Counts[et][p] += 1; // Kalman Filter Failed Eff
						}
						
						CrisTotalCounts[et][p]+=1; // Counting the number of events in each et-p bin
						
					}	
				}
				for(int pT = 0 ; pT < size_mom_bin-1 ; pT++){
					if(pT_truth > pT_bin[pT] && pT_truth <= pT_bin[pT+1]){
						h1_dca2d_pT_et_bins[et][pT] -> Fill( dca2d  );
					}				
				}
			}
		}
	}
	// -------------------------------------------------------------
	// Doing fits
	TCanvas ** c_fits_p  = new TCanvas*[size_eta_bin-1];
	TCanvas ** c_fits_th = new TCanvas*[size_eta_bin-1];
	TCanvas ** c_fits_ph = new TCanvas*[size_eta_bin-1];
	TCanvas ** c_fits_dca2d = new TCanvas*[size_eta_bin-1];
	TCanvas ** c_fits_dca2d_pT = new TCanvas*[size_eta_bin-1];
	
	ofstream outfile;
	ofstream outfile1;
	
	if(isDefault) outfile.open("default_params.csv");
	else outfile.open("params.csv");
	outfile<<"etarange,prange,dp_p_p,error_dp_p_p,dth_th_p,error_dth_th_p,dph_ph_p,error_dph_ph_p,dca2d,error_dca2d,"
		"dca2d_v_pT,error_dca2d_v_pT,InEfficiency,error_InEfficiency,KF_InEfficiency,error_KF_InEfficiency,"
		"GlobalKFInEff,error_GlobalKFInEff,EventsInBin,A1,A2,A1A2,Integral,PseudoGlobalKFInEff,"
		"error_PseudoGlobalKFInEff,OverFlowUnderFlowFlag\n";

outfile1.open("Chi2NDF_DoubleGaus.csv");

outfile1<<"etarange,prange,Chi2_dpp,NDF_dpp,Chi2_dth,NDF_dth,Chi2_dph,NDF_dph,Chi2_dca2d,NDF_dca2d\n";

	for(int et = 0 ; et < size_eta_bin-1 ; et++){
		c_fits_p [et] = new TCanvas(Form("c_fits_p_%i" ,et),Form("dp/p  , %.1f<eta<%.1f",eta_bin[et],eta_bin[et+1]),1000,800);	c_fits_p [et] -> Divide(3,2);
		c_fits_th[et] = new TCanvas(Form("c_fits_th_%i",et),Form("dtheta, %.1f<eta<%.1f",eta_bin[et],eta_bin[et+1]),1000,800);	c_fits_th[et] -> Divide(3,2);
		c_fits_ph[et] = new TCanvas(Form("c_fits_ph_%i",et),Form("dphi  , %.1f<eta<%.1f",eta_bin[et],eta_bin[et+1]),1000,800);	c_fits_ph[et] -> Divide(3,2);
		c_fits_dca2d[et] = new TCanvas(Form("c_fits_dca2d_%i",et),Form("dca2d  , %.1f<eta<%.1f",eta_bin[et],eta_bin[et+1]),1000,800);	c_fits_dca2d[et] -> Divide(3,2);
		c_fits_dca2d_pT[et] = new TCanvas(Form("c_fits_dca2d_pT_%i",et),Form("dca2d  , %.1f<eta<%.1f",eta_bin[et],eta_bin[et+1]),1000,800);	c_fits_dca2d_pT[et] -> Divide(3,2);

		for(int p = 0 ; p < size_mom_bin-1 ; p++){
			
	//******************** Momentum resolutions *************************************** //

			c_fits_p [et] -> cd(p+1);
			//h1_dpp_p_et_bins[et][p] -> Draw();	h1_dpp_p_et_bins[et][p] -> Fit(Form("f_gaus_dpp_%i_%i",et,p),"RQ");
			//width_dpp[et][p] = f_gaus_dpp[et][p] -> GetParameter(2);
			//error_dpp[et][p] = (f_gaus_dpp[et][p] -> GetParError(2))*(f_gaus_dpp[et][p] -> GetChisquare())/(f_gaus_dpp[et][p] -> GetNDF());
			h1_dpp_p_et_bins[et][p]->GetXaxis()->SetRangeUser(-0.5, 0.5);
			double Meandpp = h1_dpp_p_et_bins[et][p]->GetMean();
			double Sigmadpp = h1_dpp_p_et_bins[et][p]->GetRMS();
			double Maxdpp = h1_dpp_p_et_bins[et][p]->GetMaximum();
			h1_dpp_p_et_bins[et][p]->GetXaxis()->SetRangeUser(-1, 1);

			f_gaus_dpp[et][p] -> SetParameters(Maxdpp*0.8, Meandpp, Sigmadpp*0.5, Maxdpp*0.2,Sigmadpp*2);
			
			f_gaus_dpp[et][p] -> SetParLimits(0, 0, Maxdpp); // Make sure there the amplitude is not negative

			f_gaus_dpp[et][p] -> SetParLimits(1, Meandpp - 0.2*Sigmadpp, 0.2*Meandpp + Sigmadpp); // Mean range		

			f_gaus_dpp[et][p] -> SetParLimits(2, 0.1*Sigmadpp, 50*Sigmadpp); // Make sure there the sigma1 is not negative

			f_gaus_dpp[et][p] -> SetParLimits(3, 0, Maxdpp); // Make sure the amplitude of the second gaussian is also not negative

			f_gaus_dpp[et][p] -> SetParLimits(4, 0.1*Sigmadpp, 50*Sigmadpp); // Make sure there the sigma2 is not negative


			h1_dpp_p_et_bins[et][p] -> Draw();	h1_dpp_p_et_bins[et][p] -> Fit(Form("f_gaus_dpp_%i_%i",et,p),"Q", "", -5*Sigmadpp, 5*Sigmadpp);

			// Calculation of momentum resolution

			TF1 *ppgaus1 = new TF1(Form("ppgaus1_dpp_%i_%i", et, p), "gaus", -1, 1);
			TF1 *ppgaus2 = new TF1(Form("ppgaus2_dpp_%i_%i", et, p), "gaus", -1, 1);

			
			ppgaus1->SetParameters(f_gaus_dpp[et][p]->GetParameter(0), f_gaus_dpp[et][p]->GetParameter(1), f_gaus_dpp[et][p]->GetParameter(2));
			ppgaus2->SetParameters(f_gaus_dpp[et][p]->GetParameter(3), f_gaus_dpp[et][p]->GetParameter(1), f_gaus_dpp[et][p]->GetParameter(4));

			
			double A1pp = ppgaus1->Integral(-0.3, 0.3); // Calculating the Area of the gaussian func 1
			double s1pp = f_gaus_dpp[et][p] -> GetParameter(2);
			double error_s1pp = (f_gaus_dpp[et][p] -> GetParError(2))*(f_gaus_dpp[et][p] -> GetChisquare())/(f_gaus_dpp[et][p] -> GetNDF());

			double A2pp = ppgaus2->Integral(-0.3, 0.3); // Calculating the Area of the gaussian func 2
			double s2pp = f_gaus_dpp[et][p] -> GetParameter(4);	
			double error_s2pp = (f_gaus_dpp[et][p] -> GetParError(4))*(f_gaus_dpp[et][p] -> GetChisquare())/(f_gaus_dpp[et][p] -> GetNDF());

			// Calculating the Ratio of (A1+A2)/Integral() a measure of how good the reconstructions are 

			A1_inf[et][p] = ppgaus1->Integral(-1, 1); // The range of -1 to 1 for the histogram
			A2_inf[et][p] = ppgaus2->Integral(-1, 1); 
			Entry_inf[et][p] = h1_dpp_p_et_bins[et][p]->Integral("width"); // This does not include over-flow and under-flow bins

			
			// Checking for the overflow and underflow and poor reconstructions

			int NBins = h1_dpp_p_et_bins[et][p]->GetNbinsX();
			double entries_outside = h1_dpp_p_et_bins[et][p]->GetBinContent(-1) + h1_dpp_p_et_bins[et][p]->Integral(0, (int)NBins*1/4) + h1_dpp_p_et_bins[et][p]->Integral((int)NBins*3/4, NBins) + h1_dpp_p_et_bins[et][p]->GetBinContent(NBins + 1);

			OverFlowUnderFlowFlag[et][p] = (entries_outside/h1_dpp_p_et_bins[et][p]->GetEntries()>=0.1);	// compare with all Entries including the overflow and underflow	


			Chi2_dpp[et][p] = f_gaus_dpp[et][p] -> GetChisquare();
			NDF_dpp[et][p] = f_gaus_dpp[et][p] -> GetNDF();

			if(ppgaus2->Eval(0) < 0.02*ppgaus1->Eval(0) || s2pp > 10*Sigmadpp || ppgaus1->Eval(0) < 0.02*ppgaus2->Eval(0) || s1pp > 10*Sigmadpp || (ppgaus1->Eval(0) + ppgaus2->Eval(0)) < 0.6*Maxdpp) 
			{ 
				TF1 *ppgaus_refit = new TF1("gaus", "gaus", -1, 1);
				h1_dpp_p_et_bins[et][p] -> Fit(ppgaus_refit,"Q", "", -5*Sigmadpp, 5*Sigmadpp);	
				if(Chi2_dpp[et][p] < ppgaus_refit -> GetChisquare()){				
				A1pp = ppgaus_refit->Integral(-0.3, 0.3);
				s1pp = ppgaus_refit -> GetParameter(2);
				error_s1pp = (ppgaus_refit -> GetParError(2))*(ppgaus_refit -> GetChisquare())/(ppgaus_refit -> GetNDF());
				A2pp = 0.; s2pp = 0;
				f_gaus_dpp[et][p]->SetParameters(ppgaus_refit -> GetParameter(0), ppgaus_refit -> GetParameter(1), ppgaus_refit -> GetParameter(2), 0, 0);

				A1_inf[et][p] = ppgaus_refit->Integral(-1, 1);
				A2_inf[et][p] = 0.;
				Entry_inf[et][p] = h1_dpp_p_et_bins[et][p]->Integral("width");

				Chi2_dpp[et][p] = ppgaus_refit -> GetChisquare();
				NDF_dpp[et][p] = ppgaus_refit -> GetNDF();
				}

			}

			

			width_dpp[et][p] = (A1pp*s1pp + A2pp*s2pp)/(A1pp + A2pp);
			error_dpp[et][p] = (A1pp*error_s1pp + A2pp*error_s2pp)/(A1pp + A2pp);

			f_gaus1_dpp[et][p]->SetParameters(f_gaus_dpp[et][p]->GetParameter(0), f_gaus_dpp[et][p]->GetParameter(1), f_gaus_dpp[et][p]->GetParameter(2));
			f_gaus2_dpp[et][p]->SetParameters(f_gaus_dpp[et][p]->GetParameter(3), f_gaus_dpp[et][p]->GetParameter(1), f_gaus_dpp[et][p]->GetParameter(4));



			sigma1_dpp[et][p] = s1pp;
			sigma2_dpp[et][p] = s2pp;

			A1_dpp[et][p] = A1pp;
			A2_dpp[et][p] = A2pp;




// ********************* end of Momentum Resolution ****************************************************//


			// Reconstruction Efficiency 

			mean_dpp[et][p] = f_gaus_dpp[et][p] -> GetParameter(1);
			double dpp_leftend = mean_dpp[et][p] - 5*width_dpp[et][p];
			double LeftXBin = h1_dpp_p_et_bins[et][p]->GetXaxis()->FindBin(dpp_leftend);
			double dpp_rightend = mean_dpp[et][p] + 5*width_dpp[et][p];
			double RightXBin = h1_dpp_p_et_bins[et][p]->GetXaxis()->FindBin(dpp_rightend);
			Yield[et][p] = h1_dpp_p_et_bins[et][p]->Integral(LeftXBin, RightXBin);
			Eff[et][p] = Yield[et][p]/TotalCounts[et][p];


// ******************************** Theta resolution ****************************************************** //


			c_fits_th[et] -> cd(p+1);
			//h1_dth_p_et_bins[et][p] -> Draw();	h1_dth_p_et_bins[et][p] -> Fit(Form("f_gaus_dth_%i_%i",et,p),"RQ");
			//width_dth[et][p] = f_gaus_dth[et][p] -> GetParameter(2);
			//error_dth[et][p] = (f_gaus_dth[et][p] -> GetParError(2))*(f_gaus_dth[et][p] -> GetChisquare())/(f_gaus_dth[et][p] -> GetNDF());

			double Maxdth = h1_dth_p_et_bins[et][p]->GetMaximum();
			double Meandth = h1_dth_p_et_bins[et][p]->GetMean();
			double Sigmadth = h1_dth_p_et_bins[et][p]->GetRMS();


			f_gaus_dth[et][p] -> SetParameters(Maxdth*0.8, 0.0, Sigmadth*0.5, Maxdth*0.2, Sigmadth*2); // setting the mean to 0.0
			
			f_gaus_dth[et][p] -> SetParLimits(0, 0, Maxdth); // Make sure there the amplitude is not negative		
			
			f_gaus_dth[et][p] -> SetParLimits(1, -0.01, 0.01); // Mean fit range

			f_gaus_dth[et][p] -> SetParLimits(2, 0, 20*Sigmadth); // Make sure there the sigma1 is not negative

			f_gaus_dth[et][p] -> SetParLimits(3, 0, Maxdth); // Make sure the amplitude of the second gaussian is also not negative

			f_gaus_dth[et][p] -> SetParLimits(4, 0, 20*Sigmadth); // Make sure there the sigma2 is not negative

			h1_dth_p_et_bins[et][p] -> Draw();	h1_dth_p_et_bins[et][p] -> Fit(Form("f_gaus_dth_%i_%i",et,p),"RQ");

			// Calculation of theta resolution

			TF1 *thgaus1 = new TF1(Form("thgaus1_dth_%i_%i", et, p), "gaus", -0.1, 0.1);
			TF1 *thgaus2 = new TF1(Form("thgaus2_dth_%i_%i", et, p), "gaus", -0.1, 0.1);
			thgaus1->SetParameters(f_gaus_dth[et][p]->GetParameter(0), f_gaus_dth[et][p]->GetParameter(1), f_gaus_dth[et][p]->GetParameter(2));
			thgaus2->SetParameters(f_gaus_dth[et][p]->GetParameter(3), f_gaus_dth[et][p]->GetParameter(1), f_gaus_dth[et][p]->GetParameter(4));

			
			double A1th = thgaus1->Integral(-0.1, 0.1); // Calculating the Area of the gaussian func 1
			double s1th = f_gaus_dth[et][p] -> GetParameter(2);
			double error_s1th = (f_gaus_dth[et][p] -> GetParError(2))*(f_gaus_dth[et][p] -> GetChisquare())/(f_gaus_dth[et][p] -> GetNDF());

			double A2th = thgaus2->Integral(-0.1, 0.1); // Calculating the Area of the gaussian func 1
			double s2th = f_gaus_dth[et][p] -> GetParameter(4);	
			double error_s2th = (f_gaus_dth[et][p] -> GetParError(4))*(f_gaus_dth[et][p] -> GetChisquare())/(f_gaus_dth[et][p] -> GetNDF());

			//if(thgaus2->Eval(0) < 0.02*thgaus1->Eval(0) || s2th > h1_dth_p_et_bins[et][p]->GetXaxis()->GetXmax()) { A2th = 0.; s2th = 0; }
			//if(thgaus1->Eval(0) < 0.02*thgaus2->Eval(0) || s2th > h1_dth_p_et_bins[et][p]->GetXaxis()->GetXmax()) { A1th = 0.; s1th = 0; } 
			Chi2_dth[et][p] = f_gaus_dth[et][p] -> GetChisquare();
			NDF_dth[et][p] = f_gaus_dth[et][p] -> GetNDF();

			if(thgaus2->Eval(0) < 0.02*thgaus1->Eval(0) || s2th > h1_dth_p_et_bins[et][p]->GetXaxis()->GetXmax() || thgaus1->Eval(0) < 0.02*thgaus2->Eval(0) || s2th > h1_dth_p_et_bins[et][p]->GetXaxis()->GetXmax()) 
			{ 
				TF1 *thgaus_refit = new TF1("thgaus", "gaus", -0.1, 0.1);
				h1_dth_p_et_bins[et][p] -> Fit(thgaus_refit,"RQ");	
				A1th = thgaus_refit->Integral(-0.1, 0.1);
				s1th = thgaus_refit -> GetParameter(2);
				error_s1th = (thgaus_refit -> GetParError(2))*(thgaus_refit -> GetChisquare())/(thgaus_refit -> GetNDF());
				A2th = 0.; s2th = 0;
				f_gaus_dth[et][p]->SetParameters(thgaus_refit -> GetParameter(0), thgaus_refit -> GetParameter(1), thgaus_refit -> GetParameter(2), 0, 0);

				Chi2_dth[et][p] = thgaus_refit -> GetChisquare();
				NDF_dth[et][p] = thgaus_refit -> GetNDF();

			}		
			
			width_dth[et][p] = (A1th*s1th + A2th*s2th)/(A1th + A2th);
			error_dth[et][p] = (A1th*error_s1th + A2th*error_s2th)/(A1th + A2th);

			f_gaus1_dth[et][p]->SetParameters(f_gaus_dth[et][p]->GetParameter(0), f_gaus_dth[et][p]->GetParameter(1), f_gaus_dth[et][p]->GetParameter(2));
			f_gaus2_dth[et][p]->SetParameters(f_gaus_dth[et][p]->GetParameter(3), f_gaus_dth[et][p]->GetParameter(1), f_gaus_dth[et][p]->GetParameter(4));



			sigma1_dth[et][p] = s1th;
			sigma2_dth[et][p] = s2th;

			A1_dth[et][p] = A1th;
			A2_dth[et][p] = A2th;


// ******************************** end of theta Resolution **************************************//



// ********************************** Phi resolution ************************************************//

			c_fits_ph[et] -> cd(p+1);
			//h1_dph_p_et_bins[et][p] -> Draw();	h1_dph_p_et_bins[et][p] -> Fit(Form("f_gaus_dph_%i_%i",et,p),"RQ");

			double Meandph = h1_dph_p_et_bins[et][p]->GetMean();
			double Maxdph = h1_dph_p_et_bins[et][p]->GetMaximum();
			double Sigmadph = h1_dph_p_et_bins[et][p]->GetRMS();
			
			f_gaus_dph[et][p] -> SetParameters(Maxdph*0.8, 0.0, Sigmadph*0.5, Maxdph*0.2, Sigmadph*2);
			
			f_gaus_dph[et][p] -> SetParLimits(0, 0, Maxdph); // Make sure there the amplitude is not negative
			f_gaus_dph[et][p] -> SetParLimits(1, -0.01, 0.01);

			f_gaus_dph[et][p] -> SetParLimits(2, 0, 50*Sigmadph); // Make sure there the sigma1 is not negative

			f_gaus_dph[et][p] -> SetParLimits(3, 0, Maxdph); // Make sure the amplitude of the second gaussian is also not negative

			f_gaus_dph[et][p] -> SetParLimits(4, 0, 50*Sigmadph); // Make sure there the sigma2 is not negative

			h1_dph_p_et_bins[et][p] -> Draw();	h1_dph_p_et_bins[et][p] -> Fit(Form("f_gaus_dph_%i_%i",et,p),"Q", "", -0.05, 0.05);

			// Calculation of phi resolution

			TF1 *phgaus1 = new TF1(Form("phgaus1_dph_%i_%i", et, p), "gaus", -0.1, 0.1);
			TF1 *phgaus2 = new TF1(Form("phgaus2_dph_%i_%i", et, p), "gaus", -0.1, 0.1);
			phgaus1->SetParameters(f_gaus_dph[et][p]->GetParameter(0), f_gaus_dph[et][p]->GetParameter(1), f_gaus_dph[et][p]->GetParameter(2));
			phgaus2->SetParameters(f_gaus_dph[et][p]->GetParameter(3), f_gaus_dph[et][p]->GetParameter(1), f_gaus_dph[et][p]->GetParameter(4));

			
			double A1ph = phgaus1->Integral(-0.1, 0.1); // Calculating the Area of the gaussian func 1
			double s1ph = f_gaus_dph[et][p] -> GetParameter(2);
			double error_s1ph = (f_gaus_dph[et][p] -> GetParError(2))*(f_gaus_dph[et][p] -> GetChisquare())/(f_gaus_dph[et][p] -> GetNDF());

			double A2ph = phgaus2->Integral(-0.1, 0.1); // Calculating the Area of the gaussian func 1
			double s2ph = f_gaus_dph[et][p] -> GetParameter(4);	
			double error_s2ph = (f_gaus_dph[et][p] -> GetParError(4))*(f_gaus_dph[et][p] -> GetChisquare())/(f_gaus_dph[et][p] -> GetNDF());

			Chi2_dph[et][p] = f_gaus_dph[et][p] -> GetChisquare();
			NDF_dph[et][p] = f_gaus_dph[et][p] -> GetNDF();

			/*if(phgaus2->Eval(0) < 0.02*phgaus1->Eval(0) || s2ph > h1_dph_p_et_bins[et][p]->GetXaxis()->GetXmax() || phgaus1->Eval(0) < 0.02*phgaus2->Eval(0) || s2ph > h1_dph_p_et_bins[et][p]->GetXaxis()->GetXmax() || (phgaus2->Eval(0) + phgaus1->Eval(0)) < 0.6*Maxdph) 
			{ 
				TF1 *phgaus_refit = new TF1("phgaus", "gaus", -0.1, 0.1);
				h1_dph_p_et_bins[et][p] -> Fit(phgaus_refit,"Q", "", -0.05, 0.05);	
				A1ph = phgaus_refit->Integral(-0.1, 0.1);
				s1ph = phgaus_refit -> GetParameter(2);
				error_s1ph = (phgaus_refit -> GetParError(2))*(phgaus_refit -> GetChisquare())/(phgaus_refit -> GetNDF());
				A2ph = 0.; s2ph = 0;
				f_gaus_dph[et][p]->SetParameters(phgaus_refit -> GetParameter(0), phgaus_refit -> GetParameter(1), phgaus_refit -> GetParameter(2), 0, 0);

				Chi2_dph[et][p] = phgaus_refit -> GetChisquare();
				NDF_dph[et][p] = phgaus_refit -> GetNDF();

			}*/		
			
			width_dph[et][p] = (A1ph*s1ph + A2ph*s2ph)/(A1ph + A2ph);
			error_dph[et][p] = (A1ph*error_s1ph + A2ph*error_s2ph)/(A1ph + A2ph);

			f_gaus1_dph[et][p]->SetParameters(f_gaus_dph[et][p]->GetParameter(0), f_gaus_dph[et][p]->GetParameter(1), f_gaus_dph[et][p]->GetParameter(2));
			f_gaus2_dph[et][p]->SetParameters(f_gaus_dph[et][p]->GetParameter(3), f_gaus_dph[et][p]->GetParameter(1), f_gaus_dph[et][p]->GetParameter(4));



			sigma1_dph[et][p] = s1ph;
			sigma2_dph[et][p] = s2ph;

			A1_dph[et][p] = A1ph;
			A2_dph[et][p] = A2ph;

// ****************************** end of Phi Resolution **************************************** //



// ***************************** dca2d resolution ********************************************* //

			c_fits_dca2d[et] -> cd(p+1);
			// Set +-5 RMS() of distribution for fitting. This is to set the Fitting dynamic 
			double leftend = h1_dca2d_p_et_bins[et][p]->GetMean() - 3*h1_dca2d_p_et_bins[et][p]->GetRMS();
			double rightend = h1_dca2d_p_et_bins[et][p]->GetMean() + 3*h1_dca2d_p_et_bins[et][p]->GetRMS();

			// Setting parameters, This is done by looking into each histograms and finding the combinations of initial parameters that made the ROOT::Fit converge

			f_gaus_dca2d[et][p] -> SetParameters(h1_dca2d_p_et_bins[et][p]->GetMaximum()*0.8, h1_dca2d_p_et_bins[et][p]->GetMean(), h1_dca2d_p_et_bins[et][p]->GetRMS()*0.5, h1_dca2d_p_et_bins[et][p]->GetMaximum()*0.2, h1_dca2d_p_et_bins[et][p]->GetRMS()*3);
			
			f_gaus_dca2d[et][p] -> SetParLimits(0, 0, h1_dca2d_p_et_bins[et][p]->GetMaximum()); // Make sure there the amplitude is not negative
			f_gaus_dca2d[et][p] -> SetParLimits(1, -0.01, 0.01);

			f_gaus_dca2d[et][p] -> SetParLimits(2, 0, 50*h1_dca2d_p_et_bins[et][p]->GetRMS()); // Make sure there the sigma1 is not negative

			f_gaus_dca2d[et][p] -> SetParLimits(3, 0, h1_dca2d_p_et_bins[et][p]->GetMaximum()); // Make sure the amplitude of the second gaussian is also not negative

			f_gaus_dca2d[et][p] -> SetParLimits(4, 0, 50*h1_dca2d_p_et_bins[et][p]->GetRMS()); // Make sure there the sigma2 is not negative

			h1_dca2d_p_et_bins[et][p] -> Draw();	h1_dca2d_p_et_bins[et][p] -> Fit(Form("f_gaus_dca2d_%i_%i",et,p),"Q","",leftend, rightend);
			
			// Calculation of dca2d resolution

			TF1 *gaus1 = new TF1(Form("gaus1_dca2d_%i_%i", et, p), "gaus", leftend, rightend);
			TF1 *gaus2 = new TF1(Form("gaus2_dca2d_%i_%i", et, p), "gaus", leftend, rightend);
			gaus1->SetParameters(f_gaus_dca2d[et][p]->GetParameter(0), f_gaus_dca2d[et][p]->GetParameter(1), f_gaus_dca2d[et][p]->GetParameter(2));
			gaus2->SetParameters(f_gaus_dca2d[et][p]->GetParameter(3), f_gaus_dca2d[et][p]->GetParameter(1), f_gaus_dca2d[et][p]->GetParameter(4));

			
			double A1 = gaus1->Integral(leftend, rightend); // Calculating the Area of the gaussian func 1
			double s1 = f_gaus_dca2d[et][p] -> GetParameter(2);
			double error_s1 = (f_gaus_dca2d[et][p] -> GetParError(2))*(f_gaus_dca2d[et][p] -> GetChisquare())/(f_gaus_dca2d[et][p] -> GetNDF());
			

			double A2 = gaus2->Integral(leftend, rightend); // Calculating the Area of the gaussian func 1
			double s2 = f_gaus_dca2d[et][p] -> GetParameter(4);	
			double error_s2 = (f_gaus_dca2d[et][p] -> GetParError(4))*(f_gaus_dca2d[et][p] -> GetChisquare())/(f_gaus_dca2d[et][p] -> GetNDF());
			//cout<< error_s1 << " " << A1 << " " << error_s2 << " " << A2 << endl;
			//if(gaus2->Eval(0) < 0.02*gaus1->Eval(0)|| s2 > h1_dca2d_p_et_bins[et][p]->GetXaxis()->GetXmax()) { A2 = 0.; s2 = 0; }
			//if(gaus1->Eval(0) < 0.02*gaus2->Eval(0)|| s2 > h1_dca2d_p_et_bins[et][p]->GetXaxis()->GetXmax()) { A1 = 0.; s1 = 0; } 

			Chi2_dca2d[et][p] = f_gaus_dca2d[et][p] -> GetChisquare();
			NDF_dca2d[et][p] = f_gaus_dca2d[et][p] -> GetNDF();

			if(gaus2->Eval(0) < 0.02*gaus1->Eval(0) || s2 > h1_dca2d_p_et_bins[et][p]->GetXaxis()->GetXmax() || gaus1->Eval(0) < 0.02*gaus2->Eval(0) || s2 > h1_dca2d_p_et_bins[et][p]->GetXaxis()->GetXmax()) 
			{ 
				TF1 *dca2dgaus_refit = new TF1("dca2dgaus", "gaus", -0.3, 0.3);
				h1_dca2d_p_et_bins[et][p] -> Fit(dca2dgaus_refit,"Q", "", leftend, rightend);	
				A1 = dca2dgaus_refit->Integral(leftend, rightend);
				s1 = dca2dgaus_refit -> GetParameter(2);
				error_s1 = (dca2dgaus_refit -> GetParError(2))*(dca2dgaus_refit -> GetChisquare())/(dca2dgaus_refit -> GetNDF());
				A2 = 0.; s2 = 0;
				f_gaus_dca2d[et][p]->SetParameters(dca2dgaus_refit -> GetParameter(0), dca2dgaus_refit -> GetParameter(1), dca2dgaus_refit -> GetParameter(2), 0, 0);

				Chi2_dca2d[et][p] = dca2dgaus_refit -> GetChisquare();
				NDF_dca2d[et][p] = dca2dgaus_refit -> GetNDF();

			}


			width_dca2d[et][p] = (A1*s1 + A2*s2)/(A1 + A2);
			error_dca2d[et][p] = (A1*error_s1 + A2*error_s2)/(A1 + A2);	
			//width_dca2d[et][p] = h1_dca2d_p_et_bins[et][p] -> GetRMS(); //currently getting the RMS as the sigma
			//error_dca2d[et][p] = h1_dca2d_p_et_bins[et][p] -> GetRMSError(); 

			f_gaus1_dca2d[et][p]->SetParameters(f_gaus_dca2d[et][p]->GetParameter(0), f_gaus_dca2d[et][p]->GetParameter(1), f_gaus_dca2d[et][p]->GetParameter(2));
			f_gaus2_dca2d[et][p]->SetParameters(f_gaus_dca2d[et][p]->GetParameter(3), f_gaus_dca2d[et][p]->GetParameter(1), f_gaus_dca2d[et][p]->GetParameter(4));





			sigma1_dca2d[et][p] = s1;
			sigma2_dca2d[et][p] = s2;

			A1_dca2d[et][p] = A1;
			A2_dca2d[et][p] = A2;

// ************************* end of dca2d Resolution ********************************** //


// ********************** dca2d resolution vs pT *************************************************** //

			c_fits_dca2d_pT[et] -> cd(p+1);
			// Set +-5 RMS() of distribution for fitting. This is to set the Fitting dynamic 
			double leftend_pT = h1_dca2d_pT_et_bins[et][p]->GetMean() - 5*h1_dca2d_pT_et_bins[et][p]->GetRMS();
			double rightend_pT = h1_dca2d_pT_et_bins[et][p]->GetMean() + 5*h1_dca2d_pT_et_bins[et][p]->GetRMS();
			
			if(h1_dca2d_pT_et_bins[et][p] ->GetEntries() > 20)
			{
				h1_dca2d_pT_et_bins[et][p] -> Draw();	h1_dca2d_pT_et_bins[et][p] -> Fit(Form("f_gaus_dca2d_pT_%i_%i",et,p),"Q","",leftend_pT, rightend_pT);
			width_dca2d_pT[et][p] = f_gaus_dca2d_pT[et][p] -> GetParameter(2);
			error_dca2d_pT[et][p] = (f_gaus_dca2d_pT[et][p] -> GetParError(2))*(f_gaus_dca2d_pT[et][p] -> GetChisquare())/(f_gaus_dca2d_pT[et][p] -> GetNDF());			
			}
			else
			{
				width_dca2d_pT[et][p] = 0.0;
				error_dca2d_pT[et][p] = 0.0;
			}
			//width_dca2d[et][p] = h1_dca2d_p_et_bins[et][p] -> GetRMS(); //currently getting the RMS as the sigma
			//error_dca2d[et][p] = h1_dca2d_p_et_bins[et][p] -> GetRMSError(); 
			
// *************************** end of dca2d vs pT resolution ***************************************** //

			
			// Pseudo InEfficiency Calculation
			float PseudoInEff = 1 - CrisEff[et][p]/CrisTotalCounts[et][p]; // This is better since directly can use it in obj_fun.py
			error_CrisEff[et][p] = sqrt(PseudoInEff*(1-PseudoInEff)/CrisTotalCounts[et][p]);

			// error KF_InEfficiency Calculation
			float KF_InEff = KF_Failed_Counts[et][p]/CrisTotalCounts[et][p];
			error_KF_Eff[et][p] = sqrt(KF_InEff*(1-KF_InEff)/CrisTotalCounts[et][p]);

			// Global KF_InEfficiency Calculation
			float Global_KF_InEff = (float)T->GetEntries("trackID<0")/nEntries;
			float error_Global_KF_InEff = sqrt(Global_KF_InEff*(Global_KF_InEff)/nEntries);

			//float PseudoWeightKFGlobalIneff = pow((Chi2_dpp[et][p]/NDF_dpp[et][p])*(Entry_inf[et][p]/(A1_inf[et][p] + A2_inf[et][p]) - 1), 2);
			float PseudoWeightKFGlobalIneff = pow((Chi2_dpp[et][p]/NDF_dpp[et][p])*(Entry_inf[et][p]/(A1_inf[et][p] + A2_inf[et][p]) - 1), 2);
			PseudoGlobalKFInEff[et][p] = Global_KF_InEff*PseudoWeightKFGlobalIneff;
			error_PseudoGlobalKFInEff[et][p] = error_Global_KF_InEff*PseudoWeightKFGlobalIneff;

			outfile
			<<eta_bin[et]<< " - "<<eta_bin[et+1]<<"," 	// etarange
			<<mom_bin[p]<<" - "<<mom_bin[p+1]<<"," 		// prange
			<<width_dpp[et][p]*100<< "," 			// dp_p_p
			<< error_dpp[et][p]*100<<","			// error_dp_p_p
			<<width_dth[et][p]*1000<< ","			// dth_th_p
			<< error_dth[et][p]*1000<<"," 			// error_dth_th_p
			<< width_dph[et][p]*1000<< ","			// dph_ph_p
			<< error_dph[et][p]*1000<<","			// error_dph_ph_p
			<<width_dca2d[et][p]*10000<<","			// dca2d [ums]
			<<error_dca2d[et][p]*10000<<","			// error_dca2d [ums]
			<<width_dca2d_pT[et][p]*10000<<","		// dca2d vs pT [ums]
			<<error_dca2d_pT[et][p]*10000<<","		// error_dca2d vs pT [ums]
			<<PseudoInEff<<","				// InEfficiency (As discussed on June 15 2021)
			<<error_CrisEff[et][p]<<","			// error Efficiency
			<<KF_Failed_Counts[et][p]/CrisTotalCounts[et][p]<<","		// KF Failed In Efficiency; 
			<<error_KF_Eff[et][p]<<","			// error KF In Efficiency
			<<Global_KF_InEff<<","				// Global KF InEff
			<<error_Global_KF_InEff<<","			// error Global KF InEff
			<<CrisTotalCounts[et][p]<<","			// Events in Each Bin
			<<A1_inf[et][p]<<","				// A1 for dpp until -1, 1
			<<A2_inf[et][p]<<","				// A2 for dpp until -1, 1
			<<A1_inf[et][p] + A2_inf[et][p]<<","		// A1 + A2
			<<Entry_inf[et][p]<<","				// TotalIntegral dpp
			<<PseudoGlobalKFInEff[et][p]*10E6<<","		// PseudoGlobalKFInEff KFIneff*((Chi2/NDF)*(Integral/(A1 + A2) -1)^2)
			<<error_PseudoGlobalKFInEff[et][p]*10E6<<","	// error in PseudoGlobalKFInEff
			<<OverFlowUnderFlowFlag[et][p]<<"\n";		// Overflow and Underflow Flag All has to be False
			// ----


			outfile1
			<<eta_bin[et]<< " - "<<eta_bin[et+1]<<"," 	// etarange
			<<mom_bin[p]<<" - "<<mom_bin[p+1]<<"," 		// prange
			<<Chi2_dpp[et][p]<<","
			<<NDF_dpp[et][p]<<","
			<<Chi2_dth[et][p]<<","
			<<NDF_dth[et][p]<<","
			<<Chi2_dph[et][p]<<","
			<<NDF_dph[et][p]<<","
			<<Chi2_dca2d[et][p]<<","
			<<NDF_dca2d[et][p]<<"\n";

			h1_dpp_v_p_et_bins[et] -> SetBinContent(p +1,width_dpp[et][p]*100. );
			h1_dth_v_p_et_bins[et] -> SetBinContent(p +1,width_dth[et][p]*1000.);
			h1_dph_v_p_et_bins[et] -> SetBinContent(p +1,width_dph[et][p]*1000.);
			h1_dca2d_v_p_et_bins[et] -> SetBinContent(p +1,width_dca2d[et][p]*10000);
			h1_dca2d_v_pT_et_bins[et] -> SetBinContent(p +1,width_dca2d[et][p]*10000);
			
			h1_dpp_v_et_p_bins[ p] -> SetBinContent(et+1,width_dpp[et][p]*100. );
			h1_dth_v_et_p_bins[ p] -> SetBinContent(et+1,width_dth[et][p]*1000.);
			h1_dph_v_et_p_bins[ p] -> SetBinContent(et+1,width_dph[et][p]*1000.);
			h1_dca2d_v_et_p_bins[ p] -> SetBinContent(et+1,width_dca2d[et][p]*10000);
			h1_dca2d_v_et_pT_bins[ p] -> SetBinContent(et+1,width_dca2d[et][p]*10000);

			h1_dpp_v_p_et_bins[et] -> SetBinError  (p +1,error_dpp[et][p]*100. );
			h1_dth_v_p_et_bins[et] -> SetBinError  (p +1,error_dth[et][p]*1000.);
			h1_dph_v_p_et_bins[et] -> SetBinError  (p +1,error_dph[et][p]*1000.);
			h1_dca2d_v_p_et_bins[et] -> SetBinError  (p +1,error_dca2d[et][p]*10000);
			h1_dca2d_v_pT_et_bins[et] -> SetBinError  (p +1,error_dca2d[et][p]*10000);

			h1_dpp_v_et_p_bins[ p] -> SetBinError  (et+1,error_dpp[et][p]*100. );
			h1_dth_v_et_p_bins[ p] -> SetBinError  (et+1,error_dth[et][p]*1000.);
			h1_dph_v_et_p_bins[ p] -> SetBinError  (et+1,error_dph[et][p]*1000.);
			h1_dca2d_v_et_p_bins[ p] -> SetBinError  (et+1,error_dca2d[et][p]*10000);
			h1_dca2d_v_et_pT_bins[ p] -> SetBinError  (et+1,error_dca2d[et][p]*10000);
		}
	}
outfile.close();
outfile1.close();

	// -------------------------------------------------------------
	// Updating table with width values
	ofstream updated_tab;
	if(update_tab){
		updated_tab.open("tables/"+tab_name);
		for(int et = 0 ; et < size_eta_bin-1 ; et++){
			for(int p = 0 ; p < size_mom_bin-1 ; p++){
				updated_tab << width_dpp[et][p];
				if(p == size_mom_bin-2) updated_tab << "\n";
				else updated_tab << "\t";
			}
		}
		updated_tab << "\n";
		for(int et = 0 ; et < size_eta_bin-1 ; et++){
			for(int p = 0 ; p < size_mom_bin-1 ; p++){
				updated_tab << width_dth[et][p];
				if(p == size_mom_bin-2) updated_tab << "\n";
				else updated_tab << "\t";
			}
		}
		updated_tab << "\n";
		for(int et = 0 ; et < size_eta_bin-1 ; et++){
			for(int p = 0 ; p < size_mom_bin-1 ; p++){
				updated_tab << width_dph[et][p];
				if(p == size_mom_bin-2) updated_tab << "\n";
				else updated_tab << "\t";
			}
		}
		updated_tab.close();
	}

	// -------------------------------------------------------------
	// Plotting histograms
	TCanvas * c1 = new TCanvas("c1","c1",1300,900);
	c1 -> Divide(3,2);

	c1 -> cd(1);
	h1_dpp_v_p_et_bins[0] -> Draw();
	for(int et = 0 ; et < size_eta_bin-1 ; et++) h1_dpp_v_p_et_bins[et] -> Draw("same");
	c1 -> cd(2);
	h1_dth_v_p_et_bins[0] -> Draw();
	for(int et = 0 ; et < size_eta_bin-1 ; et++) h1_dth_v_p_et_bins[et] -> Draw("same");
	c1 -> cd(3);
	h1_dph_v_p_et_bins[0] -> Draw();
	for(int et = 0 ; et < size_eta_bin-1 ; et++) h1_dph_v_p_et_bins[et] -> Draw("same");
	c1 -> cd(4);
	h1_dpp_v_et_p_bins[0] -> Draw();
	for(int p = 0 ; p < size_mom_bin-1 ; p++) h1_dpp_v_et_p_bins[p] -> Draw("same");
	c1 -> cd(5);
	h1_dth_v_et_p_bins[0] -> Draw();
	for(int p = 0 ; p < size_mom_bin-1 ; p++) h1_dth_v_et_p_bins[p] -> Draw("same");
	c1 -> cd(6);
	h1_dph_v_et_p_bins[0] -> Draw();
	for(int p = 0 ; p < size_mom_bin-1 ; p++) h1_dph_v_et_p_bins[p] -> Draw("same");
	TLegend * leg1 = new TLegend(0.50,0.6,0.70,0.85);
	leg1 -> SetLineColor(0);
	for(int et = 0 ; et < size_eta_bin-1 ; et++) leg1 -> AddEntry(h1_dph_v_p_et_bins[et],Form("%.1f < |#eta| < %.1f",eta_bin[et],eta_bin[et+1]));
	c1 -> cd(3);
	leg1 -> Draw("same");
	TLegend * leg2 = new TLegend(0.20,0.6,0.40,0.85);
	leg2 -> SetLineColor(0);
	for(int p = 0 ; p < size_mom_bin-1 ; p++) leg2 -> AddEntry(h1_dph_v_et_p_bins[p],Form("%.1f < p < %.1f GeV/c",mom_bin[p],mom_bin[p+1]));
	c1 -> cd(6);
	leg2 -> Draw("same");
	//c1->Print("out.pdf");

	// -------------------------------------------------------------
	// Saving fits to pdf
/*	TString out_pdf_name = "fits_AllS_"+partic+"_B_"+Form("%.1f",Bfield)+"T.pdf";

	for(int et = 0 ; et < size_eta_bin-1 ; et++){
		TString fname = out_pdf_name;
		//if(et == 0) fname+="(";
		//else if(et == size_eta_bin-2) fname+=")";
		c_fits_p[et] -> Print(fname);
	}
*/
	// -------------------------------------------------------------

	TFile * Fout = new TFile("histos_"+partic+"_B_"+Form("%.1f",Bfield)+"T.root","recreate");
	for(int et = 0 ; et < size_eta_bin-1 ; et++){
		h1_dpp_v_p_et_bins[et] -> Write(Form("h1_dpp_v_p_et_bins_%i",et));
		h1_dth_v_p_et_bins[et] -> Write(Form("h1_dth_v_p_et_bins_%i",et));
		h1_dph_v_p_et_bins[et] -> Write(Form("h1_dph_v_p_et_bins_%i",et));
		h1_dca2d_v_p_et_bins[et] -> Write(Form("h1_dca2d_v_p_et_bins_%i",et));
		h1_dca2d_v_pT_et_bins[et] -> Write(Form("h1_dca2d_v_pT_et_bins_%i",et));
	}
	for(int p = 0 ; p < size_mom_bin-1 ; p++){
		h1_dpp_v_et_p_bins[p] -> Write(Form("h1_dpp_v_et_p_bins_%i",p));
		h1_dth_v_et_p_bins[p] -> Write(Form("h1_dth_v_et_p_bins_%i",p));
		h1_dph_v_et_p_bins[p] -> Write(Form("h1_dph_v_et_p_bins_%i",p));
		h1_dca2d_v_et_p_bins[p] -> Write(Form("h1_dca2d_v_et_p_bins_%i",p));
		h1_dca2d_v_et_pT_bins[p] -> Write(Form("h1_dca2d_v_et_pT_bins_%i",p));
		
	}

	for(int et = 0 ; et < size_eta_bin-1 ; et++){
	  for(int p = 0 ; p < size_mom_bin-1 ; p++){

	    h1_dpp_p_et_bins[et][p]->Write((Form("h1_dpp_p_et_bins_%i_%i",et,p)));
	    h1_dth_p_et_bins[et][p]->Write((Form("h1_dth_p_et_bins_%i_%i",et,p)));
	    h1_dph_p_et_bins[et][p]->Write((Form("h1_dph_p_et_bins_%i_%i",et,p)));
	    h1_dca2d_p_et_bins[et][p]->Write((Form("h1_dca2d_p_et_bins_%i_%i",et,p)));
	    h1_dca2d_pT_et_bins[et][p]->Write((Form("h1_dca2d_pT_et_bins_%i_%i",et,p)));

// Below block of code to output pdfs of the 4 resolutions along with fitted functions and Chi2/NDF


            gStyle->SetOptFit(1111);
            TCanvas *temp_th = new TCanvas(Form("dthCanvas_eta_%i_p_%i", et, p), Form("dthCanvas_eta_%i_p_%i", et, p));
            temp_th->cd();
            h1_dth_p_et_bins[et][p]->Draw();
            f_gaus1_dth[et][p]->Draw("SAME");
            f_gaus2_dth[et][p]->Draw("SAME");
            TLatex thChi2 = TLatex(0.0005, 0.7*h1_dth_p_et_bins[et][p]->GetMaximum(), Form("#chi^{2} = %.2f", f_gaus_dth[et][p] -> GetChisquare()));
            TLatex thNDF = TLatex(0.0005, 0.65*h1_dth_p_et_bins[et][p]->GetMaximum(), Form("NDF = %.2d", f_gaus_dth[et][p] -> GetNDF()));
            thChi2.Draw("SAME");
            thNDF.Draw("SAME");
            TLatex ths1A1 = TLatex(0.0005/2, 0.6*h1_dth_p_et_bins[et][p]->GetMaximum(), Form("#sigma_{1} = %.1e; A1 = %.1e", sigma1_dth[et][p], A1_dth[et][p]));
            TLatex ths2A2 = TLatex(0.0005/2, 0.55*h1_dth_p_et_bins[et][p]->GetMaximum(), Form("#sigma_{2} = %.1e; A2 = %.1e", sigma2_dth[et][p], A2_dth[et][p]));
            ths1A1.Draw("SAME");
            ths2A2.Draw("SAME");

            if(et==0 && p==0) temp_th->Print("dth.pdf(");
            else if(et==size_eta_bin-2 && p==size_mom_bin-2) temp_th->Print("dth.pdf)");
            else temp_th->Print("dth.pdf");



            gStyle->SetOptFit(1111);
            TCanvas *temp_p = new TCanvas(Form("dppCanvas_eta_%i_p_%i", et, p), Form("dppCanvas_eta_%i_p_%i", et, p));
            temp_p->cd();
	    h1_dpp_p_et_bins[et][p]->GetXaxis()->SetRangeUser(-8*h1_dpp_p_et_bins[et][p]->GetRMS(), 8*h1_dpp_p_et_bins[et][p]->GetRMS());
            h1_dpp_p_et_bins[et][p]->Draw();
	    f_gaus1_dpp[et][p]->Draw("SAME");
	    f_gaus2_dpp[et][p]->Draw("SAME");
	    TLatex Chi2 = TLatex(h1_dpp_p_et_bins[et][p]->GetRMS()/2, 0.7*h1_dpp_p_et_bins[et][p]->GetMaximum(), Form("#chi^{2} = %.2f", f_gaus_dpp[et][p] -> GetChisquare()));
	    TLatex NDF = TLatex(h1_dpp_p_et_bins[et][p]->GetRMS()/2, 0.65*h1_dpp_p_et_bins[et][p]->GetMaximum(), Form("NDF = %.2d", f_gaus_dpp[et][p] -> GetNDF()));
	    Chi2.Draw("SAME");
	    NDF.Draw("SAME");
            TLatex pps1A1 = TLatex(h1_dpp_p_et_bins[et][p]->GetRMS()/4, 0.6*h1_dpp_p_et_bins[et][p]->GetMaximum(), Form("#sigma_{1} = %.1e; A1 = %.1e", sigma1_dpp[et][p], A1_dpp[et][p]));
            TLatex pps2A2 = TLatex(h1_dpp_p_et_bins[et][p]->GetRMS()/4, 0.55*h1_dpp_p_et_bins[et][p]->GetMaximum(), Form("#sigma_{2} = %.1e; A2 = %.1e", sigma2_dpp[et][p], A2_dpp[et][p]));
            pps1A1.Draw("SAME");
            pps2A2.Draw("SAME");
            if(et==0 && p==0) temp_p->Print("dpp.pdf(");
            else if(et==size_eta_bin-2 && p==size_mom_bin-2) temp_p->Print("dpp.pdf)");
            else temp_p->Print("dpp.pdf");


	    gStyle->SetOptFit(1111);
	    TCanvas *temp1 = new TCanvas(Form("dphCanvas_eta_%i_p_%i", et, p), Form("dphCanvas_eta_%i_p_%i", et, p));
            temp1->cd();
	    h1_dph_p_et_bins[et][p]->GetXaxis()->SetRangeUser(-6*h1_dph_p_et_bins[et][p]->GetRMS(), 6*h1_dph_p_et_bins[et][p]->GetRMS());
            h1_dph_p_et_bins[et][p]->Draw();
	    f_gaus1_dph[et][p]->Draw("SAME");
	    f_gaus2_dph[et][p]->Draw("SAME");
	    TLatex phChi2 = TLatex(h1_dph_p_et_bins[et][p]->GetRMS()/2, 0.7*h1_dph_p_et_bins[et][p]->GetMaximum(), Form("#chi^{2} = %.2f", f_gaus_dph[et][p] -> GetChisquare()));
	    TLatex phNDF = TLatex(h1_dph_p_et_bins[et][p]->GetRMS()/2, 0.65*h1_dph_p_et_bins[et][p]->GetMaximum(), Form("NDF = %.2d", f_gaus_dph[et][p] -> GetNDF()));
	    phChi2.Draw("SAME");
	    phNDF.Draw("SAME");
            TLatex phs1A1 = TLatex(h1_dph_p_et_bins[et][p]->GetRMS()/4, 0.6*h1_dph_p_et_bins[et][p]->GetMaximum(), Form("#sigma_{1} = %.1e; A1 = %.1e", sigma1_dph[et][p], A1_dph[et][p]));
            TLatex phs2A2 = TLatex(h1_dph_p_et_bins[et][p]->GetRMS()/4, 0.55*h1_dph_p_et_bins[et][p]->GetMaximum(), Form("#sigma_{2} = %.1e; A2 = %.1e", sigma2_dph[et][p], A2_dph[et][p]));
            phs1A1.Draw("SAME");
            phs2A2.Draw("SAME");
            if(et==0 && p==0) temp1->Print("dph.pdf(");
            else if(et==size_eta_bin-2 && p==size_mom_bin-2) temp1->Print("dph.pdf)");
            else temp1->Print("dph.pdf");
	    
	    gStyle->SetOptFit(1111);
	    TCanvas *temp = new TCanvas(Form("dca2dCanvas_eta_%i_p_%i", et, p), Form("dca2dCanvas_eta_%i_p_%i", et, p));
            temp->cd();
            h1_dca2d_p_et_bins[et][p]->Draw();
	    f_gaus1_dca2d[et][p]->Draw("SAME");
	    f_gaus2_dca2d[et][p]->Draw("SAME");
	    TLatex dca2dChi2 = TLatex(0.004, 0.7*h1_dca2d_p_et_bins[et][p]->GetMaximum(), Form("#chi^{2} = %.2f", f_gaus_dca2d[et][p] -> GetChisquare()));
	    TLatex dca2dNDF = TLatex(0.004, 0.65*h1_dca2d_p_et_bins[et][p]->GetMaximum(), Form("NDF = %.2d", f_gaus_dca2d[et][p] -> GetNDF()));
	    dca2dChi2.Draw("SAME");
	    dca2dNDF.Draw("SAME");
            TLatex dca2ds1A1 = TLatex(0.004/2, 0.6*h1_dca2d_p_et_bins[et][p]->GetMaximum(), Form("#sigma_{1} = %.1e; A1 = %.1e", sigma1_dca2d[et][p], A1_dca2d[et][p]));
            TLatex dca2ds2A2 = TLatex(0.004/2, 0.55*h1_dca2d_p_et_bins[et][p]->GetMaximum(), Form("#sigma_{2} = %.1e; A2 = %.1e", sigma2_dca2d[et][p], A2_dca2d[et][p]));
            dca2ds1A1.Draw("SAME");
            dca2ds2A2.Draw("SAME");
            if(et==0 && p==0) temp->Print("dca2d.pdf(");
            else if(et==size_eta_bin-2 && p==size_mom_bin-2) temp->Print("dca2d.pdf)");
            else temp->Print("dca2d.pdf");
	    
//end of the pdf creation
	  }
	}

	c1 -> Write("c1");
	Fout -> Close();

}
// ============================================================================================================================================
void prettyTH1F( TH1F * h1 , int color , int marker , float min , float max ){
	h1 -> SetLineWidth(2);
	h1 -> SetLineColor(color);
	h1 -> SetMarkerStyle(marker);
	h1 -> SetMarkerColor(color);

	h1 -> SetMinimum(min);
	h1 -> SetMaximum(max);

	h1 -> GetXaxis() -> CenterTitle();
	h1 -> GetXaxis() -> SetNdivisions(107); // to draw less tick marks
	h1 -> GetYaxis() -> CenterTitle();
	h1 -> GetYaxis() -> SetNdivisions(107); // to draw less tick marks

	h1 -> SetMinimum(0.001);
}

