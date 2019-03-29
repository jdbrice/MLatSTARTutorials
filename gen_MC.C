TH1 * hGlauber = NULL;
TH1 * hPt = NULL;
TH1 * hEta = NULL;

TH1 * hmc_Mult = NULL;
TH1 * hmc_Pt   = NULL;
TH1 * hmc_Eta  = NULL;

TH2 * hMult = NULL;
TH2 * hMultTOF = NULL;
TH2 * hMultBEMC = NULL;
TH1 * hrc_Mult = NULL;
TH1 * hrc_Pt   = NULL;
TH1 * hrc_Eta  = NULL;

TRandom3 r;

TF1 * fPtRes = NULL;
TF1 * fCrystalBall = NULL;
TF1 * fMultEff = NULL;
TF1 * fVzEff = NULL;
TF1 * fZDC = NULL;

// the tree
TNtuple *ntuple = NULL;

void gen_event(){

	// sample the mult
	size_t mult = hGlauber->GetRandom(); 
	hmc_Mult -> Fill( mult );
	// if ( r.Uniform() < 0.1 )
	// 	mult += hGlauber->GetRandom();

	float vz = r.Gaus( 0, 35 );

	size_t reco_mult = 0;
	size_t tof_mult = 0;
	size_t bemc_mult = 0;

	size_t zdcx = r.Rndm() * 30000 + 40000;
	for ( size_t i = 0; i < mult; i++ ){

		// sample pT particles
		double pt = hPt -> GetRandom() / 4;;
		double eta = hEta -> GetRandom();

		hmc_Pt->Fill( pt );
		hmc_Eta->Fill( eta );


		// blur the momentum 
		float ptRes = fPtRes->Eval(pt);
		float rndCrystalBall = fCrystalBall->GetRandom() * 100;
		float recoPt = pt * (1 + rndCrystalBall * ptRes  );


		if ( r.Binomial( 1, fZDC->Eval( zdcx ) ) == 0 ) continue;
		if ( r.Binomial( 1, fVzEff->Eval( vz ) ) == 0 ) continue;
		
		if ( recoPt > 0.075 && r.Binomial( 1, 0.8 ) == 1 && r.Binomial( 1, fMultEff->Eval((float)mult) ) == 1 ){
			hrc_Pt -> Fill( recoPt );
			hrc_Eta -> Fill( eta );
			reco_mult++;
		}

		if ( pt > 0.2 && r.Binomial( 1, 0.5 ) == 1 ) {
			tof_mult++;
		}

		if ( pt > 1.0 && r.Binomial( 1, 0.90 ) == 1 ){
			bemc_mult ++;
		}
	}

	tof_mult += r.Poisson( 4 );

	hrc_Mult->Fill( reco_mult );
	hMult->Fill( mult, reco_mult );
	hMultTOF->Fill( mult, tof_mult );
	hMultBEMC->Fill( mult, bemc_mult );


	ntuple->Fill( reco_mult, tof_mult, bemc_mult, vz, zdcx, mult );

}


double CrystalBall( double x, double N, double mu, double sig, double n, double alpha ){

	double A = pow( n/fabs(alpha), n) * exp(-alpha*alpha/2.0);
	double B = n/fabs(alpha) - fabs(alpha);
	double norm = (x-mu)/sig;

	if(norm > -alpha) {
		return N * exp(-0.5*norm*norm);
	} else {
		return N * A * pow(B-norm, -n);
	}
}
double CrystalBall( double *x, double *par ){
	return CrystalBall( x[0], par[0], par[1], par[2], par[3], par[4] );
}


void gen_MC( unsigned long int n_events = 500000 ) {

	TFile * fGlauber = new TFile( "fit.root" );
	hGlauber = (TH1*)fGlauber->Get( "hRefMultSim" );

	TFile * fJPsi = new TFile( "CERES_jpsi_mumu.root" );
	hPt = (TH1*)fJPsi->Get("FullAcc_l1PtMc");
	hEta = (TH1*)fJPsi->Get("FullAcc_l1Eta");
	hEta->GetXaxis()->SetRangeUser(-1, 1);

	TFile * fOUT = new TFile( "simMult.root", "RECREATE" );

	// hGlauber->Draw();

	hmc_Pt = new TH1F( "hmc_Pt", ";p_{T} (GeV/c)", 500, 0, 10 );
	hmc_Eta = new TH1F( "hmc_Eta", ";#eta", 500, -5, 5 );
	hmc_Mult = new TH1F( "hmc_Mult", ";mult", 400, 0, 400 );

	hrc_Pt = new TH1F( "hrc_Pt", ";p_{T} (GeV/c)", 500, 0, 10 );
	hrc_Eta = new TH1F( "hrc_Eta", ";#eta", 500, -5, 5 );
	hrc_Mult = new TH1F( "hrc_Mult", ";mult", 400, 0, 400 );
	hMult = new TH2F( "hMult", ";mult MC ; mult RC", 400, 0, 400, 400, 0, 400 );
	hMultTOF = new TH2F( "hMultTof", ";mult MC ; mult RC", 400, 0, 400, 400, 0, 400 );
	hMultBEMC = new TH2F( "hMultBEMC", ";mult MC ; mult RC", 400, 0, 400, 400, 0, 400 );

	fPtRes = new TF1( "fPtRes", "sqrt( [0]*[0]*x*x + [1]*[1])" );
	fPtRes->SetParameters( 6.0e-3, 8.3e-3 );
	fPtRes->SetRange( 0, 100 );
	fPtRes->Write();

	fCrystalBall = new TF1( "fCrystalBall", CrystalBall, -1, 1, 5 );
	fCrystalBall->SetParameters( 1, -1e-3, 0.01, 1.29, 1.75 );
	fCrystalBall->SetNpx(5000);
	fCrystalBall->Write();


	fMultEff = new TF1( "fMultEff", "pol1" );
	fMultEff -> SetParameters( 1, -0.0025 );
	fMultEff -> SetRange( 0, 1000 );
	fMultEff -> Write();


	fVzEff = new TF1( "fVzEff", "[0] - [1] * x * x + [2] * x*x*x" );
	fVzEff -> SetParameters( 1, 0.00005, 0.00000009 );
	fVzEff -> SetRange( -200, 200 );
	fVzEff -> Write();

	fZDC  = new TF1( "fZDC", "1 - (0.05 / 10000)*x" );
	fZDC->SetRange(0, 100000);
	fZDC->Write();


	ntuple = new TNtuple("mult","MULT","reco_mult:tof_mult:bemc_mult:vz:zdcx:mult");

	for ( unsigned long int i = 0; i < n_events; i++ ){
		gen_event();
	}

	ntuple->Write();
	fOUT->Write();
}