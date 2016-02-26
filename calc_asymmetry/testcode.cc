// Program to test LHAPDF6 PDF behaviour by writing out their values at lots of x and Q points
// Note: the OpenMP directives are there as an example. In fact, in this case OpenMP slows things
// down because of the need to make the stream operations critical!
#include "LHAPDF/LHAPDF.h"

#include <iostream>
#include <fstream>
#include <cmath>

#include <vector>

#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TProfile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"

#define LogDebug(message) std::cout << "DEBUG: "<< "  " << __FUNCTION__ << "  " << __LINE__ << " :::  " << message << std::endl

using namespace LHAPDF;
using namespace std;

const double G_AHAT = 1; //ahat

const double G_YRANGEMIN = -0.1;
const double G_YRANGEMAX =  0.1;

const double G_Q2 = 10;
const double G_CL_Scale = 2.;

string G_binning_type = "pT";

vector<TH1D*> G_v_A_LL_replicas;

vector<double> G_v_weighting_factors;

TGraphErrors *G_g_run13_data;



double quadSum(double x, double y)
{
	return sqrt(x*x+y*y);
}

double getMean(vector<double> v)
{
	int n = v.size();
	double mean = 0;
	BOOST_FOREACH(double x, v)
	{
		mean += x;
	}
	return mean/n;
}

double getMean(vector<double> v,  vector<double> vw)
{
	int n = v.size();
	double mean = 0;
	double norminal = 0;
	//BOOST_FOREACH(double x, v)
	for(int i=0;i<n;i++)
	{
		double x = v[i];
		double w = vw[i];
		mean += x*w;
		norminal += w;
	}
	return mean/norminal;
}

double getStdDev(vector<double> v, double mean)
{
	int n = v.size();
	double std_dev = 0;
	BOOST_FOREACH(double x, v)
	{
		std_dev += (x-mean)*(x-mean);
	}

	return sqrt(std_dev/(n-1));
}

double getStdDev(vector<double> v, double mean, vector<double> vw)
{
	int n = v.size();
	double std_dev = 0;
	double norminal = 0;
	for(int i=0;i<n;i++)
	{
		double x = v[i];
		double w = vw[i];

		std_dev += (x-mean)*(x-mean)*w;
		norminal += w;
	}
	norminal *= 1.*(n-1)/n;

	//http://stats.stackexchange.com/questions/6534/how-do-i-calculate-a-weighted-standard-deviation-in-excel
	return sqrt(std_dev/norminal);
}

#include "boost/tuple/tuple.hpp"
#include <algorithm>
#include <stdlib.h>     /* exit, EXIT_FAILURE */
typedef tuple<double,double> tdd;
tdd getLowHighCLBoundry(vector<double>v)
{
	if(v.size()!=100)
		exit (EXIT_FAILURE);

	int cl = 0;
	if(fabs(G_CL_Scale - 1.00) < 0.001) cl = 68;
	if(fabs(G_CL_Scale - 1.64) < 0.001) cl = 90;
	if(fabs(G_CL_Scale - 2.00) < 0.001) cl = 95;

	std::sort(v.begin(),v.end());
	double Low = 0;
	double High = 0;
	if(cl == 68)
	{
		Low = v[15];
		High = v[84];
	}
	if(cl == 90)
	{
		Low = v[4];
		High = v[95];
	}
	if(cl == 95)
	{
		Low = 0.5*(v[1]+v[2]);
		High = 0.5*(v[98]+v[97]);
	}
	return make_tuple(Low,High);
}

TGraphErrors* getDataPlot()
{
	int n = 3;
	double x[] = {1.12,2.79,5.25};
	double y[] = {0.003,0.007,0.053};
	double ex[] = {0,0,0};
	double ey[] = {0.014,0.016,0.029};
	//double ey[] = {0.005,0.005,0.005};
	TGraphErrors*gr = new TGraphErrors(n,x,y,ex,ey);
	return gr;
}

void get_Bjorken_x_from_pT(double & x1, double & x2, double pT)
{
	x1 = 0.045 + 0.017*pT;
	x2 = 0.0013 + 0.000063*pT;
}

//TH1D *plotPDF(const string setname = "NNPDFpol11_100", const int imem = 0, double minlogx = -3, double maxlogx = 0)
//{
//	char *hname = Form("plotPDF_%s_%d_Q2=%2.1fGeV",setname.data(),imem,G_Q2);
//	const int nbin = 100;
//	const double dbin = (maxlogx-minlogx)/nbin;
//	TH1D *h_result = new TH1D(hname,hname,nbin,minlogx,maxlogx);
//	const PDF* pdf = mkPDF(setname, imem);
//	const int pid = 21;
//	for(double logx=minlogx;logx<maxlogx;logx+=dbin)
//	{
//		double x = pow(10,logx);
//		double xf = pdf->xfxQ2(pid, x, G_Q2);
//		h_result->Fill(logx+dbin/2.,xf);
//	}
//
//	delete pdf;
//
//	h_result->GetXaxis()->SetTitle("Log(x)");
//	h_result->GetYaxis()->SetTitle("xf(x)");
//
//	return h_result;
//}

TH2D *getAsymmetryEbE(
		const int imem = 0,
		const string binning_type = "pT",
		const string polset = "NNPDFpol11_100",
		const string unpolset = "NNPDF23_nlo_as_0119",
		const string pythia_name = "./input/pythia_input.root"
		)
{
	// Load pythia file
	TFile *fin = TFile::Open(pythia_name.data(),"read");
	TTree *T = (TTree*)fin->Get("T");
	float x1 = 0;
	float x2 = 0;
	float pt = 0;
	float y_com = 0;
	T->SetBranchAddress("x1",&x1);
	T->SetBranchAddress("x2",&x2);
	T->SetBranchAddress("pt",&pt);
	T->SetBranchAddress("y_com",&y_com);


	double min_bin = 1.;
	double max_bin = 10;
	int nbin = 100;

	if(binning_type == "pT")
	{}
	if(binning_type == "y")
	{
		min_bin = 1.2;
		max_bin = 2.2;
	}
	if(binning_type == "None")
	{
		min_bin = 0;
		max_bin = 1;
		nbin = 1;
	}

	double dbin = (max_bin-min_bin)/nbin;

	char *hname = Form("A_LL_%s_pol=%s-unpol=%s-Q2=%2.1f",binning_type.data(),polset.data(),unpolset.data(),G_Q2);

	TH2D *h_var_all = new TH2D(hname,hname,
			nbin,min_bin,max_bin,
			2000,-0.1,0.1);

	const int pid = 21;

	const PDF* polpdf = mkPDF(polset, imem);
	const PDF* unpolpdf = mkPDF(unpolset, 0);

	long nentries = T->GetEntries();


	for(long ientry = 0; ientry<nentries; ientry++)
	{
		T->GetEntry(ientry);

		double x1dg = polpdf->xfxQ2(pid, x1, G_Q2);
		double x2dg = polpdf->xfxQ2(pid, x2, G_Q2);
		double x1g = unpolpdf->xfxQ2(pid, x1, G_Q2);
		double x2g = unpolpdf->xfxQ2(pid, x2, G_Q2);

		double asymmetry = G_AHAT*x1dg*x2dg/x1g/x2g;

		if(binning_type == "pT")
			h_var_all->Fill(pt,asymmetry);
		else if(binning_type == "y")
			h_var_all->Fill(y_com,asymmetry);
		else if(binning_type == "None")
			h_var_all->Fill(0.5,asymmetry);
	}

	delete polpdf;
	delete unpolpdf;

	if(binning_type == "pT")
		h_var_all->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	else if(binning_type == "y")
		h_var_all->GetXaxis()->SetTitle("|y|");
	h_var_all->GetYaxis()->SetTitle("A_{LL}");

	return h_var_all;
}

TH1D *getAsymmetry(
		const int imem = 0,
		const string polset = "NNPDFpol11_100",
		const string unpolset = "NNPDF23_nlo_as_0119"
		)
{
	const double min_bin = 0;
	const double max_bin = 10;
	const int nbin = 100;
	const double dbin = (max_bin-min_bin)/nbin;

	char *hname = Form("pol=%s-unpol=%s-Q2=%2.1f",polset.data(),unpolset.data(),G_Q2);

	TH1D *h_result = new TH1D(hname,hname,nbin,min_bin,max_bin);

	const int pid = 21;

	const PDF* polpdf = mkPDF(polset, imem);
	const PDF* unpolpdf = mkPDF(unpolset, 0);

	for(double pT=min_bin;pT<max_bin;pT+=dbin)
	{
		double x1 = 0;
		double x2 = 0;
		get_Bjorken_x_from_pT(x1,x2,pT);
		double x1dg = polpdf->xfxQ2(pid, x1, G_Q2);
		double x2dg = polpdf->xfxQ2(pid, x2, G_Q2);
		double x1g = unpolpdf->xfxQ2(pid, x1, G_Q2);
		double x2g = unpolpdf->xfxQ2(pid, x2, G_Q2);

		double asymmetry = G_AHAT*x1dg*x2dg/x1g/x2g;

		h_result->Fill(pT+dbin/2.,asymmetry);
	}

	delete polpdf;
	delete unpolpdf;

	h_result->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	h_result->GetYaxis()->SetTitle("A_{LL}");

	return h_result;
}

vector<TH1D*> getAsymmetryReplicaVector(const int nmem = 100)
{
	vector<TH1D*> vrep;

	for(int imem = 0;imem<=nmem;imem++)
	{
		TH1D* h_temp = getAsymmetry(imem);
		vrep.push_back(h_temp);
	}

	return vrep;
}


vector<TH1D*> getAsymmetryEbEReplicaVector(
		const int nmem = 100,
		const string binning_type = "pT"
		)
{
	vector<TH1D*> vrep;

	for(int imem = 1;imem<=nmem;imem++)
	{
		TH1D* h_temp = getAsymmetryEbE(imem,binning_type)->ProfileX()->ProjectionX();
		vrep.push_back(h_temp);
	}

	return vrep;
}

// if g_data is NULL, fill the vector with 1.
vector<double> getWeithtingFactors(vector<TH1D*> v_hists, TGraphErrors* g_data)
{
	vector<double> vw;

	//BOOST_FOREACH(TH1D* hist, v_hists)
	for(int j=0;j<v_hists.size();j++)
	{
		TH1D* hist = v_hists[j];
		double wf = 1;
		if(g_data)
		{
			for(int i=0;i<3;i++)
			{
				double x, mean;
				g_data->GetPoint(i,x,mean);
				double error = g_data->GetErrorY(i);
				double theory = hist->GetBinContent(hist->FindBin(x));

				wf *= TMath::Gaus(theory,mean,error);
				//DEBUG
				cout<<Form("DEBUG: x=%f; \t mean=%f; \t error=%f; \t theory=%f; \t wf=%f; \n",x,mean,error,theory,wf);
			}
		}
		//DEBUG
		cout<<Form("DEBUG: wf=%f; \n",wf);
		vw.push_back(wf);
	}
	return vw;
}


TGraphAsymmErrors* getAsymmetryEbEUncertaintyTGraphErrors(
		const string uncertainty_method = "Counting",
		const double nsigma = 1.,
		TFile* f = NULL,
		const bool do_reweighting = false
		)
{
	TH2D* h2d_temp = getAsymmetryEbE(0,G_binning_type);
	TH1D* h_template = h2d_temp->ProfileX()->ProjectionX();;

	//Show middle step
	if(f)
	{
		h2d_temp->SetName("h2d_asymmetry_pT_y");
		f->cd();
		h2d_temp->Write();
	}

	double xmin = h_template->GetXaxis()->GetXmin();
	double xmax = h_template->GetXaxis()->GetXmax();
	int nbin = h_template->GetXaxis()->GetNbins();

	vector<double> tge_x;
	vector<double> tge_y;
	vector<double> tge_exl;
	vector<double> tge_exh;
	vector<double> tge_eyl;
	vector<double> tge_eyh;

	for(int ibin=1;ibin<=nbin;ibin++)
	{
		vector<double> v_temporary;
		BOOST_FOREACH(TH1D* h, G_v_A_LL_replicas)
		{
			v_temporary.push_back(h->GetBinContent(ibin));
		}
		double mean = h_template->GetBinContent(ibin);
		tge_x.push_back(h_template->GetXaxis()->GetBinCenter(ibin));
		tge_exl.push_back(0);
		tge_exh.push_back(0);
		if(uncertainty_method == "StdDeV")
		{
			double std_dev;
			if(do_reweighting)
			{
				mean = getMean(v_temporary,G_v_weighting_factors);
				std_dev = getStdDev(v_temporary,mean,G_v_weighting_factors);
			}
			else
			{
				mean = getMean(v_temporary);
				std_dev = getStdDev(v_temporary,mean);
			}
			tge_eyl.push_back(nsigma*std_dev);
			tge_eyh.push_back(nsigma*std_dev);
			cout<<Form("DEBUG: ibin: %d; \t mean: %f; \t std_dev: %f; \n",ibin,mean,std_dev);
		}
		if(uncertainty_method == "Counting")
		{
			tdd low_high_boundry = getLowHighCLBoundry(v_temporary);
			tge_eyl.push_back(mean - get<0>(low_high_boundry));
			tge_eyh.push_back(get<1>(low_high_boundry) - mean);
		}
		if(G_binning_type == "None")
		{
			const char* hname = "h_A_LL_dist_Combine";
			TH1D *h_all_dist = new TH1D(hname,hname,100,-0.02,0.02);
			cout<<"====="<<endl;
			cout<<mean<<endl;
			cout<<"-----"<<endl;
			BOOST_FOREACH(double var, v_temporary)
			{
				cout<<var<<endl;
				h_all_dist->Fill(var);
			}
			cout<<"====="<<endl;
			if(f)
			{
				f->cd();
				h_all_dist->Write();
			}

		}
		tge_y.push_back(mean);
	}

	TGraphAsymmErrors* tge = new TGraphAsymmErrors(tge_x.size(),
			&tge_x[0],&tge_y[0],
			&tge_exl[0],&tge_exh[0],
			&tge_eyl[0],&tge_eyh[0]);

	return tge;
}

TH1D* estimateOneValue(
		const double x1 = 0.064,
		const double x2 = 0.0019,
		const string polset = "NNPDFpol11_100",
		const string unpolset = "NNPDF23_nlo_as_0119"
		)
{
	char *hname = Form("pol=%s-unpol=%s-Q2=%2.1f",polset.data(),unpolset.data(),G_Q2);
	TH1D *h_result = new TH1D(hname,hname,100,-0.02,0.02);

	double central_value = -999;

	const int pid = 21;

	PDF* unpolpdf = mkPDF(unpolset, 0);
	PDF* polpdf = NULL;

	const int nmem_pol = 100;
	for(int imem_pol = 0;imem_pol<=nmem_pol;imem_pol++)
	{
		polpdf = mkPDF(polset, imem_pol);
		double x1dg = polpdf->xfxQ2(pid, x1, G_Q2);
		double x2dg = polpdf->xfxQ2(pid, x2, G_Q2);
		double x1g = unpolpdf->xfxQ2(pid, x1, G_Q2);
		double x2g = unpolpdf->xfxQ2(pid, x2, G_Q2);

		double asymmetry = G_AHAT*x1dg*x2dg/x1g/x2g;

		h_result->Fill(asymmetry);

		if(imem_pol==0) central_value = asymmetry;
	}

	char* htitle = Form("%s-Central=%f;A_{LL}^{J/#psi}(x1=%f,x2=%f);",hname,central_value,x1,x2);
	h_result->SetTitle(htitle);

	return h_result;
}

void initStuff()
{
	// G_g_run13_data for the result
	G_g_run13_data = getDataPlot();
	G_g_run13_data->SetMarkerStyle(21);
	G_g_run13_data->SetMarkerColor(kBlack);
	G_g_run13_data->SetLineColor(kBlack);
	if(G_binning_type == "pT")
		G_g_run13_data->SetTitle(";p_{T} GeV; A_{LL}");
	if(G_binning_type == "y")
		G_g_run13_data->SetTitle(";|y|; A_{LL}");
	if(G_binning_type == "None")
		G_g_run13_data->SetTitle(";One Bin; A_{LL}");

	// A_LL calculated by replicas
	G_v_A_LL_replicas = getAsymmetryEbEReplicaVector(100,G_binning_type);//HARD_CODED

	// 
	G_v_weighting_factors = getWeithtingFactors(G_v_A_LL_replicas,G_g_run13_data);
}

//---------------------------------------------------------------------

const double G_Min_Logx = -3;
const double G_Max_Logx = -0;
const int G_N_Points = 100;

TH1D *plotPDF(
		LHAPDF::PDF* pdf
		)
{
	TH1D *h = new TH1D("h","h",G_N_Points,G_Min_Logx,G_Max_Logx);
	double d_logx = (G_Max_Logx - G_Min_Logx)/G_N_Points;
	for(double logx = G_Min_Logx+d_logx/2.; logx<G_Max_Logx; logx += d_logx)
	{
		double x = pow(10.,logx);
		double xf = pdf->xfxQ2(21,x,G_Q2);
		h->Fill(logx,xf);
	}

	return h;
}

TGraphErrors *getLogxXfErrorTH1D(
		string pdf_name = "NNPDFpol11_100",
		double n_scale = 1.,
		const bool do_reweighting = false
		)
{
	const LHAPDF::PDFSet set(pdf_name);

	TH1D* h_central = plotPDF(set.mkPDF(0));

	vector<TH1D*> v_pdf_dist_replicas;
	for(int i=1;i<set.size();i++)
	{
		v_pdf_dist_replicas.push_back(plotPDF(set.mkPDF(i)));
	}

	double xmin = h_central->GetXaxis()->GetXmin();
	double xmax = h_central->GetXaxis()->GetXmax();
	int nbin = h_central->GetXaxis()->GetNbins();

	vector<double> tge_x;
	vector<double> tge_y;
	vector<double> tge_xe;
	vector<double> tge_ye;

	for(int ibin=1;ibin<=nbin;ibin++)
	{
		vector<double> v_temporary;
		BOOST_FOREACH(TH1D* h,v_pdf_dist_replicas)
		{
			v_temporary.push_back(h->GetBinContent(ibin));
		}

		double mean;
		double std_dev;

		if(do_reweighting)
		{
			mean = getMean(v_temporary,G_v_weighting_factors);
			std_dev = getStdDev(v_temporary,mean,G_v_weighting_factors);
		}
		else
		{
			mean = getMean(v_temporary);
			std_dev = getStdDev(v_temporary,mean);
		}

		tge_x.push_back(h_central->GetXaxis()->GetBinCenter(ibin));
		tge_xe.push_back(0);
		tge_y.push_back(mean);
		tge_ye.push_back(n_scale*std_dev);
	}

	TGraphErrors* tge = NULL;

	tge = new TGraphErrors(nbin,&tge_x[0],&tge_y[0],&tge_xe[0],&tge_ye[0]);
	tge->SetName("tge_xf_logx_th1d");
	tge->SetTitle(";log10(x);xf");
	return tge;
}

//---------------------------------------------------------------------


int main(int argc, char* argv[]) {
	if (argc >= 2) {
		G_binning_type = argv[1];
	}

	initStuff();

	TFile *tfile = TFile::Open("./output/out.root","RECREATE");

	TH1D* h_central_EbE = getAsymmetryEbE(0,G_binning_type)->ProfileX()->ProjectionX();
	TGraphAsymmErrors* g_uncertainty_EbE = getAsymmetryEbEUncertaintyTGraphErrors("StdDeV",G_CL_Scale,tfile);

	TLegend* leg = NULL;

	// Main result
	if(true)
	{
		TCanvas *c_main = new TCanvas("c_main","c_main",100,100,1000,800);
		G_g_run13_data->Draw("ap");

		h_central_EbE->SetLineWidth(5);
		h_central_EbE->SetLineColor(46);
		h_central_EbE->GetYaxis()->SetRangeUser(G_YRANGEMIN,G_YRANGEMAX);
		h_central_EbE->Draw("same");
		h_central_EbE->SetName("h_central_EbE");

		g_uncertainty_EbE->SetFillColor(46);
		g_uncertainty_EbE->SetFillStyle(3001);
		g_uncertainty_EbE->Draw("4");
		g_uncertainty_EbE->SetName("g_uncertainty_EbE");

		TGraphAsymmErrors* g_uncertainty_EbE_rw = getAsymmetryEbEUncertaintyTGraphErrors("StdDeV",G_CL_Scale,tfile,true);
		g_uncertainty_EbE_rw->SetFillColor(41);
		g_uncertainty_EbE_rw->SetFillStyle(3004);
		g_uncertainty_EbE_rw->Draw("p");
		g_uncertainty_EbE_rw->SetName("g_uncertainty_EbE_rw");

		leg = new TLegend(0.1,0.6,0.5,0.9);
		leg->AddEntry(h_central_EbE,"Central","l");
		leg->AddEntry(g_uncertainty_EbE,"Uncertianty","f");
		leg->AddEntry(g_uncertainty_EbE_rw,"Re-weighted","lep");
		leg->Draw();


		c_main->Update();
		tfile->cd();
		c_main->Write();
		//c_main->SaveAs("compare.png");
		h_central_EbE->Write();
		g_uncertainty_EbE->Write();
		g_uncertainty_EbE_rw->Write();
	}

	if(true)// show all replicas
	{
		TCanvas *c_show_all_replicas = new TCanvas("c_show_all_replicas","c_show_all_replicas",100,100,1000,800);
		h_central_EbE->Draw();

		BOOST_FOREACH(TH1* h_temp, G_v_A_LL_replicas)
		{
			h_temp->SetLineColor(46);
			h_temp->SetLineStyle(2);
			h_temp->Draw("same");
		}

		g_uncertainty_EbE->SetFillColor(4);
		g_uncertainty_EbE->SetFillStyle(0);
		g_uncertainty_EbE->Draw("4");

		h_central_EbE->Draw("same");
		G_g_run13_data->Draw("samep");
		c_show_all_replicas->Update();

		tfile->cd();
		c_show_all_replicas->Write();
	}

	// Test re-weighting
	if(true)
	{
		TCanvas *c_rw = new TCanvas("c_rw","c_rw",100,100,1000,800);
		TGraphErrors *g_pdf_og = getLogxXfErrorTH1D("NNPDFpol11_100",1.,NULL); // Origninal plot
		g_pdf_og->SetMarkerStyle(20);
		g_pdf_og->SetFillColor(41);
		g_pdf_og->SetFillStyle(3001);
		g_pdf_og->SetName("g_pdf_og");
		g_pdf_og->Draw("ap4");

		TGraphErrors *g_pdf_rw = getLogxXfErrorTH1D("NNPDFpol11_100",1.,G_g_run13_data);// Re-weighted plot
		g_pdf_rw->SetMarkerStyle(4);
		g_pdf_rw->SetFillColor(1);
		g_pdf_rw->SetFillStyle(3005);
		g_pdf_rw->SetName("g_pdf_rw");
		g_pdf_rw->Draw("p");

		TLegend* leg = new TLegend(0.5,0.1,0.9,0.3);
		leg->AddEntry(g_pdf_og,"NNPDFpol1.1 origin","lf");
		leg->AddEntry(g_pdf_rw,"NNPDFpol1.1 re-weighted with PHENIX data","lep");
		leg->Draw();

		c_rw->Update();

		tfile->cd();
		c_rw->Write();
		g_pdf_og->Write();
		g_pdf_rw->Write();

	}


	// One Value
	//TH1D* h_one_value = estimateOneValue();
	//tfile->cd();
	//h_one_value->Write();

	//tfile->Close();

	return 0;
}
