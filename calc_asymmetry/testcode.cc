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

using namespace LHAPDF;
using namespace std;

const double G_AHAT = 1; //ahat

const double G_YRANGEMIN = -0.1;
const double G_YRANGEMAX =  0.1;

const double G_Q2 = 10;
const double G_CL_Scale = 2.;

string G_binning_type = "pT";

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
	TGraphErrors*gr = new TGraphErrors(n,x,y,ex,ey);
	return gr;
}

void get_Bjorken_x_from_pT(double & x1, double & x2, double pT)
{
	x1 = 0.045 + 0.017*pT;
	x2 = 0.0013 + 0.000063*pT;
}

TH1D *plotPDF(const string setname = "NNPDFpol11_100", const int imem = 0, double minlogx = -3, double maxlogx = 0)
{
	char *hname = Form("plotPDF_%s_%d_Q2=%2.1fGeV",setname.data(),imem,G_Q2);
	const int nbin = 100;
	const double dbin = (maxlogx-minlogx)/nbin;
	TH1D *h_result = new TH1D(hname,hname,nbin,minlogx,maxlogx);
	const PDF* pdf = mkPDF(setname, imem);
	const int pid = 21;
	for(double logx=minlogx;logx<maxlogx;logx+=dbin)
	{
		double x = pow(10,logx);
		double xf = pdf->xfxQ2(pid, x, G_Q2);
		h_result->Fill(logx+dbin/2.,xf);
	}

	delete pdf;

	h_result->GetXaxis()->SetTitle("Log(x)");
	h_result->GetYaxis()->SetTitle("xf(x)");

	return h_result;
}

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

vector<TH1D*> getAsymmetryUncertaintyVector(const int nmem = 100)
{
	vector<TH1D*> v_eigen;

	for(int imem = 0;imem<=nmem;imem++)
	{
		TH1D* h_temp = getAsymmetry(imem);
		v_eigen.push_back(h_temp);
	}

	return v_eigen;
}


vector<TH1D*> getAsymmetryEbEUncertaintyVector(
		const int nmem = 100,
		const string binning_type = "pT"
		)
{
	vector<TH1D*> v_eigen;

	for(int imem = 0;imem<=nmem;imem++)
	{
		TH1D* h_temp = getAsymmetryEbE(imem,binning_type)->ProfileX()->ProjectionX();
		v_eigen.push_back(h_temp);
	}

	return v_eigen;
}


TGraphAsymmErrors* getAsymmetryEbEUncertaintyTGraphErrors(
		const int nmem = 100,
		const string uncertainty_method = "Counting",
		const double nsigma = 1.,
		TFile* f = NULL)
{
	vector<TH1D*> v_eigen = getAsymmetryEbEUncertaintyVector(nmem,G_binning_type);
	TH1D* h_template = v_eigen[0];

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
		for(int irep = 1; irep < v_eigen.size(); irep++)
		{
			TH1D* h = v_eigen[irep];
			v_temporary.push_back(h->GetBinContent(ibin));
		}
		double mean = h_template->GetBinContent(ibin);
		tge_x.push_back(h_template->GetXaxis()->GetBinCenter(ibin));
		tge_exl.push_back(0);
		tge_exh.push_back(0);
		tge_y.push_back(mean);
		if(uncertainty_method == "StdDeV")
		{
			double std_dev = getStdDev(v_temporary,mean);
			tge_eyl.push_back(nsigma*std_dev);
			tge_eyh.push_back(nsigma*std_dev);
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
	}

	//TGraphErrors* tge = new TGraphErrors(tge_x.size(),&tge_x[0],&tge_y[0],&tge_ex[0],&tge_ey[0]);
	TGraphAsymmErrors* tge = new TGraphAsymmErrors(tge_x.size(),
			&tge_x[0],&tge_y[0],
			&tge_exl[0],&tge_exh[0],
			&tge_eyl[0],&tge_eyh[0]);

	return tge;
}

TGraphErrors* getAsymmetryUncertaintyTGraphErrors(const int nmem = 100, const int nsigma = 2)
{
	vector<TH1D*> v_eigen = getAsymmetryUncertaintyVector(nmem);
	TH1D* h_template = v_eigen[0];

	double xmin = h_template->GetXaxis()->GetXmin();
	double xmax = h_template->GetXaxis()->GetXmax();
	int nbin = h_template->GetXaxis()->GetNbins();

	vector<double> tge_x;
	vector<double> tge_y;
	vector<double> tge_ex;
	vector<double> tge_ey;

	for(int ibin=1;ibin<=nbin;ibin++)
	{
		vector<double> v_temporary;
		for(int irep = 1; irep < v_eigen.size(); irep++)
		{
			TH1D* h = v_eigen[irep];
			v_temporary.push_back(h->GetBinContent(ibin));
		}
		double mean = h_template->GetBinContent(ibin);
		double std_dev = getStdDev(v_temporary,mean);

		tge_x.push_back(h_template->GetXaxis()->GetBinCenter(ibin));
		tge_ex.push_back(0);
		tge_y.push_back(mean);
		tge_ey.push_back(nsigma*std_dev);
	}

	TGraphErrors* tge = new TGraphErrors(tge_x.size(),&tge_x[0],&tge_y[0],&tge_ex[0],&tge_ey[0]);

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


int main(int argc, char* argv[]) {
	if (argc >= 2) {
		G_binning_type = argv[1];
	}

	TFile *tfile = TFile::Open("./output/out.root","RECREATE");

	// Test Code
	//TH1D* h_test = plotPDF("NNPDF23_nlo_as_0124");
	//TH1D* h_test = plotPDF("MSTW2008lo68cl_nf3");
	//TH1D* h_test = plotPDF("NNPDF23_nlo_as_0119");
	TH1D* h_test = plotPDF("NNPDF23_nlo_as_0119",0);
	tfile->cd();
	h_test->Write();

	//
	TCanvas *c_asymmetry_tge = new TCanvas("c_asymmetry_tge","c_asymmetry_tge",100,100,1000,800);

	TGraphErrors *dataplot = getDataPlot();
	dataplot->SetMarkerStyle(21);
	dataplot->SetMarkerColor(kBlack);
	dataplot->SetLineColor(kBlack);
	if(G_binning_type == "pT")
		dataplot->SetTitle(";p_{T} GeV; A_{LL}");
	if(G_binning_type == "y")
		dataplot->SetTitle(";|y|; A_{LL}");
	if(G_binning_type == "None")
		dataplot->SetTitle(";One Bin; A_{LL}");
	dataplot->Draw("ap");

	//TH1D* h_central = getAsymmetry();
	//h_central->SetLineColor(41);
	//h_central->SetLineWidth(5);
	//h_central->GetYaxis()->SetRangeUser(G_YRANGEMIN,G_YRANGEMAX);
	//h_central->Draw("same");

	//TGraphErrors* g_uncertainty = getAsymmetryUncertaintyTGraphErrors();
	//g_uncertainty->SetFillColor(41);
	//g_uncertainty->SetFillStyle(3004);
	//g_uncertainty->Draw("4");

	TH1D* h_central_EbE = getAsymmetryEbE(0,G_binning_type)->ProfileX()->ProjectionX();
	h_central_EbE->SetLineWidth(5);
	h_central_EbE->SetLineColor(46);
	h_central_EbE->GetYaxis()->SetRangeUser(G_YRANGEMIN,G_YRANGEMAX);
	h_central_EbE->Draw("same");
	tfile->cd();
	h_central_EbE->SetName("h_central_EbE");
	h_central_EbE->Write();

	TGraphAsymmErrors* g_uncertainty_EbE = getAsymmetryEbEUncertaintyTGraphErrors(100,"Counting",G_CL_Scale,tfile);
	g_uncertainty_EbE->SetFillColor(46);
	g_uncertainty_EbE->SetFillStyle(3005);
	g_uncertainty_EbE->Draw("4");
	tfile->cd();
	g_uncertainty_EbE->SetName("g_uncertainty_EbE");
	g_uncertainty_EbE->Write();

	TGraphAsymmErrors* g_uncertainty_EbE_StdDev = getAsymmetryEbEUncertaintyTGraphErrors(100,"StdDeV",G_CL_Scale);
	g_uncertainty_EbE_StdDev->SetFillColor(46);
	g_uncertainty_EbE_StdDev->SetFillStyle(3005);
	g_uncertainty_EbE_StdDev->Draw("4");
	tfile->cd();
	g_uncertainty_EbE_StdDev->SetName("g_uncertainty_EbE_StdDev");
	g_uncertainty_EbE_StdDev->Write();

	c_asymmetry_tge->Update();
	tfile->cd();
	c_asymmetry_tge->Write();
	c_asymmetry_tge->SaveAs("compare.png");

	if(true)// <pT> method result
	{
		//TCanvas *c_asymmetry_vector = new TCanvas("c_asymmetry_vector","c_asymmetry_vector",100,100,1000,800);
		//h_central->Draw();

		//vector<TH1D*> v_uncertainty = getAsymmetryUncertaintyVector();
		//BOOST_FOREACH(TH1* h_temp, v_uncertainty)
		//{
		//	h_temp->SetLineColor(41);
		//	h_temp->SetLineStyle(2);
		//	h_temp->Draw("same");
		//}

		//g_uncertainty->SetFillColor(4);
		//g_uncertainty->SetFillStyle(0);
		//g_uncertainty->Draw("4");

		//h_central->Draw("same");
		//dataplot->Draw("samep");
		//c_asymmetry_vector->Update();

		//tfile->cd();
		//c_asymmetry_vector->Write();
		//c_asymmetry_vector->SaveAs("c_asymmetry_vector.png");
	}

	if(true)// Event-by-event plot all results
	{
		TCanvas *c_asymmetry_EbE_vector = new TCanvas("c_asymmetry_EbE_vector","c_asymmetry_EbE_vector",100,100,1000,800);
		h_central_EbE->Draw();

		vector<TH1D*> v_uncertainty_EbE = getAsymmetryEbEUncertaintyVector(100,G_binning_type);
		BOOST_FOREACH(TH1* h_temp, v_uncertainty_EbE)
		{
			h_temp->SetLineColor(46);
			h_temp->SetLineStyle(2);
			h_temp->Draw("same");
		}

		g_uncertainty_EbE->SetFillColor(4);
		g_uncertainty_EbE->SetFillStyle(0);
		g_uncertainty_EbE->Draw("4");

		h_central_EbE->Draw("same");
		dataplot->Draw("samep");
		c_asymmetry_EbE_vector->Update();

		tfile->cd();
		c_asymmetry_EbE_vector->Write();
		c_asymmetry_EbE_vector->SaveAs("c_asymmetry_EbE_vector.png");
	}

	// One Value
	TH1D* h_one_value = estimateOneValue();
	tfile->cd();
	h_one_value->Write();

	tfile->Close();

	return 0;
}
