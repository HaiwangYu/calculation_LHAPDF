#include <vector>

#include "math.h"

#include <boost/random.hpp>
#include "boost/tuple/tuple.hpp"

#include "LHAPDF/LHAPDF.h"

#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TProfile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TGraphErrors.h"

using namespace std;
using namespace LHAPDF;
using namespace boost;

typedef tuple<double,double> tdd;

const double G_Min_Logx = -5;
const double G_Max_Logx = -0;
const int G_N_Points = 100;

const double G_Q2 = 4;
const double G_CL = 1.;

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
	//double mean = getMean(v);
	double std_dev = 0;
	BOOST_FOREACH(double x, v)
	{
		std_dev += (x-mean)*(x-mean);
	}

	return sqrt(std_dev/(n-1));
}

tdd getMedianErrStandard(
		double x,
		//string pdf_name = "NNPDFpol11_100",
		LHAPDF::PDFSet set,
		double cl_level = 68
		)
{
	//const LHAPDF::PDFSet set(pdf_name);
	const size_t nmem = set.size()-1;
	const vector<LHAPDF::PDF*> pdfs = set.mkPDFs();
	vector<double> xgAll, xuAll;
	vector<string> pdftypes;
	for (size_t imem = 0; imem <= nmem; imem++) {
		xgAll.push_back(pdfs[imem]->xfxQ2(21,x,G_Q2)); // gluon distribution
		pdftypes.push_back(pdfs[imem]->type()); // PdfType of member
	}

	//set._checkPdfType(pdftypes);

	// Calculate PDF uncertainty on gluon distribution.
	const LHAPDF::PDFUncertainty xgErr = set.uncertainty(xgAll, cl_level); // -1 => same C.L. as set
	return make_tuple(
			xgErr.central,
			xgErr.errsymm
			);
}

TGraphErrors *getLogxXfErrorStandard(
		string pdf_name = "NNPDFpol11_100",
		double cl_level = 1.
		)
{
	vector<double>a_logx;
	vector<double>a_logxe;
	vector<double>a_xf;
	vector<double>a_xfe;

	double d_logx = (G_Max_Logx - G_Min_Logx)/G_N_Points;
	for(double logx = G_Min_Logx+d_logx/2.; logx<G_Max_Logx; logx += d_logx)
	{
		a_logx.push_back(logx);
		a_logxe.push_back(0.0);
		double x = pow(10.,logx);
		const LHAPDF::PDFSet set(pdf_name);
		tdd xf_err = getMedianErrStandard(x,set,cl_level);
		a_xf.push_back(get<0>(xf_err));
		a_xfe.push_back(get<1>(xf_err));
		printf("logx = %e, xf = %e, err = %e \n",logx, get<0>(xf_err),get<1>(xf_err));
	}

	int N = a_logx.size();

	TGraphErrors *g = new TGraphErrors(N,&a_logx[0],&a_xf[0],&a_logxe[0],&a_xfe[0]);
	g->SetName("tge_logx_xf_standard");
	g->SetTitle(";log10(x);xf");
	return g;
}

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
		double n_scale = 1.
		)
{
	const LHAPDF::PDFSet set(pdf_name);
	//const vector<LHAPDF::PDF*> pdfs_replicas = set.mkPDFs();
	vector<TH1D*> v_th1d_replicas;
	//BOOST_FOREACH(LHAPDF::PDF*pdf, pdfs_replicas)
	for(int i=0;i<set.size();i++)
	{
		v_th1d_replicas.push_back(plotPDF(set.mkPDF(i)));
	}
	TH1D* h_central = v_th1d_replicas[0];

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

		for(int irep=1;irep<v_th1d_replicas.size();irep++)
		{
			TH1D* h = v_th1d_replicas[irep];
			v_temporary.push_back(h->GetBinContent(ibin));
		}
		double mean = h_central->GetBinContent(ibin);
		double std_dev = getStdDev(v_temporary,mean);

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

TH1D *plotPDFRatio(
		const LHAPDF::PDF* pdf1,
		const LHAPDF::PDF* pdf2
		)
{
	TH1D *h = new TH1D("h","h",G_N_Points,G_Min_Logx,G_Max_Logx);
	double d_logx = (G_Max_Logx - G_Min_Logx)/G_N_Points;
	for(double logx = G_Min_Logx+d_logx/2.; logx<G_Max_Logx; logx += d_logx)
	{
		double x = pow(10.,logx);
		double xf1 = pdf1->xfxQ2(21,x,G_Q2);
		double xf2 = pdf2->xfxQ2(21,x,G_Q2);

		h->Fill(logx,xf1/xf2);
	}

	return h;
}
TGraphErrors *getLogxXfRatioErrorTH1D(
		string pdf1_name = "NNPDFpol11_100",
		string pdf2_name = "NNPDF23_nlo_as_0119",
		double n_scale = 1.
		)
{
	const LHAPDF::PDFSet set(pdf1_name);
	//const vector<LHAPDF::PDF*> pdf1s = set.mkPDFs();
	const LHAPDF::PDF *pdf2 = mkPDF(pdf2_name,0);
	vector<TH1D*> v_th1d_replicas;
	//BOOST_FOREACH(LHAPDF::PDF*pdf1, pdf1s)
	for(int i=0;i<set.size();i++)
	{
		v_th1d_replicas.push_back(plotPDFRatio(set.mkPDF(i),pdf2));
	}
	TH1D* h_central = v_th1d_replicas[0];

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

		//BOOST_FOREACH(TH1D* h, v_th1d_replicas)
		for(int irep=1;irep<v_th1d_replicas.size();irep++)
		{
			TH1D* h = v_th1d_replicas[irep];
			v_temporary.push_back(h->GetBinContent(ibin));
		}
		double mean = h_central->GetBinContent(ibin);
		double std_dev = getStdDev(v_temporary,mean);

		tge_x.push_back(h_central->GetXaxis()->GetBinCenter(ibin));
		tge_xe.push_back(0);
		tge_y.push_back(mean);
		tge_ye.push_back(n_scale*std_dev);
	}

	TGraphErrors* tge = NULL;

	tge = new TGraphErrors(nbin,&tge_x[0],&tge_y[0],&tge_xe[0],&tge_ye[0]);
	tge->SetName("tge_xf_ratio_logx_th1d");
	tge->SetTitle(";log10(x);xf1/xf2");
	return tge;
}

int main(int argc, char* argv[]) {
	TFile *fout = TFile::Open("out.root","recreate");

	TCanvas *c0 = new TCanvas("c0","c0",100,100,1000,800);
	//TGraphErrors *tge_lhapdf_standard = getLogxXfErrorStandard("NNPDFpol11_100",95);
	//tge_lhapdf_standard->SetFillColor(41);
	//tge_lhapdf_standard->SetFillStyle(3004);
	//tge_lhapdf_standard->SetMarkerColor(41);
	//tge_lhapdf_standard->SetMarkerStyle(20);
	//tge_lhapdf_standard->Draw("ap4");
	//tge_lhapdf_standard->Write();

	if(true)
	{
		//const char* set_name = "NNPDFpol10_100";
		//const char* set_name = "NNPDFpol11_100";
		const char* set_name = "NNPDF23_nlo_as_0119";
		TGraphErrors *tge_lhapdf_th1d = getLogxXfErrorTH1D(set_name,G_CL);
		tge_lhapdf_th1d->SetFillColor(40);
		tge_lhapdf_th1d->SetFillStyle(3005);
		tge_lhapdf_th1d->SetMarkerColor(40);
		tge_lhapdf_th1d->SetMarkerStyle(20);
		tge_lhapdf_th1d->GetYaxis()->SetRangeUser(-1.15,0.75);
		tge_lhapdf_th1d->SetTitle(Form("%s: Q2 = %1.f, C.L. scale = %1.f",set_name,G_Q2,G_CL));
		tge_lhapdf_th1d->Draw("ap");
		tge_lhapdf_th1d->Write();
	}
	else
	{
		//const char* set1_name = "NNPDFpol10_100";
		const char* set1_name = "NNPDFpol11_100";
		const char* set2_name = "NNPDF23_nlo_as_0119";
		//const char* set2_name = "NNPDF23_nnlo_as_0119";
		//const char* set2_name = "NNPDF23_nlo_as_0119_mc";
		TGraphErrors *tge_ratio_th1d = getLogxXfRatioErrorTH1D(set1_name,set2_name,G_CL);
		tge_ratio_th1d->SetFillColor(41);
		tge_ratio_th1d->SetFillStyle(3001);
		tge_ratio_th1d->SetMarkerColor(41);
		tge_ratio_th1d->SetMarkerStyle(20);
		tge_ratio_th1d->GetYaxis()->SetRangeUser(-0.22,0.42);
		tge_ratio_th1d->SetTitle(Form("%s/%s: Q2 = %1.f, C.L. scale = %1.f",set1_name,set2_name,G_Q2,G_CL));
		tge_ratio_th1d->Draw("ap");
		tge_ratio_th1d->Write();
	}


	c0->Write();
	c0->SaveAs("c0.pdf");

	return 0;
}

