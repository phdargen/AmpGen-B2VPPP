void Fig4_1()
{
//=========Macro generated from canvas: c1/c1
//=========  (Tue Mar 26 22:55:10 2024) by ROOT version 6.18/00
   TCanvas *c1 = new TCanvas("c1", "c1",0,0,700,500);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1->SetHighLightColor(2);
   c1->Range(-0.3802469,-0.1908124,2.335802,1.001765);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetTickx(1);
   c1->SetTicky(1);
   c1->SetLeftMargin(0.14);
   c1->SetRightMargin(0.05);
   c1->SetTopMargin(0.05);
   c1->SetBottomMargin(0.16);
   c1->SetFrameLineWidth(2);
   c1->SetFrameBorderMode(0);
   c1->SetFrameLineWidth(2);
   c1->SetFrameBorderMode(0);
   
   Double_t Graph0_fx1[40] = {
   0,
   0.3202563,
   0.4529108,
   0.5547002,
   0.6405126,
   0.7161149,
   0.7844645,
   0.8473185,
   0.9058216,
   0.9607689,
   1.012739,
   1.06217,
   1.1094,
   1.154701,
   1.198289,
   1.240347,
   1.281025,
   1.320451,
   1.358732,
   1.395965,
   1.43223,
   1.467599,
   1.502135,
   1.535895,
   1.568929,
   1.601282,
   1.632993,
   1.664101,
   1.694637,
   1.724633,
   1.754116,
   1.783112,
   1.811643,
   1.839732,
   1.867399,
   1.894662,
   1.921538,
   1.948043,
   1.974192,
   2};
   Double_t Graph0_fy1[40] = {
   0,
   0,
   0,
   0,
   0,
   0,
   7.663856e-08,
   2.692931e-05,
   0.0002140915,
   0.0009378543,
   0.004216748,
   0.0169114,
   0.02972601,
   0.04204034,
   0.05681312,
   0.07812016,
   0.1085612,
   0.1433415,
   0.1769023,
   0.20782,
   0.2364126,
   0.263331,
   0.2892212,
   0.3146646,
   0.3401603,
   0.366106,
   0.3927848,
   0.4203725,
   0.4489627,
   0.4785964,
   0.5092912,
   0.5410683,
   0.5739831,
   0.6081634,
   0.643864,
   0.6815255,
   0.7217482,
   0.7649172,
   0.8104675,
   0.8564874};
   TGraph *graph = new TGraph(40,Graph0_fx1,Graph0_fy1);
   graph->SetName("Graph0");
   graph->SetTitle("; #sqrt{#it{s}} [GeV]  ; #sqrt{#it{s}} / #it{m_{0} #it{#Gamma(s)}} [GeV]");
   graph->SetFillStyle(1000);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#0000ff");
   graph->SetLineColor(ci);
   graph->SetLineWidth(2);
   graph->SetMarkerStyle(3);
   graph->SetMarkerSize(0.5);
   
   TH1F *Graph_Graph01 = new TH1F("Graph_Graph01","",100,0,2.2);
   Graph_Graph01->SetMinimum(0);
   Graph_Graph01->SetMaximum(0.9421361);
   Graph_Graph01->SetDirectory(0);
   Graph_Graph01->SetStats(0);
   Graph_Graph01->SetLineWidth(2);
   Graph_Graph01->SetMarkerStyle(8);
   Graph_Graph01->SetMarkerSize(0.5);
   Graph_Graph01->GetXaxis()->SetTitle(" #sqrt{#it{s}} [GeV]  ");
   Graph_Graph01->GetXaxis()->SetNdivisions(505);
   Graph_Graph01->GetXaxis()->SetLabelFont(132);
   Graph_Graph01->GetXaxis()->SetLabelSize(0.065);
   Graph_Graph01->GetXaxis()->SetTitleSize(0.08);
   Graph_Graph01->GetXaxis()->SetTitleOffset(0.9);
   Graph_Graph01->GetXaxis()->SetTitleFont(132);
   Graph_Graph01->GetYaxis()->SetTitle(" #sqrt{#it{s}} / #it{m_{0} #it{#Gamma(s)}} [GeV]");
   Graph_Graph01->GetYaxis()->SetLabelFont(132);
   Graph_Graph01->GetYaxis()->SetLabelSize(0.065);
   Graph_Graph01->GetYaxis()->SetTitleSize(0.08);
   Graph_Graph01->GetYaxis()->SetTitleOffset(0.8);
   Graph_Graph01->GetYaxis()->SetTitleFont(132);
   Graph_Graph01->GetZaxis()->SetLabelFont(132);
   Graph_Graph01->GetZaxis()->SetLabelSize(0.06);
   Graph_Graph01->GetZaxis()->SetTitleSize(0.072);
   Graph_Graph01->GetZaxis()->SetTitleOffset(1.2);
   Graph_Graph01->GetZaxis()->SetTitleFont(132);
   graph->SetHistogram(Graph_Graph01);
   
   graph->Draw("apc");
   
   TPaveText *pt = new TPaveText(0.19,0.82,0.34,0.9,"BRNDC");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetLineWidth(2);
   pt->SetTextAlign(12);
   pt->SetTextFont(132);
   pt->SetTextSize(0.08);
   TText *pt_LaTex = pt->AddText("LHCb");
   pt->Draw();
   
   pt = new TPaveText(0.19,0.72,0.34,0.8,"BRNDC");
   pt->SetFillColor(0);
   pt->SetLineColor(0);
   pt->SetLineWidth(2);
   pt->SetTextAlign(12);
   pt->SetTextFont(132);
   pt->SetTextSize(0.08);
   pt_LaTex = pt->AddText("K_{1}(1270)^{+}");
   pt->Draw();
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
