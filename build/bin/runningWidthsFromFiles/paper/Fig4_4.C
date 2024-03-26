void Fig4_4()
{
//=========Macro generated from canvas: c1_n4/c1_n4
//=========  (Tue Mar 26 22:55:10 2024) by ROOT version 6.18/00
   TCanvas *c1_n4 = new TCanvas("c1_n4", "c1_n4",0,0,700,500);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1_n4->SetHighLightColor(2);
   c1_n4->Range(-0.3802469,-0.08491086,2.335802,0.445782);
   c1_n4->SetFillColor(0);
   c1_n4->SetBorderMode(0);
   c1_n4->SetBorderSize(2);
   c1_n4->SetTickx(1);
   c1_n4->SetTicky(1);
   c1_n4->SetLeftMargin(0.14);
   c1_n4->SetRightMargin(0.05);
   c1_n4->SetTopMargin(0.05);
   c1_n4->SetBottomMargin(0.16);
   c1_n4->SetFrameLineWidth(2);
   c1_n4->SetFrameBorderMode(0);
   c1_n4->SetFrameLineWidth(2);
   c1_n4->SetFrameBorderMode(0);
   
   Double_t Graph0_fx4[40] = {
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
   Double_t Graph0_fy4[40] = {
   0,
   0,
   0,
   0,
   5.716151e-07,
   0.000106181,
   0.0005236566,
   0.001364172,
   0.002704571,
   0.004663605,
   0.007875402,
   0.01540341,
   0.02333639,
   0.03065264,
   0.0380873,
   0.04655984,
   0.0562557,
   0.06685487,
   0.07802283,
   0.08963942,
   0.1016795,
   0.1141318,
   0.1269807,
   0.1402052,
   0.1537801,
   0.1676782,
   0.1818708,
   0.1963293,
   0.2110253,
   0.2259311,
   0.24102,
   0.2562666,
   0.2716467,
   0.2871373,
   0.3027169,
   0.3183653,
   0.3340635,
   0.349794,
   0.3655404,
   0.381134};
   TGraph *graph = new TGraph(40,Graph0_fx4,Graph0_fy4);
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
   
   TH1F *Graph_Graph04 = new TH1F("Graph_Graph04","",100,0,2.2);
   Graph_Graph04->SetMinimum(0);
   Graph_Graph04->SetMaximum(0.4192474);
   Graph_Graph04->SetDirectory(0);
   Graph_Graph04->SetStats(0);
   Graph_Graph04->SetLineWidth(2);
   Graph_Graph04->SetMarkerStyle(8);
   Graph_Graph04->SetMarkerSize(0.5);
   Graph_Graph04->GetXaxis()->SetTitle(" #sqrt{#it{s}} [GeV]  ");
   Graph_Graph04->GetXaxis()->SetNdivisions(505);
   Graph_Graph04->GetXaxis()->SetLabelFont(132);
   Graph_Graph04->GetXaxis()->SetLabelSize(0.065);
   Graph_Graph04->GetXaxis()->SetTitleSize(0.08);
   Graph_Graph04->GetXaxis()->SetTitleOffset(0.9);
   Graph_Graph04->GetXaxis()->SetTitleFont(132);
   Graph_Graph04->GetYaxis()->SetTitle(" #sqrt{#it{s}} / #it{m_{0} #it{#Gamma(s)}} [GeV]");
   Graph_Graph04->GetYaxis()->SetLabelFont(132);
   Graph_Graph04->GetYaxis()->SetLabelSize(0.065);
   Graph_Graph04->GetYaxis()->SetTitleSize(0.08);
   Graph_Graph04->GetYaxis()->SetTitleOffset(0.8);
   Graph_Graph04->GetYaxis()->SetTitleFont(132);
   Graph_Graph04->GetZaxis()->SetLabelFont(132);
   Graph_Graph04->GetZaxis()->SetLabelSize(0.06);
   Graph_Graph04->GetZaxis()->SetTitleSize(0.072);
   Graph_Graph04->GetZaxis()->SetTitleOffset(1.2);
   Graph_Graph04->GetZaxis()->SetTitleFont(132);
   graph->SetHistogram(Graph_Graph04);
   
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
   pt_LaTex = pt->AddText("K_{2}^{*}(1430)^{+}");
   pt->Draw();
   c1_n4->Modified();
   c1_n4->cd();
   c1_n4->SetSelected(c1_n4);
}
