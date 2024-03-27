void Fig4_5()
{
//=========Macro generated from canvas: c1_n5/c1_n5
//=========  (Wed Mar 27 12:40:45 2024) by ROOT version 6.18/00
   TCanvas *c1_n5 = new TCanvas("c1_n5", "c1_n5",0,0,700,500);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1_n5->SetHighLightColor(2);
   c1_n5->Range(-0.3802469,-0.09168305,2.335802,0.4813361);
   c1_n5->SetFillColor(0);
   c1_n5->SetBorderMode(0);
   c1_n5->SetBorderSize(2);
   c1_n5->SetTickx(1);
   c1_n5->SetTicky(1);
   c1_n5->SetLeftMargin(0.14);
   c1_n5->SetRightMargin(0.05);
   c1_n5->SetTopMargin(0.05);
   c1_n5->SetBottomMargin(0.16);
   c1_n5->SetFrameLineWidth(2);
   c1_n5->SetFrameBorderMode(0);
   c1_n5->SetFrameLineWidth(2);
   c1_n5->SetFrameBorderMode(0);
   
   Double_t Graph0_fx5[40] = {
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
   Double_t Graph0_fy5[40] = {
   0,
   0,
   0,
   0,
   0.0002164321,
   0.004613326,
   0.01118664,
   0.01857153,
   0.02618888,
   0.03378679,
   0.04135452,
   0.04950227,
   0.05906538,
   0.07012558,
   0.08315176,
   0.09933756,
   0.1194183,
   0.1412464,
   0.1624176,
   0.1822147,
   0.2006371,
   0.217828,
   0.2339333,
   0.2490767,
   0.2633599,
   0.2768673,
   0.2896695,
   0.3018267,
   0.3133908,
   0.3244075,
   0.3349171,
   0.3449554,
   0.3545546,
   0.3637439,
   0.3725496,
   0.380996,
   0.3891048,
   0.3968965,
   0.4043895,
   0.4115319};
   TGraph *graph = new TGraph(40,Graph0_fx5,Graph0_fy5);
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
   
   TH1F *Graph_Graph05 = new TH1F("Graph_Graph05","",100,0,2.2);
   Graph_Graph05->SetMinimum(0);
   Graph_Graph05->SetMaximum(0.4526851);
   Graph_Graph05->SetDirectory(0);
   Graph_Graph05->SetStats(0);
   Graph_Graph05->SetLineWidth(2);
   Graph_Graph05->SetMarkerStyle(8);
   Graph_Graph05->SetMarkerSize(0.5);
   Graph_Graph05->GetXaxis()->SetTitle(" #sqrt{#it{s}} [GeV]  ");
   Graph_Graph05->GetXaxis()->SetNdivisions(505);
   Graph_Graph05->GetXaxis()->SetLabelFont(132);
   Graph_Graph05->GetXaxis()->SetLabelSize(0.065);
   Graph_Graph05->GetXaxis()->SetTitleSize(0.08);
   Graph_Graph05->GetXaxis()->SetTitleOffset(0.9);
   Graph_Graph05->GetXaxis()->SetTitleFont(132);
   Graph_Graph05->GetYaxis()->SetTitle(" #sqrt{#it{s}} / #it{m_{0} #it{#Gamma(s)}} [GeV]");
   Graph_Graph05->GetYaxis()->SetLabelFont(132);
   Graph_Graph05->GetYaxis()->SetLabelSize(0.065);
   Graph_Graph05->GetYaxis()->SetTitleSize(0.08);
   Graph_Graph05->GetYaxis()->SetTitleOffset(0.8);
   Graph_Graph05->GetYaxis()->SetTitleFont(132);
   Graph_Graph05->GetZaxis()->SetLabelFont(132);
   Graph_Graph05->GetZaxis()->SetLabelSize(0.06);
   Graph_Graph05->GetZaxis()->SetTitleSize(0.072);
   Graph_Graph05->GetZaxis()->SetTitleOffset(1.2);
   Graph_Graph05->GetZaxis()->SetTitleFont(132);
   graph->SetHistogram(Graph_Graph05);
   
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
   pt_LaTex = pt->AddText("K^{*}(1680)^{+}");
   pt->Draw();
   c1_n5->Modified();
   c1_n5->cd();
   c1_n5->SetSelected(c1_n5);
}
