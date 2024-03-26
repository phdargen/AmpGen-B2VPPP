void Fig4_5()
{
//=========Macro generated from canvas: c1_n5/c1_n5
//=========  (Tue Mar 26 22:55:10 2024) by ROOT version 6.18/00
   TCanvas *c1_n5 = new TCanvas("c1_n5", "c1_n5",0,0,700,500);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1_n5->SetHighLightColor(2);
   c1_n5->Range(-0.4657055,-0.09168305,2.860762,0.4813361);
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
   0.3922323,
   0.5547002,
   0.6793662,
   0.7844645,
   0.877058,
   0.9607689,
   1.037749,
   1.1094,
   1.176697,
   1.240347,
   1.300887,
   1.358732,
   1.414214,
   1.467599,
   1.519109,
   1.568929,
   1.617215,
   1.664101,
   1.709701,
   1.754116,
   1.797434,
   1.839732,
   1.88108,
   1.921538,
   1.961161,
   2,
   2.038099,
   2.075498,
   2.112235,
   2.148345,
   2.183857,
   2.218801,
   2.253203,
   2.287087,
   2.320477,
   2.353394,
   2.385856,
   2.417882,
   2.44949};
   Double_t Graph0_fy5[40] = {
   0,
   0,
   0,
   0.001988774,
   0.01118664,
   0.02237248,
   0.03378679,
   0.04527256,
   0.05906538,
   0.07633543,
   0.09933756,
   0.130309,
   0.1624176,
   0.1915895,
   0.217828,
   0.2416184,
   0.2633599,
   0.2833526,
   0.3018267,
   0.318965,
   0.3349171,
   0.349808,
   0.3637439,
   0.3768163,
   0.3891048,
   0.4006792,
   0.4115319,
   0.4115319,
   0.4115319,
   0.4115319,
   0.4115319,
   0.4115319,
   0.4115319,
   0.4115319,
   0.4115319,
   0.4115319,
   0.4115319,
   0.4115319,
   0.4115319,
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
   
   TH1F *Graph_Graph05 = new TH1F("Graph_Graph05","",100,0,2.694439);
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
