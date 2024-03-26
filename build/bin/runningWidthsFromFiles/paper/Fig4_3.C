void Fig4_3()
{
//=========Macro generated from canvas: c1_n3/c1_n3
//=========  (Tue Mar 26 22:55:10 2024) by ROOT version 6.18/00
   TCanvas *c1_n3 = new TCanvas("c1_n3", "c1_n3",0,0,700,500);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1_n3->SetHighLightColor(2);
   c1_n3->Range(-0.3802469,-0.1597846,2.335802,0.8388694);
   c1_n3->SetFillColor(0);
   c1_n3->SetBorderMode(0);
   c1_n3->SetBorderSize(2);
   c1_n3->SetTickx(1);
   c1_n3->SetTicky(1);
   c1_n3->SetLeftMargin(0.14);
   c1_n3->SetRightMargin(0.05);
   c1_n3->SetTopMargin(0.05);
   c1_n3->SetBottomMargin(0.16);
   c1_n3->SetFrameLineWidth(2);
   c1_n3->SetFrameBorderMode(0);
   c1_n3->SetFrameLineWidth(2);
   c1_n3->SetFrameBorderMode(0);
   
   Double_t Graph0_fx3[40] = {
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
   Double_t Graph0_fy3[40] = {
   0,
   0,
   0,
   0,
   3.851992e-05,
   0.0008210656,
   0.001990965,
   0.003306941,
   0.004686378,
   0.006177104,
   0.008215306,
   0.01375264,
   0.02626542,
   0.04452328,
   0.06722408,
   0.09346171,
   0.1224675,
   0.1531851,
   0.1846671,
   0.2163648,
   0.2479462,
   0.2791833,
   0.3099123,
   0.3400157,
   0.369411,
   0.3980432,
   0.4258786,
   0.4528991,
   0.4790994,
   0.5044832,
   0.5290613,
   0.5528494,
   0.5758669,
   0.5981358,
   0.6196796,
   0.6405229,
   0.6606906,
   0.680208,
   0.6990998,
   0.7172152};
   TGraph *graph = new TGraph(40,Graph0_fx3,Graph0_fy3);
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
   
   TH1F *Graph_Graph03 = new TH1F("Graph_Graph03","",100,0,2.2);
   Graph_Graph03->SetMinimum(0);
   Graph_Graph03->SetMaximum(0.7889367);
   Graph_Graph03->SetDirectory(0);
   Graph_Graph03->SetStats(0);
   Graph_Graph03->SetLineWidth(2);
   Graph_Graph03->SetMarkerStyle(8);
   Graph_Graph03->SetMarkerSize(0.5);
   Graph_Graph03->GetXaxis()->SetTitle(" #sqrt{#it{s}} [GeV]  ");
   Graph_Graph03->GetXaxis()->SetNdivisions(505);
   Graph_Graph03->GetXaxis()->SetLabelFont(132);
   Graph_Graph03->GetXaxis()->SetLabelSize(0.065);
   Graph_Graph03->GetXaxis()->SetTitleSize(0.08);
   Graph_Graph03->GetXaxis()->SetTitleOffset(0.9);
   Graph_Graph03->GetXaxis()->SetTitleFont(132);
   Graph_Graph03->GetYaxis()->SetTitle(" #sqrt{#it{s}} / #it{m_{0} #it{#Gamma(s)}} [GeV]");
   Graph_Graph03->GetYaxis()->SetLabelFont(132);
   Graph_Graph03->GetYaxis()->SetLabelSize(0.065);
   Graph_Graph03->GetYaxis()->SetTitleSize(0.08);
   Graph_Graph03->GetYaxis()->SetTitleOffset(0.8);
   Graph_Graph03->GetYaxis()->SetTitleFont(132);
   Graph_Graph03->GetZaxis()->SetLabelFont(132);
   Graph_Graph03->GetZaxis()->SetLabelSize(0.06);
   Graph_Graph03->GetZaxis()->SetTitleSize(0.072);
   Graph_Graph03->GetZaxis()->SetTitleOffset(1.2);
   Graph_Graph03->GetZaxis()->SetTitleFont(132);
   graph->SetHistogram(Graph_Graph03);
   
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
   pt_LaTex = pt->AddText("K^{*}(1410)^{+}");
   pt->Draw();
   c1_n3->Modified();
   c1_n3->cd();
   c1_n3->SetSelected(c1_n3);
}
