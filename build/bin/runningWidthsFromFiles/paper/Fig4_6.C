void Fig4_6()
{
//=========Macro generated from canvas: c1_n6/c1_n6
//=========  (Wed Mar 27 12:40:45 2024) by ROOT version 6.18/00
   TCanvas *c1_n6 = new TCanvas("c1_n6", "c1_n6",0,0,700,500);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1_n6->SetHighLightColor(2);
   c1_n6->Range(-0.3802469,-0.05372456,2.335802,0.282054);
   c1_n6->SetFillColor(0);
   c1_n6->SetBorderMode(0);
   c1_n6->SetBorderSize(2);
   c1_n6->SetTickx(1);
   c1_n6->SetTicky(1);
   c1_n6->SetLeftMargin(0.14);
   c1_n6->SetRightMargin(0.05);
   c1_n6->SetTopMargin(0.05);
   c1_n6->SetBottomMargin(0.16);
   c1_n6->SetFrameLineWidth(2);
   c1_n6->SetFrameBorderMode(0);
   c1_n6->SetFrameLineWidth(2);
   c1_n6->SetFrameBorderMode(0);
   
   Double_t Graph0_fx6[40] = {
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
   Double_t Graph0_fy6[40] = {
   0,
   0,
   0,
   0,
   0,
   0,
   6.725847e-05,
   0.002521044,
   0.007451535,
   0.0138996,
   0.02128643,
   0.02925038,
   0.03755889,
   0.04605893,
   0.05464807,
   0.0632569,
   0.07183804,
   0.08035912,
   0.08879813,
   0.09714032,
   0.1053761,
   0.1134995,
   0.1215072,
   0.1293978,
   0.1371712,
   0.1448283,
   0.1523707,
   0.1598005,
   0.1671201,
   0.1743321,
   0.1814393,
   0.1884445,
   0.1953506,
   0.2021603,
   0.2088766,
   0.2155022,
   0.2220397,
   0.2284919,
   0.2348612,
   0.24115};
   TGraph *graph = new TGraph(40,Graph0_fx6,Graph0_fy6);
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
   
   TH1F *Graph_Graph06 = new TH1F("Graph_Graph06","",100,0,2.2);
   Graph_Graph06->SetMinimum(0);
   Graph_Graph06->SetMaximum(0.265265);
   Graph_Graph06->SetDirectory(0);
   Graph_Graph06->SetStats(0);
   Graph_Graph06->SetLineWidth(2);
   Graph_Graph06->SetMarkerStyle(8);
   Graph_Graph06->SetMarkerSize(0.5);
   Graph_Graph06->GetXaxis()->SetTitle(" #sqrt{#it{s}} [GeV]  ");
   Graph_Graph06->GetXaxis()->SetNdivisions(505);
   Graph_Graph06->GetXaxis()->SetLabelFont(132);
   Graph_Graph06->GetXaxis()->SetLabelSize(0.065);
   Graph_Graph06->GetXaxis()->SetTitleSize(0.08);
   Graph_Graph06->GetXaxis()->SetTitleOffset(0.9);
   Graph_Graph06->GetXaxis()->SetTitleFont(132);
   Graph_Graph06->GetYaxis()->SetTitle(" #sqrt{#it{s}} / #it{m_{0} #it{#Gamma(s)}} [GeV]");
   Graph_Graph06->GetYaxis()->SetLabelFont(132);
   Graph_Graph06->GetYaxis()->SetLabelSize(0.065);
   Graph_Graph06->GetYaxis()->SetTitleSize(0.08);
   Graph_Graph06->GetYaxis()->SetTitleOffset(0.8);
   Graph_Graph06->GetYaxis()->SetTitleFont(132);
   Graph_Graph06->GetZaxis()->SetLabelFont(132);
   Graph_Graph06->GetZaxis()->SetLabelSize(0.06);
   Graph_Graph06->GetZaxis()->SetTitleSize(0.072);
   Graph_Graph06->GetZaxis()->SetTitleOffset(1.2);
   Graph_Graph06->GetZaxis()->SetTitleFont(132);
   graph->SetHistogram(Graph_Graph06);
   
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
   pt_LaTex = pt->AddText("K_{2}(1770)^{+}");
   pt->Draw();
   c1_n6->Modified();
   c1_n6->cd();
   c1_n6->SetSelected(c1_n6);
}
