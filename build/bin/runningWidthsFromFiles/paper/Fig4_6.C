void Fig4_6()
{
//=========Macro generated from canvas: c1_n6/c1_n6
//=========  (Tue Mar 26 22:55:10 2024) by ROOT version 6.18/00
   TCanvas *c1_n6 = new TCanvas("c1_n6", "c1_n6",0,0,700,500);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1_n6->SetHighLightColor(2);
   c1_n6->Range(-0.4657055,-0.0665317,2.860762,0.3492915);
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
   Double_t Graph0_fy6[40] = {
   0,
   0,
   0,
   0,
   6.725847e-05,
   0.004748207,
   0.0138996,
   0.02521312,
   0.03755889,
   0.05034732,
   0.0632569,
   0.07610766,
   0.08879813,
   0.1012719,
   0.1134995,
   0.1254672,
   0.1371712,
   0.1486137,
   0.1598005,
   0.1707394,
   0.1814393,
   0.1919098,
   0.2021603,
   0.2122006,
   0.2220397,
   0.2316867,
   0.24115,
   0.2504378,
   0.2595577,
   0.268517,
   0.2773224,
   0.2859806,
   0.2944974,
   0.2986366,
   0.2986366,
   0.2986366,
   0.2986366,
   0.2986366,
   0.2986366,
   0.2986366};
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
   
   TH1F *Graph_Graph06 = new TH1F("Graph_Graph06","",100,0,2.694439);
   Graph_Graph06->SetMinimum(0);
   Graph_Graph06->SetMaximum(0.3285003);
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
