void Fig4_8()
{
//=========Macro generated from canvas: c1_n8/c1_n8
//=========  (Tue Mar 26 22:55:10 2024) by ROOT version 6.18/00
   TCanvas *c1_n8 = new TCanvas("c1_n8", "c1_n8",0,0,700,500);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1_n8->SetHighLightColor(2);
   c1_n8->Range(2.684548,-0.1852393,6.441434,0.9725064);
   c1_n8->SetFillColor(0);
   c1_n8->SetBorderMode(0);
   c1_n8->SetBorderSize(2);
   c1_n8->SetTickx(1);
   c1_n8->SetTicky(1);
   c1_n8->SetLeftMargin(0.14);
   c1_n8->SetRightMargin(0.05);
   c1_n8->SetTopMargin(0.05);
   c1_n8->SetBottomMargin(0.16);
   c1_n8->SetFrameLineWidth(2);
   c1_n8->SetFrameBorderMode(0);
   c1_n8->SetFrameLineWidth(2);
   c1_n8->SetFrameBorderMode(0);
   
   Double_t Graph0_fx8[40] = {
   3.464102,
   3.551814,
   3.637412,
   3.721042,
   3.802833,
   3.882901,
   3.961352,
   4.038278,
   4.113767,
   4.187895,
   4.260733,
   4.332347,
   4.402796,
   4.472136,
   4.540417,
   4.607686,
   4.673987,
   4.739361,
   4.803845,
   4.867474,
   4.930283,
   4.992302,
   5.053559,
   5.114083,
   5.173899,
   5.233031,
   5.291503,
   5.349335,
   5.406549,
   5.463163,
   5.519197,
   5.574668,
   5.629592,
   5.683986,
   5.737863,
   5.79124,
   5.844129,
   5.896544,
   5.948497,
   6};
   Double_t Graph0_fy8[40] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   4.348895e-05,
   0.001848964,
   0.006371084,
   0.01363598,
   0.02361551,
   0.03624706,
   0.05144499,
   0.06911035,
   0.08913631,
   0.1114122,
   0.1358266,
   0.1622687,
   0.1906298,
   0.2208048,
   0.2526922,
   0.2861947,
   0.3212195,
   0.3576777,
   0.3954852,
   0.4345623,
   0.4748331,
   0.5162265,
   0.5586745,
   0.6021133,
   0.6464825,
   0.6917253,
   0.7377879,
   0.7846196,
   0.8314719};
   TGraph *graph = new TGraph(40,Graph0_fx8,Graph0_fy8);
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
   
   TH1F *Graph_Graph08 = new TH1F("Graph_Graph08","",100,3.210512,6.25359);
   Graph_Graph08->SetMinimum(0);
   Graph_Graph08->SetMaximum(0.9146191);
   Graph_Graph08->SetDirectory(0);
   Graph_Graph08->SetStats(0);
   Graph_Graph08->SetLineWidth(2);
   Graph_Graph08->SetMarkerStyle(8);
   Graph_Graph08->SetMarkerSize(0.5);
   Graph_Graph08->GetXaxis()->SetTitle(" #sqrt{#it{s}} [GeV]  ");
   Graph_Graph08->GetXaxis()->SetNdivisions(505);
   Graph_Graph08->GetXaxis()->SetLabelFont(132);
   Graph_Graph08->GetXaxis()->SetLabelSize(0.065);
   Graph_Graph08->GetXaxis()->SetTitleSize(0.08);
   Graph_Graph08->GetXaxis()->SetTitleOffset(0.9);
   Graph_Graph08->GetXaxis()->SetTitleFont(132);
   Graph_Graph08->GetYaxis()->SetTitle(" #sqrt{#it{s}} / #it{m_{0} #it{#Gamma(s)}} [GeV]");
   Graph_Graph08->GetYaxis()->SetLabelFont(132);
   Graph_Graph08->GetYaxis()->SetLabelSize(0.065);
   Graph_Graph08->GetYaxis()->SetTitleSize(0.08);
   Graph_Graph08->GetYaxis()->SetTitleOffset(0.8);
   Graph_Graph08->GetYaxis()->SetTitleFont(132);
   Graph_Graph08->GetZaxis()->SetLabelFont(132);
   Graph_Graph08->GetZaxis()->SetLabelSize(0.06);
   Graph_Graph08->GetZaxis()->SetTitleSize(0.072);
   Graph_Graph08->GetZaxis()->SetTitleOffset(1.2);
   Graph_Graph08->GetZaxis()->SetTitleFont(132);
   graph->SetHistogram(Graph_Graph08);
   
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
   pt_LaTex = pt->AddText("X_{s}^{0}");
   pt->Draw();
   c1_n8->Modified();
   c1_n8->cd();
   c1_n8->SetSelected(c1_n8);
}
