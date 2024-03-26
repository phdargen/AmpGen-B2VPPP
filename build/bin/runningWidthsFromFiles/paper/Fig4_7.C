void Fig4_7()
{
//=========Macro generated from canvas: c1_n7/c1_n7
//=========  (Tue Mar 26 22:55:10 2024) by ROOT version 6.18/00
   TCanvas *c1_n7 = new TCanvas("c1_n7", "c1_n7",0,0,700,500);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1_n7->SetHighLightColor(2);
   c1_n7->Range(2.991955,-0.2494333,5.26736,1.309525);
   c1_n7->SetFillColor(0);
   c1_n7->SetBorderMode(0);
   c1_n7->SetBorderSize(2);
   c1_n7->SetTickx(1);
   c1_n7->SetTicky(1);
   c1_n7->SetLeftMargin(0.14);
   c1_n7->SetRightMargin(0.05);
   c1_n7->SetTopMargin(0.05);
   c1_n7->SetBottomMargin(0.16);
   c1_n7->SetFrameLineWidth(2);
   c1_n7->SetFrameBorderMode(0);
   c1_n7->SetFrameLineWidth(2);
   c1_n7->SetFrameBorderMode(0);
   
   Double_t Graph0_fx7[40] = {
   3.464102,
   3.511885,
   3.559026,
   3.605551,
   3.651484,
   3.696846,
   3.741657,
   3.785939,
   3.829708,
   3.872983,
   3.91578,
   3.958114,
   4,
   4.041452,
   4.082483,
   4.123106,
   4.163332,
   4.203173,
   4.242641,
   4.281744,
   4.320494,
   4.358899,
   4.396969,
   4.434712,
   4.472136,
   4.50925,
   4.546061,
   4.582576,
   4.618802,
   4.654747,
   4.690416,
   4.725816,
   4.760952,
   4.795832,
   4.830459,
   4.86484,
   4.898979,
   4.932883,
   4.966555,
   5};
   Double_t Graph0_fy7[40] = {
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
   0,
   0.0008656648,
   0.004332578,
   0.01062521,
   0.01987984,
   0.03219031,
   0.04761851,
   0.06619865,
   0.08794235,
   0.112845,
   0.1408867,
   0.1720354,
   0.206251,
   0.243485,
   0.2836825,
   0.326785,
   0.3727299,
   0.4214512,
   0.4728822,
   0.5269545,
   0.5835979,
   0.6427435,
   0.7043214,
   0.768262,
   0.8344969,
   0.9029582,
   0.9735785,
   1.046292,
   1.119615};
   TGraph *graph = new TGraph(40,Graph0_fx7,Graph0_fy7);
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
   
   TH1F *Graph_Graph07 = new TH1F("Graph_Graph07","",100,3.310512,5.15359);
   Graph_Graph07->SetMinimum(0);
   Graph_Graph07->SetMaximum(1.231577);
   Graph_Graph07->SetDirectory(0);
   Graph_Graph07->SetStats(0);
   Graph_Graph07->SetLineWidth(2);
   Graph_Graph07->SetMarkerStyle(8);
   Graph_Graph07->SetMarkerSize(0.5);
   Graph_Graph07->GetXaxis()->SetTitle(" #sqrt{#it{s}} [GeV]  ");
   Graph_Graph07->GetXaxis()->SetNdivisions(505);
   Graph_Graph07->GetXaxis()->SetLabelFont(132);
   Graph_Graph07->GetXaxis()->SetLabelSize(0.065);
   Graph_Graph07->GetXaxis()->SetTitleSize(0.08);
   Graph_Graph07->GetXaxis()->SetTitleOffset(0.9);
   Graph_Graph07->GetXaxis()->SetTitleFont(132);
   Graph_Graph07->GetYaxis()->SetTitle(" #sqrt{#it{s}} / #it{m_{0} #it{#Gamma(s)}} [GeV]");
   Graph_Graph07->GetYaxis()->SetLabelFont(132);
   Graph_Graph07->GetYaxis()->SetLabelSize(0.065);
   Graph_Graph07->GetYaxis()->SetTitleSize(0.08);
   Graph_Graph07->GetYaxis()->SetTitleOffset(0.8);
   Graph_Graph07->GetYaxis()->SetTitleFont(132);
   Graph_Graph07->GetZaxis()->SetLabelFont(132);
   Graph_Graph07->GetZaxis()->SetLabelSize(0.06);
   Graph_Graph07->GetZaxis()->SetTitleSize(0.072);
   Graph_Graph07->GetZaxis()->SetTitleOffset(1.2);
   Graph_Graph07->GetZaxis()->SetTitleFont(132);
   graph->SetHistogram(Graph_Graph07);
   
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
   pt_LaTex = pt->AddText("X^{0}");
   pt->Draw();
   c1_n7->Modified();
   c1_n7->cd();
   c1_n7->SetSelected(c1_n7);
}
