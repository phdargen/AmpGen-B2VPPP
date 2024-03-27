void Fig4_8()
{
//=========Macro generated from canvas: c1_n8/c1_n8
//=========  (Wed Mar 27 12:40:45 2024) by ROOT version 6.18/00
   TCanvas *c1_n8 = new TCanvas("c1_n8", "c1_n8",0,0,700,500);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1_n8->SetHighLightColor(2);
   c1_n8->Range(3.385185,-0.101022,6.348148,0.5303656);
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
   4,
   4.063597,
   4.126214,
   4.187895,
   4.24868,
   4.308608,
   4.367714,
   4.42603,
   4.483588,
   4.540417,
   4.596543,
   4.651992,
   4.706787,
   4.760952,
   4.814508,
   4.867474,
   4.919871,
   4.971715,
   5.023024,
   5.073814,
   5.124101,
   5.173899,
   5.223222,
   5.272084,
   5.320497,
   5.368474,
   5.416026,
   5.463163,
   5.509898,
   5.556239,
   5.602197,
   5.647782,
   5.693001,
   5.737863,
   5.782378,
   5.826553,
   5.870395,
   5.913912,
   5.957112,
   6};
   Double_t Graph0_fy8[40] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0.0003332904,
   0.001664266,
   0.004031093,
   0.007436512,
   0.01187044,
   0.01731347,
   0.02374014,
   0.0311206,
   0.03942239,
   0.04861133,
   0.0586524,
   0.06951035,
   0.08115006,
   0.09353701,
   0.1066373,
   0.120418,
   0.1348472,
   0.1498941,
   0.1655289,
   0.1817232,
   0.1984494,
   0.2156816,
   0.2333945,
   0.2515642,
   0.2701677,
   0.2891834,
   0.3085901,
   0.3283682,
   0.3484986,
   0.3689634,
   0.3897451,
   0.4108277,
   0.4321953,
   0.4534511};
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
   
   TH1F *Graph_Graph08 = new TH1F("Graph_Graph08","",100,3.8,6.2);
   Graph_Graph08->SetMinimum(0);
   Graph_Graph08->SetMaximum(0.4987962);
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
