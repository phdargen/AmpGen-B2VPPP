void Fig4_7()
{
//=========Macro generated from canvas: c1_n7/c1_n7
//=========  (Wed Mar 27 12:40:45 2024) by ROOT version 6.18/00
   TCanvas *c1_n7 = new TCanvas("c1_n7", "c1_n7",0,0,700,500);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1_n7->SetHighLightColor(2);
   c1_n7->Range(3.692593,-0.08150196,5.174074,0.4294825);
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
   4,
   4.028743,
   4.057282,
   4.085622,
   4.113767,
   4.14172,
   4.169486,
   4.197069,
   4.224471,
   4.251696,
   4.278749,
   4.305631,
   4.332347,
   4.358899,
   4.38529,
   4.411523,
   4.437602,
   4.463527,
   4.489304,
   4.514932,
   4.540417,
   4.565759,
   4.590961,
   4.616026,
   4.640955,
   4.665751,
   4.690416,
   4.714952,
   4.739361,
   4.763644,
   4.787805,
   4.811844,
   4.835764,
   4.859566,
   4.883252,
   4.906824,
   4.930283,
   4.953631,
   4.97687,
   5};
   Double_t Graph0_fy7[40] = {
   0.0002839412,
   0.0009749083,
   0.002102563,
   0.003683619,
   0.005732106,
   0.008258254,
   0.01127019,
   0.01477345,
   0.01877172,
   0.02326711,
   0.02825992,
   0.03374975,
   0.0397343,
   0.04621135,
   0.05317668,
   0.06062657,
   0.06855556,
   0.07695865,
   0.08582962,
   0.09516244,
   0.1049504,
   0.1151868,
   0.1258646,
   0.1369764,
   0.1485151,
   0.160473,
   0.172843,
   0.185617,
   0.1987879,
   0.2123478,
   0.2262893,
   0.2406047,
   0.2552867,
   0.2703276,
   0.2857201,
   0.301457,
   0.3175308,
   0.3339346,
   0.3506611,
   0.3672378};
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
   
   TH1F *Graph_Graph07 = new TH1F("Graph_Graph07","",100,3.9,5.1);
   Graph_Graph07->SetMinimum(0.0002555471);
   Graph_Graph07->SetMaximum(0.4039332);
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
