void Fig4_2()
{
//=========Macro generated from canvas: c1_n2/c1_n2
//=========  (Wed Mar 27 12:40:44 2024) by ROOT version 6.18/00
   TCanvas *c1_n2 = new TCanvas("c1_n2", "c1_n2",0,0,700,500);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1_n2->SetHighLightColor(2);
   c1_n2->Range(-0.3802469,-0.1357303,2.335802,0.7125839);
   c1_n2->SetFillColor(0);
   c1_n2->SetBorderMode(0);
   c1_n2->SetBorderSize(2);
   c1_n2->SetTickx(1);
   c1_n2->SetTicky(1);
   c1_n2->SetLeftMargin(0.14);
   c1_n2->SetRightMargin(0.05);
   c1_n2->SetTopMargin(0.05);
   c1_n2->SetBottomMargin(0.16);
   c1_n2->SetFrameLineWidth(2);
   c1_n2->SetFrameBorderMode(0);
   c1_n2->SetFrameLineWidth(2);
   c1_n2->SetFrameBorderMode(0);
   
   Double_t Graph0_fx2[40] = {
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
   Double_t Graph0_fy2[40] = {
   0,
   0,
   0,
   0,
   0,
   0,
   1.852055e-07,
   1.020572e-05,
   5.45448e-05,
   0.0002177844,
   0.0008855009,
   0.004525432,
   0.01403405,
   0.02838084,
   0.04646403,
   0.06749462,
   0.09125773,
   0.1169242,
   0.1428769,
   0.1689879,
   0.1950406,
   0.2208696,
   0.2463538,
   0.2714084,
   0.2959787,
   0.3200348,
   0.3435676,
   0.3665851,
   0.3891094,
   0.4111733,
   0.4328157,
   0.4540742,
   0.4749751,
   0.495524,
   0.5156998,
   0.5354594,
   0.5547506,
   0.573526,
   0.5917525,
   0.6092438};
   TGraph *graph = new TGraph(40,Graph0_fx2,Graph0_fy2);
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
   
   TH1F *Graph_Graph02 = new TH1F("Graph_Graph02","",100,0,2.2);
   Graph_Graph02->SetMinimum(0);
   Graph_Graph02->SetMaximum(0.6701682);
   Graph_Graph02->SetDirectory(0);
   Graph_Graph02->SetStats(0);
   Graph_Graph02->SetLineWidth(2);
   Graph_Graph02->SetMarkerStyle(8);
   Graph_Graph02->SetMarkerSize(0.5);
   Graph_Graph02->GetXaxis()->SetTitle(" #sqrt{#it{s}} [GeV]  ");
   Graph_Graph02->GetXaxis()->SetNdivisions(505);
   Graph_Graph02->GetXaxis()->SetLabelFont(132);
   Graph_Graph02->GetXaxis()->SetLabelSize(0.065);
   Graph_Graph02->GetXaxis()->SetTitleSize(0.08);
   Graph_Graph02->GetXaxis()->SetTitleOffset(0.9);
   Graph_Graph02->GetXaxis()->SetTitleFont(132);
   Graph_Graph02->GetYaxis()->SetTitle(" #sqrt{#it{s}} / #it{m_{0} #it{#Gamma(s)}} [GeV]");
   Graph_Graph02->GetYaxis()->SetLabelFont(132);
   Graph_Graph02->GetYaxis()->SetLabelSize(0.065);
   Graph_Graph02->GetYaxis()->SetTitleSize(0.08);
   Graph_Graph02->GetYaxis()->SetTitleOffset(0.8);
   Graph_Graph02->GetYaxis()->SetTitleFont(132);
   Graph_Graph02->GetZaxis()->SetLabelFont(132);
   Graph_Graph02->GetZaxis()->SetLabelSize(0.06);
   Graph_Graph02->GetZaxis()->SetTitleSize(0.072);
   Graph_Graph02->GetZaxis()->SetTitleOffset(1.2);
   Graph_Graph02->GetZaxis()->SetTitleFont(132);
   graph->SetHistogram(Graph_Graph02);
   
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
   pt_LaTex = pt->AddText("K_{1}(1400)^{+}");
   pt->Draw();
   c1_n2->Modified();
   c1_n2->cd();
   c1_n2->SetSelected(c1_n2);
}
