void sign_window(TVirtualPad *p, TH2D* h, TString xaxis, TString yaxis, TString title)
{
  Double_t height = 1-p->GetTopMargin()-p->GetBottomMargin();
  Double_t width = 1-p->GetLeftMargin()-p->GetRightMargin();
  //  cout << "pad width is: " << width << endl;
  //  cout << "pad height is: " << height << endl;
  if((width == 0.) || (height == 0.)) {
    cout << "In sign_window(). Width or height of pad is 0. Cannot set text size." << endl;
    exit(-1);
  }
  Double_t fontsizex = 0.04/width; //0.06
  Double_t fontsizey = 0.04/height;
  h->SetStats(kFALSE);
  h->SetTitle(title);
  h->GetXaxis()->SetTitle(xaxis);
  h->GetYaxis()->SetTitle(yaxis);
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
//  h->GetXaxis()->SetTitleSize(fontsizex);
  h->GetXaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleSize(0.65*fontsizey);
  h->GetYaxis()->SetTitleFont(42);
//  h->GetXaxis()->SetLabelSize(fontsizex);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetTitleOffset(1.2);
//  h->GetYaxis()->SetLabelSize(fontsizey);
  h->GetYaxis()->SetLabelFont(42);
  if(p->GetLogy())
    h->GetYaxis()->SetTitleOffset(2.);
  else
    h->GetYaxis()->SetTitleOffset(2.);
  h->GetXaxis()->SetNdivisions(511, kTRUE);
  h->GetYaxis()->SetNdivisions(507, kTRUE);
}


void make_clean_canv(TCanvas *c)
{
  c->SetBorderMode(0);
  c->SetFillColor(kWhite);
  c->SetFillStyle(0);

}
void make_clean_pads(TCanvas *c, Int_t n, Int_t gridx, Int_t gridy)
{
  c->SetBorderMode(0);
  c->SetFillColor(kWhite);
  c->SetFillStyle(0);
 
  for(Int_t i=0; i<n; i++)
    {
      c->GetPad(i+1)->SetBottomMargin(0.2);
      c->GetPad(i+1)->SetTopMargin(0.05);
      c->GetPad(i+1)->SetLeftMargin(0.2);
      c->GetPad(i+1)->SetRightMargin(0.05);
      c->GetPad(i+1)->SetFrameBorderMode(0);
      c->GetPad(i+1)->SetBorderMode(0);
      c->GetPad(i+1)->SetFillColor(kWhite);
      c->GetPad(i+1)->SetGrid(gridx, gridy);
    }
}

void make_clean_pads(TVirtualPad *c, Int_t n, Int_t gridx, Int_t gridy)
{

  for(Int_t i=0; i<n; i++)
    {
      //          c->GetPad(i+1)->SetBottomMargin(0.17);
      //            c->GetPad(i+1)->SetTopMargin(0.05);
      //            c->GetPad(i+1)->SetLeftMargin(0.17);
      //            c->GetPad(i+1)->SetRightMargin(0.05);
      c->GetPad(i+1)->SetTicks(1,1);
      c->GetPad(i+1)->SetFrameBorderMode(0);
      c->GetPad(i+1)->SetBorderMode(0);
      c->GetPad(i+1)->SetFillColor(kWhite);
      c->GetPad(i+1)->SetGrid(gridx, gridy);
    }
}
