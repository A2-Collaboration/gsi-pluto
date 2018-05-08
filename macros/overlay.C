void overlay(Char_t *image,
             Float_t xmin, Float_t xmax, Float_t ymin, Float_t ymax,
             Bool_t logX=kFALSE, Bool_t logY=kTRUE)
{
   // Display image and graph for picking on top.

   TImage *img = TImage::Open(image);
   if (!img) {
      printf("Could not read the image... exit\n");
      return;
   }
   img->SetConstRatio(kTRUE);

   TCanvas *c = new TCanvas("overlay", "overlay", 800, 800);
   img->Draw("X");

   gPad->Divide(1,1);
   gPad->cd(1);
   gPad->SetFillStyle(4000);
   Float_t x[2] = {xmin,xmax};
   Float_t y[2] = {ymin,ymax};
   TGraph *g = new TGraph(2,x,y);
   g->SetTitle(image);
   g->GetXaxis()->SetRangeUser(xmin,xmax);
   g->GetYaxis()->SetRangeUser(ymin,ymax);
   gPad->SetGridx();
   gPad->SetGridy();
   if (logX) gPad->SetLogx();
   if (logY) gPad->SetLogy();
   gPad->SetFrameLineColor(2);
   gPad->SetFrameLineWidth(2);
   gPad->SetFrameFillStyle(4000);
   g->Draw("ap");

   c->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)", 0, "",
               "execEvent(Int_t,Int_t,Int_t,TObject*)");

   cout << "Adjust red frame onto picture." << endl << endl;
   cout << "Press middle mouse button to start." << endl;
   cout << "Press repeatedly left button to mark data points." << endl;
   cout << "Press again middle button to finish list." << endl << endl;
}

void execEvent(Int_t event, Int_t ix, Int_t iy, TObject *selected)
{
   TCanvas *c = (TCanvas *) gTQSender;

   if (strcmp(c->GetName(),"overlay") == 0) { // inside proper canvas
     if (strcmp(selected->IsA()->GetName(),"TFrame") == 0) { // inside frame
       if (event==1) {
         Float_t x;
         Float_t y;
         Char_t *string = selected->GetObjectInfo(ix,iy);    // get x,y as a string
	 cout << string << endl;
//          Char_t *pos = strtok(string,"x=, y");
//          sscanf(pos, "%f", &x);
//          pos = strtok(NULL,"x=, y");
//          sscanf(pos, "%f", &y);
//          printf("%f   %f\n", x, y); // left mouse button pressed
       }
       else if (event==2) {
         if (gPad->GetCrosshair() == 0) gPad->SetCrosshair(2);
         else gPad->SetCrosshair(0);
         printf("\n    x        y\n\n", event); // middle mouse button pressed
       }
     }
   }
}
