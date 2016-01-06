void run(){
    //gSystem->Load("./StLPVPlotMaker_cc.so");
    gROOT->LoadMacro("./StLPVPlotMaker.cc++");
    StLPVPlotMaker* maker = new StLPVPlotMaker();

    maker->Init(7, 10);
    maker->Compute();
    maker->DrawSeparate();
    maker->DrawDiff();
}
