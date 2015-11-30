void run(){
    gSystem->Load("./StLPVPlotMaker_cc.so");
    StLPVPlotMaker* maker = new StLPVPlotMaker();

    maker->Init(4, 10);
    maker->Compute();
    maker->DrawSeparate();
}
