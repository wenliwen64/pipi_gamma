// Basically this class is used to compute/correct gamma results 
class StLPVPlotMaker{
public:
    StLPVPlotMaker();
    ~StLPVPlotMaker();
    void Init(int , int );
    void Compute();
    void DrawSeparate();
    void DrawDiff();
private:
    double chi(double res);
    double resEventPlane(double chi);

    int mNCent;
    int mNJobs;
    double gamma[9];
    double ep_resolution[9];
    double mGamma_ss[9][10];
    double mGamma_ss_err[9][10];
    double mGamma_os[9][10];
    double mGamma_os_err[9][10];

    double mGamma_ss_final[9];
    double mGamma_ss_err_final[9];
    double mGamma_os_final[9];
    double mGamma_os_err_final[9];

    double mGamma_diff_final[9];
    double mGamma_diff_err_final[9];

    ClassDef(StLPVPlotMaker, 1)
};
