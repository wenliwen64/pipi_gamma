#include <stdio.h>
#include <TFile.h>
#include <fstream>
#include <iostream>
#include "TTree.h"
#include "TMath.h"
#include "TVector2.h"
#include "./StThreeVectorD.h"
#include "./StLorentzVectorD.h"
#include "./SystemOfUnits.h"

const float PI = TMath::Pi();
const float c_light = 2.99792458e+8 * meter/second;

 struct Helice {
    StThreeVectorD p;
    StThreeVectorD o;
    double B;
    double q;
    double pha;
    double phase();
    double dipAngle() {return atan2(p.z(),p.perp());}
    double curvature() {return fabs((c_light*nanosecond/meter*q*B/tesla)/(abs(p)/GeV*cos(dipAngle()))/meter);}
    double xcenter() {return o.x()-cos(phase())/curvature();}
    double ycenter() {return o.y()-sin(phase())/curvature();}
    double x(double s) { int mH = (q*B <= 0) ? 1 : -1; return o.x() + (cos(phase() + s*mH*curvature()*cos(dipAngle()))-cos(phase()))/curvature(); }
    double y(double s) { int mH = (q*B <= 0) ? 1 : -1; return o.y() + (sin(phase() + s*mH*curvature()*cos(dipAngle()))-sin(phase()))/curvature(); }
    double z(double s) { return o.z() + s*sin(dipAngle());}
    double period() { int mH = (q*B <= 0) ? 1 : -1; return fabs(2*PI/(mH*curvature()*cos(dipAngle())));}
    StThreeVector<double> at(double s) { return StThreeVector<double>(x(s), y(s), z(s));}
    double fudgePathLength(StThreeVector<double>& pp);
    double      pathLength(StThreeVector<double>& pp);
    double distance(StThreeVector<double>& pp) { return abs(this->at(pathLength(pp))-pp); }
    double pathLength(double x, double y) { StThreeVectorD temp(x, y, 0); return fudgePathLength(temp); }
    pair<double, double> pathLengths(struct Helice& h);
    StThreeVector<double> momentum(double BB) { int mH = (q*B <= 0) ? 1 : -1; double pt = GeV*fabs(c_light*nanosecond/meter*BB/tesla)/(fabs(curvature())*meter);
	return (StThreeVector<double>(pt*cos(phase()+mH*PI/2),pt*sin(phase()+mH*PI/2), pt*tan(dipAngle())));}
    void moveOrigin(double s) { StThreeVector<double> newOrigin = at(s); double newPhase = atan2(newOrigin.y()-ycenter(), newOrigin.x()-xcenter()); o = newOrigin; pha = newPhase;}
    StThreeVector<double> momentumAt(double S, double BB);// { struct Helice tmp(*this); tmp.moveOrigin(S); return tmp.momentum(BB); }
    
  };

StThreeVector<double> Helice::momentumAt(double S, double BB) {
	struct Helice tmp(*this);
	StThreeVector<double> newOrigin = tmp.at(S);
	double newPhase = atan2(newOrigin.y()-tmp.ycenter(), newOrigin.x()-tmp.xcenter());
	tmp.o = newOrigin; 
	tmp.pha = newPhase;
	int mH = (q*B <= 0) ? 1 : -1;
        double pt = GeV*fabs(c_light*nanosecond/meter*BB/tesla)/(fabs(curvature())*meter);
	return (StThreeVector<double>(pt*cos(newPhase+mH*PI/2),pt*sin(newPhase+mH*PI/2), pt*tan(dipAngle())));	
//	return tmp.momentum(BB);
}

double Helice::phase() {
        int mH = (q*B <= 0) ? 1 : -1;
        if(p.y() == 0 && p.x() == 0) return (PI/4)*(1-2.*mH);
        else return atan2(p.y(),p.x())-mH*PI/2;
    }

double Helice::fudgePathLength(StThreeVector<double>& pp) {
    int mH = (q*B <= 0) ? 1 : -1;
    double dx = pp.x()-o.x();
    double dy = pp.y()-o.y();
    double s = atan2(dy*cos(phase()) - dx*sin(phase()), 1/curvature() + dx*cos(phase())+dy*sin(phase()))/ (mH*curvature()*cos(dipAngle()));
    return s;
}

double Helice::     pathLength(StThreeVector<double>& pp) {
    int mH = (q*B <= 0) ? 1 : -1;
    double s;
    double dx = pp.x()-o.x();
    double dy = pp.y()-o.y();
    double dz = pp.z()-o.z();
            const double MaxPrecisionNeeded = micrometer;
            const int    MaxIterations      = 100;

            double t34 = curvature()*cos(dipAngle())*cos(dipAngle());
            double t41 = sin(dipAngle())*sin(dipAngle());
            double t6, t7, t11, t12, t19;

	    s = fudgePathLength(pp);
    	        double ds = period();
                int    j, jmin = 0;
                double d, dmin = abs(at(s) - pp);
                for(j=1; j<MaxIterations; j++) {
                  if ((d = abs(at(s+j*ds) - p)) < dmin) {
                      dmin = d;
                      jmin = j;
                  }
                  else
                      break;
                }
                for(j=-1; -j<MaxIterations; j--) {
                  if ((d = abs(at(s+j*ds) - p)) < dmin) {
                      dmin = d;
                      jmin = j;
                  }
                  else
                      break;
                }
                if (jmin) s += jmin*ds;

            double sOld = s;
            for (int i=0; i<MaxIterations; i++) {
                t6  = phase()+s*mH*curvature()*cos(dipAngle());
                t7  = cos(t6);
                t11 = dx-(1/curvature())*(t7-cos(phase()));
                t12 = sin(t6);
                t19 = dy-(1/curvature())*(t12-sin(phase()));
                s  -= (t11*t12*mH*cos(dipAngle())-t19*t7*mH*cos(dipAngle()) - (dz-s*sin(dipAngle()))*sin(dipAngle()))/
                    (t12*t12*cos(dipAngle())*cos(dipAngle())+t11*t7*t34 + t7*t7*cos(dipAngle())*cos(dipAngle()) + t19*t12*t34+t41);
                if (fabs(sOld-s) < MaxPrecisionNeeded) break;
                sOld = s;
            }

    return s;
}

//double Helice::pathLength(double x, double y) { return fudgePathLength(StThreeVector<double>(x, y, 0)); }

pair<double, double> Helice::pathLengths(struct Helice& h) {
    double s1, s2;
        double dx = h.xcenter() - xcenter();
        double dy = h.ycenter() - ycenter();
        double dd = ::sqrt(dx*dx + dy*dy);
        double r1 = 1/curvature();
        double r2 = 1/h.curvature();
        double cosAlpha = (r1*r1 + dd*dd - r2*r2)/(2*r1*dd);

        double s;
        double x, y;
        if (fabs(cosAlpha) < 1) {
            double sinAlpha = sin(acos(cosAlpha));
            x = xcenter() + r1*(cosAlpha*dx - sinAlpha*dy)/dd;
            y = ycenter() + r1*(sinAlpha*dx + cosAlpha*dy)/dd;
            s = pathLength(x, y);
            x = xcenter() + r1*(cosAlpha*dx + sinAlpha*dy)/dd;
            y = ycenter() + r1*(cosAlpha*dy - sinAlpha*dx)/dd;
            double a = pathLength(x, y);
	    StThreeVectorD temp_a = at(a);
	    StThreeVectorD temp_s = at(s);
            if (h.distance(temp_a) < h.distance(temp_s)) s = a;
        }
        else { 
            int rsign = ((r2-r1) > dd ? -1 : 1);
            x = xcenter() + rsign*r1*dx/dd;
            y = ycenter() + rsign*r1*dy/dd;
            s = pathLength(x, y);
        }
        const double MinStepSize = 10*micrometer;
        const double MinRange    = 10*centimeter;
        StThreeVectorD temp_s = at(s);
        double dmin              = h.distance(temp_s);
        double range             = max(2*dmin, MinRange);
        double ds                = range/10;
        double slast=-999999, ss, d;
        s1 = s - range/2.;
        s2 = s + range/2.;

        while (ds > MinStepSize) {
            for (ss=s1; ss<s2+ds; ss+=ds) {
		StThreeVectorD temp_ss = at(ss);
                d = h.distance(temp_ss);
                if (d < dmin) {
                    dmin = d;
                    s = ss;
                }
                slast = ss;
            }
            if (s == s1) {
                d = 0.8*(s2-s1);
                s1 -= d;
                s2 -= d;
            }
            else if (s == slast) {
                d = 0.8*(s2-s1);
                s1 += d;
                s2 += d;
            }
            else {
                s1 = s-ds;
                s2 = s+ds;
                ds /= 10;
            }
        }

        StThreeVectorD temp = at(s);
    return pair<double, double>(s, h.pathLength(temp));
}
