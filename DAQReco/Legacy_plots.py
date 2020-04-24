import os
import ROOT
from ROOT import *

def make_legend():
    output = ROOT.TLegend(0.44,0.06, 0.88,0.21, "", "brNDC")
    output.SetLineWidth(0)
    output.SetLineStyle(0)
    output.SetFillStyle(0)
    output.SetBorderSize(0)
    output.SetTextFont(62)
    output.SetTextSize(.045)
    return output

def add_CMS():
    lowX=0.10
    lowY=0.82
    lumi  = ROOT.TPaveText(lowX, lowY+0.06, lowX+0.15, lowY+0.16, "NDC")
    lumi.SetTextFont(61)
    lumi.SetTextSize(0.035)
    lumi.SetBorderSize(   0 )
    lumi.SetFillStyle(    0 )
    lumi.SetTextAlign(   12 )
    lumi.SetTextColor(    1 )
    lumi.AddText("FNAL June 2019 BTL Test Beam")
    return lumi


def main(args):

    InputFile=TFile(args.file,"read");
    data =InputFile.Get("data");

    c=ROOT.TCanvas("canvas","",0,0,600,600)
    c.cd()
    gStyle.SetOptStat(0)

    

#    X=data.Draw("gaus_mean[2]-gaus_mean[0] >> h(100,2,3)","(gaus_mean[2]-gaus_mean[0]) < 3 && (gaus_mean[2]-gaus_mean[0])> 2 && amp[2] > 100 && amp[0]> 100")

#    data->Draw("(chTime[15]-chTime[384])*1e-3 - (gaus_mean[2]-IL_50[8])","chTime[15]>0 && gaus_mean[2]>0 && (chTime[15]-chTime[384])*1e-3 - (gaus_mean[2]-IL_50[8])<-40&&(chTime[15]-chTime[384])*1e-3 - (gaus_mean[2]-IL_50[8]) > -60 && amp[2]> 50")
#    data->Draw("(chTime[15]-chTime[384])*1e-3 - (gaus_mean[2]-IL_50[8])")

    h1 = ROOT.TH1F("h1", "h1", 100, -60, -40)

    h1.GetXaxis().SetNdivisions(505)
    h1.GetYaxis().SetLabelFont(42)
    h1.GetYaxis().SetLabelOffset(0.01)
    h1.GetYaxis().SetLabelSize(0.04)
    h1.GetYaxis().SetTitleSize(0.05)
    h1.GetYaxis().SetTitleOffset(1.04)
    h1.SetTitle("")
#((chTime[6]-chTime[384])*1e-3) - ( gaus_mean[0]-IL_50[8])","chTime[6]>0 &&( ((chTime[6]-chTime[384])*1e-3) - ( gaus_mean[0]-IL_50[8]))> 970

#    data.Draw("((chTime[6]-chTime[384])*1e-3) - ( gaus_mean[0]-IL_50[8])","chTime[6]>0 &&( ((chTime[6]-chTime[384])*1e-3) - ( gaus_mean[0]-IL_50[8]))> 970")
    selection="((chTime[8]-chTime[384])*1e-3) - ( gaus_mean[2]-IL_50[8])"
    option= "chTime[8]>0 && gaus_mean[2] > 0"
#    selection="((chTime[6]-chTime[384])*1e-3) - ( gaus_mean[0]-IL_50[8])"

#    selection="(chTime[15]-chTime[384])*1e-3 - (gaus_mean[2]-IL_50[8])"
#    option= "chTime[15]>0 && gaus_mean[2]>0 && (chTime[15]-chTime[384])*1e-3 - (gaus_mean[2]-IL_50[8])<-40 && (chTime[15]-chTime[384])*1e-3 - (gaus_mean[2]-IL_50[8]) > -60 && amp[2]> 10"

    data.Project('h1',selection,option)
    

#    data.Fit('gaus',"(chTime[1]-chTime[384])*1e-3 - (gaus_mean[2]-IL_50[8])","chTime[1]>0 && gaus_mean[2]>0 && (chTime[1]-chTime[384])*1e-3 - (gaus_mean[2]-IL_50[8])<-40&&(chTime[1]-chTime[384])*1e-3 - (gaus_mean[2]-IL_50[8]) > -60 && amp[2]> 50")


#    f1 =TF1("f1", "gaus",-54,-50)
#    f1 =TF1("f1", "gaus")
#    print "\n name=",  htemp.GetName()
#    htemp.Fit("f1", '','')
#    f1.SetRange(2.6,2.61)
#    data.Fit("f1", '','')
    theFit=TF1("theFit","gaus",-50,-47)
#    h1.Fit("gaus", "", "", -54, -51.5)
    h1.Fit("theFit", "R0")
    h1.Draw()
    theFit.Draw("SAME")
    h1.GetXaxis().SetTitle(selection)
    FitParam=theFit.GetParameters()
    
    lowX=0.5
    lowY=0.6
    lumi  = ROOT.TPaveText(lowX, lowY+0.06, lowX+0.30, lowY+0.16, "NDC")
    lumi.AddText("#sigma="+str("%2.0f"%(FitParam[2]*1000))+" ps")
    lumi.Draw("same")
    
    l2=add_CMS()
    l2.Draw("same")
    
    c.SaveAs("test%s_%s.pdf"%(args.file.replace('ROOT_v4/','').replace('.root',''),args.type))


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--inputFile', '-i', action='store', dest='file', default='ROOT_v4/Tot_154.root', help='input file')
    parser.add_argument('--type', '-t', action='store', dest='type', default='', help='type of plot')
    main(parser.parse_args())

