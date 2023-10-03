#! /usr/bin/env python

import ROOT
import sys
from os import path
from DataFormats.FWLite import Events, Handle

hist = ROOT.TH1F('energy', ';HCAL energy;number of muons', 50, 0, 5)

aachen = '/net/data_cms/cmspos/HCAL/cosmic/525486FA-F98F-0849-843E-1E23C444DCB2.root'
lxplus = '/eos/cms/store/group/dpg_hcal/comm_hcal/CMS_POS/cosmic/525486FA-F98F-0849-843E-1E23C444DCB2.root'

if path.exists(aachen):
    inputfile = aachen
elif path.exists(lxplus):
    inputfile = lxplus
else:
    print('No input file found!')

events = Events(inputfile)
muonhandle = Handle('vector<reco::Muon>')
muonlabel = 'muons'

for event in events:
    event.getByLabel(muonlabel, muonhandle)
    muons = muonhandle.product()
    
    for muon in muons:
        depths = [] # use depths.append(x) to add an item to the list
        
        # Can you figure out how to extract calorimeter energy and depth information from the reco::Muon class?
        # Check the CMSSW doxygen for all its members/functions: http://cmsdoxygen.web.cern.ch/cmsdoxygen/CMSSW_10_6_1_patch2/doc/html/df/de3/classreco_1_1Muon.html
        
        if all(d in depths for d in [1,2,3,4]):
            hist.Fill(hadEnergy)

ROOT.gStyle.SetOptFit(1)
c = ROOT.TCanvas('c', 'c', 500, 450)
c.cd()
hist.Draw()
hist.Fit('landau')
c.Print('cosmics.pdf')
