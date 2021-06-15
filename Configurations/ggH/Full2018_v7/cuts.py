# cuts

supercut = ' mll > 12 \
            && Lepton_pt[0]>25 \
            && Lepton_pt[1]>13 \
            && (nLepton>=2 && Alt$(Lepton_pt[2],0)<10) \
            && abs(Lepton_eta[0])<2.5 && abs(Lepton_eta[1])<2.5 \
            && ptll>30 \
            && PuppiMET_pt > 20 \
            && (Lepton_pdgId[0]*Lepton_pdgId[1] == -11*13) \
           '
'''
## Signal regions
cuts['hww2l2v_13TeV'] = {
   'expr': 'sr',
    # Define the sub-categorization of sr
   'categories' : {
      'em_pm_0j_pt2ge20' : ' Lepton_pdgId[0]==-11 && Lepton_pt[1]>=20 && zeroJet',
      'em_mp_0j_pt2ge20' : ' Lepton_pdgId[0]==11 && Lepton_pt[1]>=20 && zeroJet',
      'me_mp_0j_pt2ge20' : ' Lepton_pdgId[0]==-13 && Lepton_pt[1]>=20 && zeroJet',
      'me_pm_0j_pt2ge20' : ' Lepton_pdgId[0]==13 && Lepton_pt[1]>=20 && zeroJet',
      #
      'em_pm_0j_pt2lt20' : ' Lepton_pdgId[0]==-11 && Lepton_pt[1]<20 && zeroJet',
      'em_mp_0j_pt2lt20' : ' Lepton_pdgId[0]==11 && Lepton_pt[1]<20 && zeroJet',
      'me_mp_0j_pt2lt20' : ' Lepton_pdgId[0]==-13 && Lepton_pt[1]<20 && zeroJet',
      'me_pm_0j_pt2lt20' : ' Lepton_pdgId[0]==13 && Lepton_pt[1]<20 && zeroJet',
      #
      'em_pm_1j_pt2ge20' : ' Lepton_pdgId[0]==-11 && Lepton_pt[1]>=20 && oneJet && Alt$(CleanJet_pt[1],0)<30',
      'em_mp_1j_pt2ge20' : ' Lepton_pdgId[0]==11 && Lepton_pt[1]>=20 && oneJet && Alt$(CleanJet_pt[1],0)<30',
      'me_mp_1j_pt2ge20' : ' Lepton_pdgId[0]==-13 && Lepton_pt[1]>=20 && oneJet && Alt$(CleanJet_pt[1],0)<30',
      'me_pm_1j_pt2ge20' : ' Lepton_pdgId[0]==13 && Lepton_pt[1]>=20 && oneJet && Alt$(CleanJet_pt[1],0)<30',
      #
      'em_pm_1j_pt2lt20' : ' Lepton_pdgId[0]==-11 && Lepton_pt[1]<20 && oneJet && Alt$(CleanJet_pt[1],0)<30',
      'em_mp_1j_pt2lt20' : ' Lepton_pdgId[0]==11 && Lepton_pt[1]<20 && oneJet && Alt$(CleanJet_pt[1],0)<30',
      'me_mp_1j_pt2lt20' : ' Lepton_pdgId[0]==-13 && Lepton_pt[1]<20 && oneJet && Alt$(CleanJet_pt[1],0)<30',
      'me_pm_1j_pt2lt20' : ' Lepton_pdgId[0]==13 && Lepton_pt[1]<20 && oneJet && Alt$(CleanJet_pt[1],0)<30',
      # FIXME fix the mjj or additional cuts to make this orthogonal to VH2j and VBF
      '2j'               : ' (mjj<65 || mjj>105) && mjj<120 && multiJet', 
   }
}
'''
## Top control regions
cuts['hww2l2v_13TeV_top']  = { 
   'expr' : 'topcr',
    # Define the sub-categorization of topcr
   'categories' : {
      '0j' : 'zeroJet',
      '1j' : 'oneJet && Alt$(CleanJet_pt[1],0)<30',
      '2j' : '(mjj<65 || mjj>105) && mjj<120 && multiJet',
   }
}

## DYtt control regions
cuts['hww2l2v_13TeV_dytt']  = { 
   'expr' : 'dycr',
   # Define the sub-categorization of dycr
   'categories' : { 
      '0j' : 'zeroJet',
      '1j' : 'oneJet && Alt$(CleanJet_pt[1],0)<30',
      '2j' : '(mjj<65 || mjj>105) && mjj<120 && multiJet',
   }
}
'''
## WW control regions
## Used only for control plots, no need to add these cuts for the fit
cuts['hww2l2v_13TeV_WW'] = {
 'expr' : 'wwcr',
 'categories' : {
   '0j' : 'zeroJet',
   '1j' : 'oneJet && Alt$(CleanJet_pt[1],0)<30',
   '2j' : '(mjj<65 || mjj>105) && mjj<120 && multiJet'
 }
}
'''
