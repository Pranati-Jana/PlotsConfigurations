
supercut = '    mll>12 \
            && Lepton_pt[0]>25 \
            && Lepton_pt[1]>10 \
            && (abs(Lepton_pdgId[1])==13 || Lepton_pt[1]>13) \
            && (nLepton>=2 && Alt$(Lepton_pt[2],0)<10) \
            && abs(Lepton_eta[0])<2.5 && abs(Lepton_eta[1])<2.5 \
            && ptll>30 \
            && PuppiMET_pt > 20 \
           '

# Some cuts

dymva0jet = 'dymva_dnn_0j  > 0.970'
dymva1jet = 'dymva_dnn_1j  > 0.970'
dymva2jet = 'dymva_dnn_2j  > 0.975'
dymvaVBF  = 'dymva_dnn_VBF > 0.975'
dymvaVH   = 'dymva_dnn_VH  > 0.830'

# Higgs Signal Regions: ee/uu * 0/1/2 jet
cuts['hww2l2v_13TeV'] = {
   'expr': 'sr && (Lepton_pdgId[0]==-Lepton_pdgId[1])' ,
   'categories' : {
       '0j_ee_pt2lt20'         : 'zeroJet && (Lepton_pdgId[0]*Lepton_pdgId[1] == -11*11) && Higgs0jet && Lepton_pt[1]< 20 && '+dymva0jet,
       '0j_mm_pt2lt20'         : 'zeroJet && (Lepton_pdgId[0]*Lepton_pdgId[1] == -13*13) && Higgs0jet && Lepton_pt[1]< 20 && '+dymva0jet,
       '0j_ee_pt2ge20'         : 'zeroJet && (Lepton_pdgId[0]*Lepton_pdgId[1] == -11*11) && Higgs0jet && Lepton_pt[1]>=20 && '+dymva0jet,
       '0j_mm_pt2ge20'         : 'zeroJet && (Lepton_pdgId[0]*Lepton_pdgId[1] == -13*13) && Higgs0jet && Lepton_pt[1]>=20 && '+dymva0jet,
       '1j_ee'                 : ' oneJet && (Lepton_pdgId[0]*Lepton_pdgId[1] == -11*11) && Higgs1jet && '+dymva1jet,
       '1j_mm'                 : ' oneJet && (Lepton_pdgId[0]*Lepton_pdgId[1] == -13*13) && Higgs1jet && '+dymva1jet,
       '2j_ee'                 : '  2jggH && (Lepton_pdgId[0]*Lepton_pdgId[1] == -11*11) && Higgs2jet && '+dymva2jet,
       '2j_mm'                 : '  2jggH && (Lepton_pdgId[0]*Lepton_pdgId[1] == -13*13) && Higgs2jet && '+dymva2jet,
   }
}

# Numerator for DY acceptance in Signal region
cuts['hww2l2v_13TeV_HAccNum'] = {
   'expr': 'sr && (Lepton_pdgId[0]==-Lepton_pdgId[1])' ,
   'categories' : {
       '0j_ee_pt2lt20' : 'zeroJet && (Lepton_pdgId[0]*Lepton_pdgId[1] == -11*11) && Higgs0jet && Lepton_pt[1]< 20 && dymva_dnn_0j > 0.8',
       '0j_mm_pt2lt20' : 'zeroJet && (Lepton_pdgId[0]*Lepton_pdgId[1] == -13*13) && Higgs0jet && Lepton_pt[1]< 20 && dymva_dnn_0j > 0.8',
       '0j_ee_pt2ge20' : 'zeroJet && (Lepton_pdgId[0]*Lepton_pdgId[1] == -11*11) && Higgs0jet && Lepton_pt[1]>=20 && dymva_dnn_0j > 0.8',
       '0j_mm_pt2ge20' : 'zeroJet && (Lepton_pdgId[0]*Lepton_pdgId[1] == -13*13) && Higgs0jet && Lepton_pt[1]>=20 && dymva_dnn_0j > 0.8',
       '1j_ee'         : ' oneJet && (Lepton_pdgId[0]*Lepton_pdgId[1] == -11*11) && Higgs1jet && dymva_dnn_1j > 0.8',
       '1j_mm'         : ' oneJet && (Lepton_pdgId[0]*Lepton_pdgId[1] == -13*13) && Higgs1jet && dymva_dnn_1j > 0.8',
       '2j_ee'         : '  2jggH && (Lepton_pdgId[0]*Lepton_pdgId[1] == -11*11) && Higgs2jet && dymva_dnn_2j > 0.8',
       '2j_mm'         : '  2jggH && (Lepton_pdgId[0]*Lepton_pdgId[1] == -13*13) && Higgs2jet && dymva_dnn_2j > 0.8',
   }
}


## Top CR: No H sel , bTag , tight DYmva
cuts['hww2l2v_13TeV_top']  = {
   'expr' : 'topcr && ZVeto * (Lepton_pdgId[0]==-Lepton_pdgId[1])',
   'categories' : {
      '0j_ee' : 'zeroJet && (Lepton_pdgId[0]*Lepton_pdgId[1] == -11*11) && '+dymva0jet,
      '0j_mm' : 'zeroJet && (Lepton_pdgId[0]*Lepton_pdgId[1] == -13*13) && '+dymva0jet,
      '1j_ee' : ' oneJet && (Lepton_pdgId[0]*Lepton_pdgId[1] == -11*11) && '+dymva1jet,
      '1j_mm' : ' oneJet && (Lepton_pdgId[0]*Lepton_pdgId[1] == -13*13) && '+dymva1jet,
      '2j_ee' : '  2jggH && (Lepton_pdgId[0]*Lepton_pdgId[1] == -11*11) && '+dymva2jet,
      '2j_mm' : '  2jggH && (Lepton_pdgId[0]*Lepton_pdgId[1] == -13*13) && '+dymva2jet,
   }
}

## WW CR: No H Sel , mll>80, tight DYMva
cuts['hww2l2v_13TeV_WW']  = {
   'expr' : 'wwcr * (Lepton_pdgId[0]==-Lepton_pdgId[1])',
   'categories' : {
      '0j_ee' : 'zeroJet && (Lepton_pdgId[0]*Lepton_pdgId[1] == -11*11) && '+dymva0jet,
      '0j_mm' : 'zeroJet && (Lepton_pdgId[0]*Lepton_pdgId[1] == -13*13) && '+dymva0jet,
      '1j_ee' : ' oneJet && (Lepton_pdgId[0]*Lepton_pdgId[1] == -11*11) && '+dymva1jet,
      '1j_mm' : ' oneJet && (Lepton_pdgId[0]*Lepton_pdgId[1] == -13*13) && '+dymva1jet,
      '2j_ee' : '  2jggH && (Lepton_pdgId[0]*Lepton_pdgId[1] == -11*11) && '+dymva2jet,
      '2j_mm' : '  2jggH && (Lepton_pdgId[0]*Lepton_pdgId[1] == -13*13) && '+dymva2jet,
   }
}

## DY Background IN with DYMVA>0.9X : Split ee/mm , No H cut !
cuts['hww2l2v_13TeV_DYin']  = {
   'expr' : 'Zpeak && bVeto',
   'categories' : {
      '0j_ee' : 'zeroJet && (Lepton_pdgId[0]*Lepton_pdgId[1] == -11*11) && '+dymva0jet,
      '0j_mm' : 'zeroJet && (Lepton_pdgId[0]*Lepton_pdgId[1] == -13*13) && '+dymva0jet,
      '0j_df' : 'zeroJet && (Lepton_pdgId[0]*Lepton_pdgId[1] == -11*13) && '+dymva0jet,
      '1j_ee' : ' oneJet && (Lepton_pdgId[0]*Lepton_pdgId[1] == -11*11) && '+dymva1jet,
      '1j_mm' : ' oneJet && (Lepton_pdgId[0]*Lepton_pdgId[1] == -13*13) && '+dymva1jet,
      '1j_df' : ' oneJet && (Lepton_pdgId[0]*Lepton_pdgId[1] == -11*13) && '+dymva1jet,
      '2j_ee' : '  2jggH && (Lepton_pdgId[0]*Lepton_pdgId[1] == -11*11) && '+dymva2jet,
      '2j_mm' : '  2jggH && (Lepton_pdgId[0]*Lepton_pdgId[1] == -13*13) && '+dymva2jet,
      '2j_df' : '  2jggH && (Lepton_pdgId[0]*Lepton_pdgId[1] == -11*13) && '+dymva2jet,
   }
}

## Acc Denominator
cuts['hww2l2v_13TeV_AccDen']  = {
   'expr' : 'sr * (Lepton_pdgId[0]==-Lepton_pdgId[1])',
   'categories' : {
      '0j_ee' : 'zeroJet && (Lepton_pdgId[0]*Lepton_pdgId[1] == -11*11) && dymva_dnn_0j > 0.8',
      '0j_mm' : 'zeroJet && (Lepton_pdgId[0]*Lepton_pdgId[1] == -13*13) && dymva_dnn_0j > 0.8',
      '1j_ee' : ' oneJet && (Lepton_pdgId[0]*Lepton_pdgId[1] == -11*11) && dymva_dnn_1j > 0.8',
      '1j_mm' : ' oneJet && (Lepton_pdgId[0]*Lepton_pdgId[1] == -13*13) && dymva_dnn_1j > 0.8',
      '2j_ee' : '  2jggH && (Lepton_pdgId[0]*Lepton_pdgId[1] == -11*11) && dymva_dnn_2j > 0.8',
      '2j_mm' : '  2jggH && (Lepton_pdgId[0]*Lepton_pdgId[1] == -13*13) && dymva_dnn_2j > 0.8',
   }
}

## WW Acc Numerator: loose dymva + WW sel
cuts['hww2l2v_13TeV_wwAcc']  = {
   'expr' : 'wwcr * (Lepton_pdgId[0]==-Lepton_pdgId[1])',
   'categories' : {
      '0j_ee' : 'zeroJet && (Lepton_pdgId[0]*Lepton_pdgId[1] == -11*11) && dymva_dnn_0j > 0.8',
      '0j_mm' : 'zeroJet && (Lepton_pdgId[0]*Lepton_pdgId[1] == -13*13) && dymva_dnn_0j > 0.8',
      '1j_ee' : ' oneJet && (Lepton_pdgId[0]*Lepton_pdgId[1] == -11*11) && dymva_dnn_1j > 0.8',
      '1j_mm' : ' oneJet && (Lepton_pdgId[0]*Lepton_pdgId[1] == -13*13) && dymva_dnn_1j > 0.8',
      '2j_ee' : '  2jggH && (Lepton_pdgId[0]*Lepton_pdgId[1] == -11*11) && dymva_dnn_2j > 0.8',
      '2j_mm' : '  2jggH && (Lepton_pdgId[0]*Lepton_pdgId[1] == -13*13) && dymva_dnn_2j > 0.8',
   }
} 