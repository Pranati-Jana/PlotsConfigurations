# cuts

supercut = '1'

#cuts['Preselection'] = '1'

#cuts['wh3l_13TeV_sssf']  = 'WH3l_flagOSSF == 0\
#                            && Alt$( CleanJet_pt[0], 0) < 30 \
#                            && Sum$( CleanJet_pt > 20. && abs(CleanJet_eta)<2.5 && Jet_btagDeepB[CleanJet_jetIdx] > 0.6321) == 0\
#                       '

#cuts['wh3l_13TeV_ossf']  = 'WH3l_flagOSSF == 1\
#                            && WH3l_ZVeto > 20\
#                            && Alt$( CleanJet_pt[0], 0) < 30 \
#                            && PuppiMET_pt > 40 \
#                            && Sum$( CleanJet_pt > 20. && abs(CleanJet_eta)<2.5 && Jet_btagDeepB[CleanJet_jetIdx] > 0.6321) == 0\
#                       '
'''
cuts['hww2l2v_13TeV_SRVBF']  = 'bVeto && CleanJet_pt[0]>=30 && CleanJet_pt[1]>=30 && abs(CleanJet_eta[0])<4.7 && abs(CleanJet_eta[1])<4.7 && CleanJet_pt[2]<30 \
                                && nCleanFatJet==0 \
                                && kd_vbf>0.8 \
                                && mjj>200 \
                                && (mth>=30 && mth<125)\
                                && mll>12 \
                                && Lepton_pt[1]>10 \
                                && (nLepton>=2 && Alt$(Lepton_pt[2],0)<10) \
                                && (abs(Lepton_pdgId[1]) == 13 || Lepton_pt[1]>13) \
                                && (Lepton_pdgId[0] * Lepton_pdgId[1] == -11*13) \
                                && ptll > 30 \
                                && PuppiMET_pt > 20 \
                                && hm > 0   \
'
'''
cuts['hww2l2v_13TeV_SRVBF'] = 'mll>12 \
            && Lepton_pt[0]>20 \
            && Lepton_pt[1]>10 \
            && (abs(Lepton_pdgId[0])==13 || Lepton_pt[0]>25) \
            && (abs(Lepton_pdgId[1])==13 || Lepton_pt[1]>13) \
            && (nLepton>=2 && Alt$(Lepton_pt[2],0)<10) \
            && abs(Lepton_eta[0])<2.5 && abs(Lepton_eta[1])<2.5 \
            && ptll>30 \
            && PuppiMET_pt > 20 \
            && hm > 0 \
            && (Lepton_pdgId[0]*Lepton_pdgId[1] == -11*13) \
'


cuts['hww2l2v_13TeV_SRWH_sssf']  = {
 'expr' : 'MinIf$( WH3l_mOSll[], WH3l_mOSll[Iteration$] > 0) > 12 \
            && Alt$(Lepton_pt[0],0)>25 \
            && Alt$(Lepton_pt[1],0)>20 \
            && Alt$(Lepton_pt[2],0)>15 \
            && (nLepton>=3 && Alt$(Lepton_pt[3],0)<10) \
            && abs(WH3l_chlll) == 1 \
            && WH3l_flagOSSF == 0\
            && Alt$( CleanJet_pt[0], 0) < 30 \
            && Sum$( CleanJet_pt > 20. && abs(CleanJet_eta)<2.5 && Jet_btagDeepB[CleanJet_jetIdx] > 0.4941) == 0\
'
}

cuts['hww2l2v_13TeV_SRWH_ossf']  = {
 'expr' : 'MinIf$( WH3l_mOSll[], WH3l_mOSll[Iteration$] > 0) > 12 \
            && Alt$(Lepton_pt[0],0)>25 \
            && Alt$(Lepton_pt[1],0)>20 \
            && Alt$(Lepton_pt[2],0)>15 \
            && (nLepton>=3 && Alt$(Lepton_pt[3],0)<10) \
            && abs(WH3l_chlll) == 1 \
            && WH3l_flagOSSF == 1\
            && WH3l_ZVeto > 20\
            && Alt$( CleanJet_pt[0], 0) < 30 \
            && PuppiMET_pt > 40 \
            && Sum$( CleanJet_pt > 20. && abs(CleanJet_eta)<2.5 && Jet_btagDeepB[CleanJet_jetIdx] > 0.4941) == 0\
'
}

cuts['hww2l2v_13TeV_SRZH']  = 'Alt$(Lepton_pt[0],0)>20 \
          && Alt$(Lepton_pt[1],0)>15 \
          && Alt$(Lepton_pt[2],0)>10 \
          && Alt$(Lepton_pt[3],0)>10 \
          && (nLepton>=4 && Alt$(Lepton_pt[4],0)<10) \
          && chllll_zh4l == 0 \
          && z0Mass_zh4l>12 \
          && abs(z0Mass_zh4l-91.1876)< 15\
          && z1Mass_zh4l < 70 && z1Mass_zh4l >10 \
          && PuppiMET_pt > 35 \
          && (Sum$(CleanJet_pt > 20. && abs(CleanJet_eta)<2.5 && Jet_btagDeepB[CleanJet_jetIdx] > 0.2217) == 0)\
'
cuts['hww2l2v_13TeV_merged'] = '( ( %s ) || ( %s ) )' %( cuts['hww2l2v_13TeV_SRWH_sssf']['expr'] , cuts['hww2l2v_13TeV_SRWH_ossf']['expr'])

cuts['hww2l2v_13TeV_SRWH']  = 'MinIf$(WH3l_mOSll[], WH3l_mOSll[Iteration$] > 0) > 12\
                                && Alt$(Lepton_pt[0],0)>25 \
                                && Alt$(Lepton_pt[1],0)>20 \
                                && Alt$(Lepton_pt[2],0)>15 \
                                && (nLepton>=3 && Alt$(Lepton_pt[3],0)<10) \
                                && abs(WH3l_chlll) == 1 \
                                && Alt$( CleanJet_pt[0], 0) < 30 \
                                && Sum$( CleanJet_pt > 20. && abs(CleanJet_eta)<2.5 && Jet_btagDeepB[CleanJet_jetIdx] > 0.4941) == 0\
'

'''
cuts['hww2l2v_13TeV_SRVH']   = 'bVeto && CleanJet_pt[0]>=30 && CleanJet_pt[1]>=30 && abs(CleanJet_eta[0])<2.4 && abs(CleanJet_eta[1])<2.4 && CleanJet_pt[2]<30 \
                                && nCleanFatJet==0 \
                                && kd_vh>0.8 \
                                && (mjj>60 && mjj<120) \
                                && (mth>=30 && mth<125)\
                                && mll>12 \
                                && Lepton_pt[1]>10 \
                                && (nLepton>=2 && Alt$(Lepton_pt[2],0)<10) \
                                && (abs(Lepton_pdgId[1]) == 13 || Lepton_pt[1]>13) \
                                && (Lepton_pdgId[0] * Lepton_pdgId[1] == -11*13) \
                                && ptll > 30 \
                                && PuppiMET_pt > 20 \
                                && hm > 0  \
'
'''

cuts['hww2l2v_13TeV_wz'] = 'MinIf$( WH3l_mOSll[], WH3l_mOSll[Iteration$] > 0) > 12 \
                         && Alt$(Lepton_pt[0],0)>25 \
                         && Alt$(Lepton_pt[1],0)>20 \
                         && Alt$(Lepton_pt[2],0)>15 \
                         && (nLepton>=3 && Alt$(Lepton_pt[3],0)<10) \
                         && abs(WH3l_chlll) == 1 \
                         && WH3l_flagOSSF == 1\
                         && PuppiMET_pt > 45\
                         && WH3l_ZVeto < 20\
                         && WH3l_mlll > 100\
                         && Alt$( CleanJet_pt[0], 0) < 30 \
                         && Sum$( CleanJet_pt > 20. && abs(CleanJet_eta)<2.5 && Jet_btagDeepB[CleanJet_jetIdx] > 0.4941) == 0\
'                        

cuts['hww2l2v_13TeV_zg'] = 'MinIf$( WH3l_mOSll[], WH3l_mOSll[Iteration$] > 0) > 12 \
                         && Alt$(Lepton_pt[0],0)>25 \
                         && Alt$(Lepton_pt[1],0)>20 \
                         && Alt$(Lepton_pt[2],0)>15 \
                         && (nLepton>=3 && Alt$(Lepton_pt[3],0)<10) \
                         && abs(WH3l_chlll) == 1 \
                         && WH3l_ZVeto < 20\
                         && PuppiMET_pt < 40\
                         && WH3l_mlll > 80\
                         && WH3l_mlll < 100\
                         && Alt$( CleanJet_pt[0], 0) < 30 \
                         && Sum$( CleanJet_pt > 20. && abs(CleanJet_eta)<2.5 && Jet_btagDeepB[CleanJet_jetIdx] > 0.4941) == 0\
'

# 11 = e
# 13 = mu
# 15 = tau

