
import os
import copy
import inspect
import ROOT

configurations = os.path.realpath(inspect.getfile(inspect.currentframe())) # this file
configurations = os.path.dirname(configurations) # Full2016
configurations = os.path.dirname(configurations) # VBF
configurations = os.path.dirname(configurations) # EFT
configurations = os.path.dirname(configurations) # Configurations

mc = [skey for skey in samples if skey not in ('Fake', 'DATA')]
mc_ggh = [skey for skey in samples if skey.startswith('H0')]
mc_ggh.append('ggH_hww')
mc_qqh = [skey for skey in samples if skey.startswith('VBF_H0')]


eleWP = 'mvaFall17V1Iso_WP90'
muWP = 'cut_Tight_HWWW'

eleFWP = 'mvaFall17V1Iso_WP90_tthmva_70'
muFWP = 'cut_Tight_HWWW_tthmva_80'
btag_algo="deepflav"


aliases['LepWPCut'] = {
    'expr': 'LepCut3l__ele_'+eleWP+'__mu_'+muWP+'*(((abs(Lepton_pdgId[0])==13 && Muon_mvaTTH[Lepton_muonIdx[0]]>0.8) || (abs(Lepton_pdgId[0])==11 && Electron_mvaTTH[Lepton_electronIdx[0]]>0.7)) && ((abs(Lepton_pdgId[1])==13 && Muon_mvaTTH[Lepton_muonIdx[1]]>0.8) || (abs(Lepton_pdgId[1])==11 && Electron_mvaTTH[Lepton_electronIdx[1]]>0.7)) && ((abs(Lepton_pdgId[2])==13 && Muon_mvaTTH[Lepton_muonIdx[2]]>0.8) || (abs(Lepton_pdgId[2])==11 && Electron_mvaTTH[Lepton_electronIdx[2]]>0.7)))',
    # 'expr': 'LepCut2l__ele_'+eleWP+'__mu_'+muWP+'*( (abs(Lepton_pdgId[0])==11 || Muon_mvaTTH[Lepton_muonIdx[0]]>0.8) && (abs(Lepton_pdgId[1])==11 || Muon_mvaTTH[Lepton_muonIdx[1]]>0.8) )',
    'samples': mc + ['DATA']
}

#aliases['LepWPCut'] = {
#    'expr': 'LepCut2l__ele_'+eleWP+'__mu_'+muWP,
#    'samples': mc + ['DATA']
#}

aliases['gstarLow'] = {
    'expr': 'Gen_ZGstar_mass >0 && Gen_ZGstar_mass < 4',
    'samples': 'VgS'
}

aliases['gstarHigh'] = {
    'expr': 'Gen_ZGstar_mass <0 || Gen_ZGstar_mass > 4',
    'samples': 'VgS'
}

aliases['embedtotal'] = {
    'expr': 'embed_total_WP90V1',  # wrt. eleWP
    'samples': 'Dyemb'
}

# Fake leptons transfer factor
aliases['fakeW'] = {
    'expr': 'fakeW_ele_'+eleFWP+'_mu_'+muFWP+'_3l',
    'samples': ['Fake']
}

# And variations - already divided by central values in formulas !
aliases['fakeWEleUp'] = {
    'expr': 'fakeW_ele_'+eleFWP+'_mu_'+muFWP+'_3lElUp/fakeW_ele_'+eleFWP+'_mu_'+muFWP+'_3l',
    'samples': ['Fake']
}
aliases['fakeWEleDown'] = {
    'expr': 'fakeW_ele_'+eleFWP+'_mu_'+muFWP+'_3lElDown/fakeW_ele_'+eleFWP+'_mu_'+muFWP+'_3l',
    'samples': ['Fake']
}
aliases['fakeWMuUp'] = {
    'expr': 'fakeW_ele_'+eleFWP+'_mu_'+muFWP+'_3lMuUp/fakeW_ele_'+eleFWP+'_mu_'+muFWP+'_3l',
    'samples': ['Fake']
}
aliases['fakeWMuDown'] = {
    'expr': 'fakeW_ele_'+eleFWP+'_mu_'+muFWP+'_3lMuDown/fakeW_ele_'+eleFWP+'_mu_'+muFWP+'_3l',
    'samples': ['Fake']
}
aliases['fakeWStatEleUp'] = {
    'expr': 'fakeW_ele_'+eleFWP+'_mu_'+muFWP+'_3lstatElUp/fakeW_ele_'+eleFWP+'_mu_'+muFWP+'_3l',
    'samples': ['Fake']
}
aliases['fakeWStatEleDown'] = {
    'expr': 'fakeW_ele_'+eleFWP+'_mu_'+muFWP+'_3lstatElDown/fakeW_ele_'+eleFWP+'_mu_'+muFWP+'_3l',
    'samples': ['Fake']
}
aliases['fakeWStatMuUp'] = {
    'expr': 'fakeW_ele_'+eleFWP+'_mu_'+muFWP+'_3lstatMuUp/fakeW_ele_'+eleFWP+'_mu_'+muFWP+'_3l',
    'samples': ['Fake']
}
aliases['fakeWStatMuDown'] = {
    'expr': 'fakeW_ele_'+eleFWP+'_mu_'+muFWP+'_3lstatMuDown/fakeW_ele_'+eleFWP+'_mu_'+muFWP+'_3l',
    'samples': ['Fake']
}

# gen-matching to prompt only (GenLepMatch2l matches to *any* gen lepton)
aliases['PromptGenLepMatch3l'] = {
    'expr': 'Alt$(Lepton_promptgenmatched[0]*Lepton_promptgenmatched[1]*Lepton_promptgenmatched[2], 0)',
    'samples': mc
}

aliases['Top_pTrw'] = {
    'expr': '(topGenPt * antitopGenPt > 0.) * (TMath::Sqrt(TMath::Exp(0.0615 - 0.0005 * topGenPt) * TMath::Exp(0.0615 - 0.0005 * antitopGenPt))) + (topGenPt * antitopGenPt <= 0.)',
    'samples': ['top']
}

aliases['nCleanGenJet'] = {
    'linesToAdd': ['.L %s/src/PlotsConfigurations/Configurations/Differential/ngenjet.cc+' % os.getenv('CMSSW_BASE')
    ],
    'class': 'CountGenJet',
    'samples': mc
}


aliases['ptllDYW_NLO'] = {
    'expr': '(0.87*(gen_ptll<10)+(0.379119+0.099744*gen_ptll-0.00487351*gen_ptll**2+9.19509e-05*gen_ptll**3-6.0212e-07*gen_ptll**4)*(gen_ptll>=10 && gen_ptll<45)+(9.12137e-01+1.11957e-04*gen_ptll-3.15325e-06*gen_ptll**2-4.29708e-09*gen_ptll**3+3.35791e-11*gen_ptll**4)*(gen_ptll>=45 && gen_ptll<200) + 1*(gen_ptll>200))',
    'samples': ['DY']
}

##### DY Z pT reweighting
aliases['getGenZpt_OTF'] = {
    'linesToAdd':['.L %s/src/PlotsConfigurations/Configurations/patches/getGenZpt.cc+' % os.getenv('CMSSW_BASE')],
    'class': 'getGenZpt',
    'samples': ['DY']
}
handle = open('%s/src/PlotsConfigurations/Configurations/patches/DYrew30.py' % os.getenv('CMSSW_BASE'),'r')
exec(handle)
handle.close()
aliases['DY_NLO_pTllrw'] = {
    'expr': '('+DYrew['2017']['NLO'].replace('x', 'getGenZpt_OTF')+')*(nCleanGenJet == 0)+1.0*(nCleanGenJet > 0)',
    'samples': ['DY']
}
aliases['DY_LO_pTllrw'] = {
    'expr': '('+DYrew['2017']['LO'].replace('x', 'getGenZpt_OTF')+')*(nCleanGenJet == 0)+1.0*(nCleanGenJet > 0)',
    'samples': ['DY']
}


# Jet bins
# using Alt$(CleanJet_pt[n], 0) instead of Sum$(CleanJet_pt >= 30) because jet pt ordering is not strictly followed in JES-varied samples

# No jet with pt > 30 GeV
aliases['zeroJet'] = {
    'expr': 'Alt$(CleanJet_pt[0], 0) < 30.'
}

aliases['oneJet'] = {
    'expr': 'Alt$(CleanJet_pt[0], 0) > 30.'
}

aliases['multiJet'] = {
    'expr': 'Alt$(CleanJet_pt[1], 0) > 30.'
}

#################### B tag ##############################3
aliases['bVeto'] = {
    'expr': 'Sum$(CleanJet_pt > 20. && abs(CleanJet_eta) < 2.5 && Jet_btagDeepB[CleanJet_jetIdx] > 0.4941) == 0'
}

aliases['bReq'] = {
    'expr': 'Sum$(CleanJet_pt > 30. && abs(CleanJet_eta) < 2.5 && Jet_btagDeepB[CleanJet_jetIdx] > 0.4941) >= 1'
}

aliases['bVetoSF'] = {
    'expr': 'TMath::Exp(Sum$(TMath::Log((CleanJet_pt>20 && abs(CleanJet_eta)<2.5)*Jet_btagSF_deepcsv_shape[CleanJet_jetIdx]+1*(CleanJet_pt<20 || abs(CleanJet_eta)>2.5))))',
    'samples': mc
}

aliases['bReqSF'] = {
    'expr': 'TMath::Exp(Sum$(TMath::Log((CleanJet_pt>30 && abs(CleanJet_eta)<2.5)*Jet_btagSF_deepcsv_shape[CleanJet_jetIdx]+1*(CleanJet_pt<30 || abs(CleanJet_eta)>2.5))))',
    'samples': mc
}

aliases['btagSF'] = {
    'expr': 'bVetoSF*bVeto',
    'samples': mc
}


for shift in ['jes', 'lf', 'hf', 'lfstats1', 'lfstats2', 'hfstats1', 'hfstats2', 'cferr1', 'cferr2']:
    for targ in ['bVeto', 'bReq']:
        alias = aliases['%sSF%sup' % (targ, shift)] = copy.deepcopy(aliases['%sSF' % targ])
        alias['expr'] = alias['expr'].replace('btagSF_shape', 'btagSF_shape_up_%s' % shift)

        alias = aliases['%sSF%sdown' % (targ, shift)] = copy.deepcopy(aliases['%sSF' % targ])
        alias['expr'] = alias['expr'].replace('btagSF_shape', 'btagSF_shape_down_%s' % shift)

    aliases['btagSF%sup' % shift] = {
        'expr': aliases['btagSF']['expr'].replace('SF', 'SF' + shift + 'up'),
        'samples': mc
    }

    aliases['btagSF%sdown' % shift] = {
        'expr': aliases['btagSF']['expr'].replace('SF', 'SF' + shift + 'down'),
        'samples': mc
    }



# Lepton centrality

aliases['lepcen1'] = {
    'expr': 'abs((Lepton_eta[0] - (CleanJet_eta[0]+CleanJet_eta[1])/2)/detajj)'
}

aliases['lepcen2'] = {
    'expr': 'abs((Lepton_eta[1] - (CleanJet_eta[0]+CleanJet_eta[1])/2)/detajj)'
}

################################################# EFT with MELA #############################################
# Couplings (gXHWW^2 = JHUXSHWWa1/JHUXSHWWaX) and cross-sections (JHUXSHWWaX) 
# https://github.com/hroskes/anomalouscouplingsconstants/blob/04fc990ad2452c79506de79474fe2c83243bb39f/constants.py#L75-L93
# https://twiki.cern.ch/twiki/bin/view/CMS/Run2MCProductionforHiggsProperties#POWHEG_ggH_Production_JHUGen_H_W

# Couplings used in JHUGen 

aliases['G2_HWW']  = {'expr': '1.133582'}
aliases['G4_HWW']  = {'expr': '1.76132'}
aliases['L1_HWW']  = {'expr': '-13752.22'}

aliases['G2_VBF']  = {'expr': '0.27196538'}
aliases['G4_VBF']  = {'expr': '0.29797901870'}
aliases['L1_VBF']  = {'expr': '-2158.21307286'}

aliases['G2_WH']   = {'expr': '0.0998956'}
aliases['G4_WH']   = {'expr': '0.1236136'}
aliases['L1_WH']   = {'expr': '-525.274'}

aliases['G2_ZH']   = {'expr': '0.112481'}
aliases['G4_ZH']   = {'expr': '0.144057'}
aliases['L1_ZH']   = {'expr': '-517.788'}

# Cross-sections : Decay 

aliases['JHUXSHWWa1']   = {'expr': '312.04019'}
aliases['JHUXSHWWa2']   = {'expr': '242.6283'}
aliases['JHUXSHWWa3']   = {'expr': '100.79135'} 
aliases['JHUXSHWWL1']   = {'expr': '1.6475531e-06'}
aliases['JHUXSHWWa1a2'] = {'expr': '1149.9181'}
aliases['JHUXSHWWa1a3'] = {'expr': '624.7195'}
aliases['JHUXSHWWa1L1'] = {'expr': '5.3585509'}

# Cross-sections : Production

aliases['JHUXSVBFa1']   = {'expr': '968.88143'}
aliases['JHUXSVBFa2']   = {'expr': '13097.831'}
aliases['JHUXSVBFa3']   = {'expr': '10910.237'}
aliases['JHUXSVBFL1']   = {'expr': '0.00020829261'}
aliases['JHUXSVBFa1a2'] = {'expr': '2207.6738'}
aliases['JHUXSVBFa1a3'] = {'expr': '1936.4327'}
aliases['JHUXSVBFa1L1'] = {'expr': '2861.7003'}

aliases['JHUXSWHa1']   = {'expr': '14813072'}
aliases['JHUXSWHa2']   = {'expr': '1.4845783e+09'}
aliases['JHUXSWHa3']   = {'expr': '9.6943028e+08'}
aliases['JHUXSWHL1']   = {'expr': '53.687994'}
aliases['JHUXSWHa1a2'] = {'expr': '7879980.3'}
aliases['JHUXSWHa1a3'] = {'expr': '29626131'}
aliases['JHUXSWHa1L1'] = {'expr': '12092167'}

aliases['JHUXSZHa1']   = {'expr': '1436880.4'}
aliases['JHUXSZHa2']   = {'expr': '1.1360424e+08'}
aliases['JHUXSZHa3']   = {'expr': '69241514'}
aliases['JHUXSZHL1']   = {'expr': '5.3610896'}
aliases['JHUXSZHa1a2'] = {'expr': '678434.94'}
aliases['JHUXSZHa1a3'] = {'expr': '2873685.9'}
aliases['JHUXSZHa1L1'] = {'expr': '1091656.8'}

# Normalisation Weights

aliases['H0PM_W']    = { 'expr': '1'}
aliases['H0M_W']     = { 'expr': '(JHUXSHWWa3/JHUXSHWWa1)'}
aliases['H0PH_W']    = { 'expr': '(JHUXSHWWa2/JHUXSHWWa1)'}
aliases['H0L1_W']    = { 'expr': '(JHUXSHWWL1/JHUXSHWWa1)'}
aliases['H0Mf05_W']  = { 'expr': '(JHUXSHWWa1a3/JHUXSHWWa1)'}
aliases['H0PHf05_W'] = { 'expr': '(JHUXSHWWa1a2/JHUXSHWWa1)'}
aliases['H0L1f05_W'] = { 'expr': '(JHUXSHWWa1L1/JHUXSHWWa1)'}

aliases['JHUXSHWWa1a2_I'] = {'expr':'(JHUXSHWWa1a2 - JHUXSHWWa1 - (G2_HWW**2)*JHUXSHWWa2)/G2_HWW'}
aliases['JHUXSHWWa1a3_I'] = {'expr':'(JHUXSHWWa1a3 - JHUXSHWWa1 - (G4_HWW**2)*JHUXSHWWa3)/G4_HWW'}
aliases['JHUXSHWWa1L1_I'] = {'expr':'(JHUXSHWWa1L1 - JHUXSHWWa1 - (L1_HWW**2)*JHUXSHWWL1)/L1_HWW'}

aliases['H0Mf05VBF_W']  = { 'expr': '(JHUXSHWWa1 + JHUXSHWWa1a3_I*G4_VBF + JHUXSHWWa3*(G4_VBF**2))/JHUXSHWWa1'}
aliases['H0PHf05VBF_W'] = { 'expr': '(JHUXSHWWa1 + JHUXSHWWa1a2_I*G2_VBF + JHUXSHWWa2*(G2_VBF**2))/JHUXSHWWa1'}
aliases['H0L1f05VBF_W'] = { 'expr': '(JHUXSHWWa1 + JHUXSHWWa1L1_I*L1_VBF + JHUXSHWWL1*(L1_VBF**2))/JHUXSHWWa1'}

aliases['H0Mf05WH_W']  = { 'expr': '(JHUXSHWWa1 + JHUXSHWWa1a3_I*G4_WH + JHUXSHWWa3*(G4_WH**2))/JHUXSHWWa1'}
aliases['H0PHf05WH_W'] = { 'expr': '(JHUXSHWWa1 + JHUXSHWWa1a2_I*G2_WH + JHUXSHWWa2*(G2_WH**2))/JHUXSHWWa1'}
aliases['H0L1f05WH_W'] = { 'expr': '(JHUXSHWWa1 + JHUXSHWWa1L1_I*L1_WH + JHUXSHWWL1*(L1_WH**2))/JHUXSHWWa1'}

aliases['H0Mf05ZH_W']  = { 'expr': '(JHUXSHWWa1 + JHUXSHWWa1a3_I*G4_ZH + JHUXSHWWa3*(G4_ZH**2))/JHUXSHWWa1'}
aliases['H0PHf05ZH_W'] = { 'expr': '(JHUXSHWWa1 + JHUXSHWWa1a2_I*G2_ZH + JHUXSHWWa2*(G2_ZH**2))/JHUXSHWWa1'}
aliases['H0L1f05ZH_W'] = { 'expr': '(JHUXSHWWa1 + JHUXSHWWa1L1_I*L1_ZH + JHUXSHWWL1*(L1_ZH**2))/JHUXSHWWa1'}

aliases['VBF_H0PM_W']    = { 'expr': '1'}
aliases['VBF_H0M_W']     = { 'expr': 'H0M_W*(JHUXSVBFa3/JHUXSVBFa1)'}
aliases['VBF_H0PH_W']    = { 'expr': 'H0PH_W*(JHUXSVBFa2/JHUXSVBFa1)'}
aliases['VBF_H0L1_W']    = { 'expr': 'H0L1_W*(JHUXSVBFL1/JHUXSVBFa1)'}
aliases['VBF_H0Mf05_W']  = { 'expr': 'H0Mf05VBF_W*(JHUXSVBFa1a3/JHUXSVBFa1)'}
aliases['VBF_H0PHf05_W'] = { 'expr': 'H0PHf05VBF_W*(JHUXSVBFa1a2/JHUXSVBFa1)'}
aliases['VBF_H0L1f05_W'] = { 'expr': 'H0L1f05VBF_W*(JHUXSVBFa1L1/JHUXSVBFa1)'}

aliases['WH_H0PM_W']    = { 'expr': '1'}
aliases['WH_H0M_W']     = { 'expr': 'H0M_W*(JHUXSWHa3/JHUXSWHa1)'}
aliases['WH_H0PH_W']    = { 'expr': 'H0PH_W*(JHUXSWHa2/JHUXSWHa1)'}
aliases['WH_H0L1_W']    = { 'expr': 'H0L1_W*(JHUXSWHL1/JHUXSWHa1)'}
aliases['WH_H0Mf05_W']  = { 'expr': 'H0Mf05WH_W*(JHUXSWHa1a3/JHUXSWHa1)'}
aliases['WH_H0PHf05_W'] = { 'expr': 'H0PHf05WH_W*(JHUXSWHa1a2/JHUXSWHa1)'}
aliases['WH_H0L1f05_W'] = { 'expr': 'H0L1f05WH_W*(JHUXSWHa1L1/JHUXSWHa1)'}

aliases['ZH_H0PM_W']    = { 'expr': '1'}
aliases['ZH_H0M_W']     = { 'expr': 'H0M_W*(JHUXSZHa3/JHUXSZHa1)'}
aliases['ZH_H0PH_W']    = { 'expr': 'H0PH_W*(JHUXSZHa2/JHUXSZHa1)'}
aliases['ZH_H0L1_W']    = { 'expr': 'H0L1_W*(JHUXSZHL1/JHUXSZHa1)'}
aliases['ZH_H0Mf05_W']  = { 'expr': 'H0Mf05ZH_W*(JHUXSZHa1a3/JHUXSZHa1)'}
aliases['ZH_H0PHf05_W'] = { 'expr': 'H0PHf05ZH_W*(JHUXSZHa1a2/JHUXSZHa1)'}
aliases['ZH_H0L1f05_W'] = { 'expr': 'H0L1f05ZH_W*(JHUXSZHa1L1/JHUXSZHa1)'}

# Get MEs for signal reweighting

mes = [ 'ME_H0PM',
    'ME_H0M', 'ME_H0M_M0', 'ME_H0M_M1', 'ME_H0M_M2', 'ME_H0M_M3', 'ME_H0Mf05', 
    'ME_H0PH','ME_H0PH_M0','ME_H0PH_M1','ME_H0PH_M2','ME_H0PH_M3','ME_H0PHf05',
    'ME_H0L1','ME_H0L1_M0','ME_H0L1_M1','ME_H0L1_M2','ME_H0L1_M3','ME_H0L1f05',
]

for me in mes:
    aliases[me] = {
    'linesToAdd': ['.L %s/EFT/VBF/Tools/getme.cc+' % configurations ],
    'class': 'GetME', 'samples': signals_rw, 'args': (me,),
}

# Constants as a function of hm 

cons = [
    'CVBF','CWH','CZH',
    'G4VBF','G4WH','G4ZH','G4VH',
    'G2VBF','G2WH','G2ZH','G2VH',
    'L1VBF','L1WH','L1ZH',
]

for con in cons:
    aliases[con] = {
    'linesToAdd': ['.L %s/EFT/VBF/Tools/getconstant.cc+' % configurations ],
    'class': 'GetConstant',
    'args': (con,),
}

# VBF KDs

aliases['kd_smvbf']     = { 'expr': '1/(1+((me_qcd_hsm*CVBF)/me_vbf_hsm))' }
aliases['kd_hmvbf']     = { 'expr': '1/(1+((me_qcd_hsm*CVBF)/(me_vbf_hm*G4VBF**2)))' }
aliases['kd_hpvbf']     = { 'expr': '1/(1+((me_qcd_hsm*CVBF)/(me_vbf_hp*G2VBF**2)))' }
aliases['kd_hlvbf']     = { 'expr': '1/(1+((me_qcd_hsm*CVBF)/(me_vbf_hl*L1VBF**2)))' }
aliases['kd_vbf']       = { 'expr': 'max(max(kd_smvbf, kd_hmvbf), max(kd_hpvbf, kd_hlvbf))' }

aliases['kd_vbf_hm']    = { 'expr': '1/(1+(me_vbf_hsm/(me_vbf_hm*G4VBF**2)))' }
aliases['kd_vbf_hp']    = { 'expr': '1/(1+(me_vbf_hsm/(me_vbf_hp*G2VBF**2)))' }
aliases['kd_vbf_hl']    = { 'expr': '1/(1+(me_vbf_hsm/(me_vbf_hl*L1VBF**2)))' }
aliases['kd_vbf_mixhm'] = { 'expr': '(me_vbf_mixhm - me_vbf_hsm - me_vbf_hm)/(2*sqrt(me_vbf_hsm*me_vbf_hm))' }
aliases['kd_vbf_mixhp'] = { 'expr': '(me_vbf_mixhp - me_vbf_hsm - me_vbf_hp)/(2*sqrt(me_vbf_hsm*me_vbf_hp))' }

# VH KDs

aliases['pjj_wh'] = { 'expr':'pjjSm_wh/pjjTr_wh' }
aliases['pjj_zh'] = { 'expr':'pjjSm_zh/pjjTr_zh' }

aliases['kd_smwh'] = { 'expr': '1/(1+((me_qcd_hsm*CWH)/(me_wh_hsm*pjj_wh)))' }
aliases['kd_smzh'] = { 'expr': '1/(1+((me_qcd_hsm*CZH)/(me_zh_hsm*pjj_zh)))' }
aliases['kd_smvh'] = { 'expr': 'max(kd_smwh, kd_smzh)' }
aliases['kd_hmwh'] = { 'expr': '1/(1+((me_qcd_hsm*CWH)/(me_wh_hm*pjj_wh*G4WH**2)))' }
aliases['kd_hmzh'] = { 'expr': '1/(1+((me_qcd_hsm*CZH)/(me_zh_hm*pjj_zh*G4ZH**2)))' }
aliases['kd_hmvh'] = { 'expr': 'max(kd_hmwh, kd_hmzh)' }
aliases['kd_hpwh'] = { 'expr': '1/(1+((me_qcd_hsm*CWH)/(me_wh_hp*pjj_wh*G2WH**2)))' }
aliases['kd_hpzh'] = { 'expr': '1/(1+((me_qcd_hsm*CZH)/(me_zh_hp*pjj_zh*G2ZH**2)))' }
aliases['kd_hpvh'] = { 'expr': 'max(kd_hpwh, kd_hpzh)' }
aliases['kd_hlwh'] = { 'expr': '1/(1+((me_qcd_hsm*CWH)/(me_wh_hl*pjj_wh*L1WH**2)))' }
aliases['kd_hlzh'] = { 'expr': '1/(1+((me_qcd_hsm*CZH)/(me_zh_hl*pjj_zh*L1ZH**2)))' }
aliases['kd_hlvh'] = { 'expr': 'max(kd_hlwh, kd_hlzh)' }
aliases['kd_vh']   = { 'expr': 'max(max(kd_smvh, kd_hmvh), max(kd_hpvh, kd_hlvh))' }

aliases['kd_wh_hm']    = { 'expr': '1/(1+(me_wh_hsm/(me_wh_hm*G4WH**2)))' }
aliases['kd_zh_hm']    = { 'expr': '1/(1+(me_zh_hsm/(me_zh_hm*G4ZH**2)))' }
aliases['kd_vh_hm']    = { 'expr': 'max(kd_wh_hm, kd_zh_hm)' }

aliases['kd_wh_hp']    = { 'expr': '1/(1+(me_wh_hsm/(me_wh_hp*G2WH**2)))' }
aliases['kd_zh_hp']    = { 'expr': '1/(1+(me_zh_hsm/(me_zh_hp*G2ZH**2)))' }
aliases['kd_vh_hp']    = { 'expr': 'max(kd_wh_hp, kd_zh_hp)' }

aliases['kd_wh_hl']    = { 'expr': '1/(1+(me_wh_hsm/(me_wh_hl*L1WH**2)))' }
aliases['kd_zh_hl']    = { 'expr': '1/(1+(me_zh_hsm/(me_zh_hl*L1ZH**2)))' }
aliases['kd_vh_hl']    = { 'expr': 'max(kd_wh_hl, kd_zh_hl)' }

aliases['me_vh_hsm']    = { 'expr': '(me_wh_hsm/meAvg_wh) + (me_zh_hsm/meAvg_zh)' }
aliases['me_vh_hm']     = { 'expr': '(me_wh_hm/meAvg_wh) + (me_zh_hm/meAvg_zh)' }
aliases['me_vh_mixhm']  = { 'expr': '((me_wh_mixhm - me_wh_hsm - me_wh_hm)/meAvg_wh) + ((me_zh_mixhm - me_zh_hsm - me_zh_hm)/meAvg_zh)' }
aliases['kd_vh_mixhm'] = { 'expr': '(me_vh_mixhm*G4VH) / (me_vh_hsm + (me_vh_hm*G4VH**2))' }

aliases['me_vh_hp']     = { 'expr': '(me_wh_hp/meAvg_wh) + (me_zh_hp/meAvg_zh)' }
aliases['me_vh_mixhp']  = { 'expr': '((me_wh_mixhp - me_wh_hsm - me_wh_hp)/meAvg_wh) + ((me_zh_mixhp - me_zh_hsm - me_zh_hp)/meAvg_zh)' }
aliases['kd_vh_mixhp'] = { 'expr': '(me_vh_mixhp*G2VH) / (me_vh_hsm + (me_vh_hp*G2VH**2))' }


# Boosted VH KDs
'''
aliases['pjj_Wh'] = { 'expr':'pjjSm_Wh/pjjTr_Wh' }
aliases['pjj_Zh'] = { 'expr':'pjjSm_Zh/pjjTr_Zh' }

aliases['kd_smWh'] = { 'expr': '1/(1+((me_QCD_hsm*CWH)/(me_Wh_hsm*pjj_Wh)))' }
aliases['kd_smZh'] = { 'expr': '1/(1+((me_QCD_hsm*CZH)/(me_Zh_hsm*pjj_Zh)))' }
aliases['kd_smVh'] = { 'expr': 'max(kd_smWh, kd_smZh)' }
aliases['kd_hmWh'] = { 'expr': '1/(1+((me_QCD_hsm*CWH)/(me_Wh_hm*pjj_Wh*G4WH**2)))' }
aliases['kd_hmZh'] = { 'expr': '1/(1+((me_QCD_hsm*CZH)/(me_Zh_hm*pjj_Zh*G4ZH**2)))' }
aliases['kd_hmVh'] = { 'expr': 'max(kd_hmWh, kd_hmZh)' }
aliases['kd_hpWh'] = { 'expr': '1/(1+((me_QCD_hsm*CWH)/(me_Wh_hp*pjj_Wh*G2WH**2)))' }
aliases['kd_hpZh'] = { 'expr': '1/(1+((me_QCD_hsm*CZH)/(me_Zh_hp*pjj_Zh*G2ZH**2)))' }
aliases['kd_hpVh'] = { 'expr': 'max(kd_hpWh, kd_hpZh)' }
aliases['kd_hlWh'] = { 'expr': '1/(1+((me_QCD_hsm*CWH)/(me_Wh_hl*pjj_Wh*L1WH**2)))' }
aliases['kd_hlZh'] = { 'expr': '1/(1+((me_QCD_hsm*CZH)/(me_Zh_hl*pjj_Zh*L1ZH**2)))' }
aliases['kd_hlVh'] = { 'expr': 'max(kd_hlWh, kd_hlZh)' }
aliases['kd_Vh']   = { 'expr': 'max(max(kd_smVh, kd_hmVh), max(kd_hpVh, kd_hlVh))' }

aliases['kd_Wh_hm']    = { 'expr': '1/(1+(me_Wh_hsm/(me_Wh_hm*G4WH**2)))' }
aliases['kd_Zh_hm']    = { 'expr': '1/(1+(me_Zh_hsm/(me_Zh_hm*G4ZH**2)))' }
aliases['kd_Vh_hm']    = { 'expr': 'max(kd_Wh_hm, kd_Zh_hm)' }

aliases['kd_Wh_hp']    = { 'expr': '1/(1+(me_Wh_hsm/(me_Wh_hp*G2WH**2)))' }
aliases['kd_Zh_hp']    = { 'expr': '1/(1+(me_Zh_hsm/(me_Zh_hp*G2ZH**2)))' }
aliases['kd_Vh_hp']    = { 'expr': 'max(kd_Wh_hp, kd_Zh_hp)' }

aliases['kd_Wh_hl']    = { 'expr': '1/(1+(me_Wh_hsm/(me_Wh_hl*L1WH**2)))' }
aliases['kd_Zh_hl']    = { 'expr': '1/(1+(me_Zh_hsm/(me_Zh_hl*L1ZH**2)))' }
aliases['kd_Vh_hl']    = { 'expr': 'max(kd_Wh_hl, kd_Zh_hl)' }

aliases['me_Vh_hsm']    = { 'expr': '(me_Wh_hsm/meAvg_wh) + (me_Zh_hsm/meAvg_zh)' }
aliases['me_Vh_hm']     = { 'expr': '(me_Wh_hm/meAvg_wh) + (me_Zh_hm/meAvg_zh)' }
aliases['me_Vh_mixhm']  = { 'expr': '((me_Wh_mixhm - me_Wh_hsm - me_Wh_hm)/meAvg_wh) + ((me_Zh_mixhm - me_Zh_hsm - me_Zh_hm)/meAvg_zh)' }
aliases['kd_Vh_mixhm'] = { 'expr': '(me_Vh_mixhm*G4VH) / (me_Vh_hsm + (me_Vh_hm*G4VH**2))' }

aliases['me_Vh_hp']     = { 'expr': '(me_Wh_hp/meAvg_wh) + (me_Zh_hp/meAvg_zh)' }
aliases['me_Vh_mixhp']  = { 'expr': '((me_Wh_mixhp - me_Wh_hsm - me_Wh_hp)/meAvg_wh) + ((me_Zh_mixhp - me_Zh_hsm - me_Zh_hp)/meAvg_zh)' }
aliases['kd_Vh_mixhp'] = { 'expr': '(me_Vh_mixhp*G2VH) / (me_Vh_hsm + (me_Vh_hp*G2VH**2))' }

################## Additional variables ##############################
'''

aliases['mV'] = { 'expr': 'FatJet_msoftdrop[CleanFatJet_jetIdx[0]]' }

aliases['j1_px'] = { 'expr': 'CleanJet_pt[0]*cos(CleanJet_phi[0])' }
aliases['j2_px'] = { 'expr': 'CleanJet_pt[1]*cos(CleanJet_phi[1])' }
aliases['j1_py'] = { 'expr': 'CleanJet_pt[0]*sin(CleanJet_phi[0])' }
aliases['j2_py'] = { 'expr': 'CleanJet_pt[1]*sin(CleanJet_phi[1])' }
aliases['ptjj'] = { 'expr': 'sqrt(pow(j1_px[0] + j2_px[0],2) + pow(j1_py[0] + j2_py[0],2))' }

###########################################################################

# ttH_MVA

# variations
eleWP_old = 'mvaFall17V1Iso_WP90'
muWP_old = 'cut_Tight_HWWW'
aliases['SFweightEleUp'] = {
    'expr': 'LepSF3l__ele_'+eleWP+'__Up',
    'samples': mc
}
aliases['SFweightEleDown'] = {
    'expr': 'LepSF3l__ele_'+eleWP+'__Do',
    'samples': mc
}
aliases['SFweightMuUp'] = {
    'expr': 'LepSF3l__mu_'+muWP+'__Up',
    'samples': mc
}
aliases['SFweightMuDown'] = {
    'expr': 'LepSF3l__mu_'+muWP+'__Do',
    'samples': mc
}

aliases['nllWOTF'] = {
    'linesToAdd': ['.L %s/src/PlotsConfigurations/Configurations/Differential/nllW.cc+' % os.getenv('CMSSW_BASE')],
    'class': 'WWNLLW',
    'args': ('central',),
    'samples': ['WW']
}

#######################
### SFs for tthMVA  ###
#######################

aliases['ttHMVA_SF_3l'] = {
    'linesToAdd': ['.L %s/src/PlotsConfigurations/Configurations/patches/compute_SF_BETA.C+' % os.getenv('CMSSW_BASE')],
    'class': 'compute_SF',
    'args' : ('2017', 3, 'total_SF'),
    'samples': mc
}

aliases['ttHMVA_SF_Up_0'] = {
    'class': 'compute_SF',
    'args' : ('2017', 3, 'single_SF_up', 0),
    'nominalOnly' : True,
    'samples': mc
}

aliases['ttHMVA_SF_Up_1'] = {
    'class': 'compute_SF',
    'args' : ('2017', 3, 'single_SF_up', 1),
    'nominalOnly' : True,
    'samples': mc
}

aliases['ttHMVA_SF_Up_2'] = {
    'class': 'compute_SF',
    'args' : ('2017', 3, 'single_SF_up', 2),
    'nominalOnly' : True,
    'samples': mc
}
aliases['ttHMVA_SF_Down_0'] = {
    'class': 'compute_SF',
    'args' : ('2017', 3, 'single_SF_down', 0),
    'nominalOnly' : True,
    'samples': mc
}

aliases['ttHMVA_SF_Down_1'] = {
    'class': 'compute_SF',
    'args' : ('2017', 3, 'single_SF_down', 1),
    'nominalOnly' : True,
    'samples': mc
}

aliases['ttHMVA_SF_Down_2'] = {
    'class': 'compute_SF',
    'args' : ('2017', 3, 'single_SF_down', 2),
    'nominalOnly' : True,
    'samples': mc
}

aliases['ttHMVA_3l_ele_SF_Up'] = {
    'expr' : '(ttHMVA_SF_Up_0[0]*(abs(Lepton_pdgId[0]) == 11) + (abs(Lepton_pdgId[0]) == 13)) *\
              (ttHMVA_SF_Up_1[0]*(abs(Lepton_pdgId[1]) == 11) + (abs(Lepton_pdgId[1]) == 13)) *\
              (ttHMVA_SF_Up_2[0]*(abs(Lepton_pdgId[2]) == 11) + (abs(Lepton_pdgId[2]) == 13))',
    'nominalOnly' : True,
    'samples' : mc
}

aliases['ttHMVA_3l_ele_SF_Down'] = {
    'expr' : '(ttHMVA_SF_Down_0[0]*(abs(Lepton_pdgId[0]) == 11) + (abs(Lepton_pdgId[0]) == 13)) *\
              (ttHMVA_SF_Down_1[0]*(abs(Lepton_pdgId[1]) == 11) + (abs(Lepton_pdgId[1]) == 13)) *\
              (ttHMVA_SF_Down_2[0]*(abs(Lepton_pdgId[2]) == 11) + (abs(Lepton_pdgId[2]) == 13))',
    'nominalOnly' : True,
    'samples' : mc
}

aliases['ttHMVA_3l_mu_SF_Up'] = {
    'expr' : '(ttHMVA_SF_Up_0[0]*(abs(Lepton_pdgId[0]) == 13) + (abs(Lepton_pdgId[0]) == 11)) *\
              (ttHMVA_SF_Up_1[0]*(abs(Lepton_pdgId[1]) == 13) + (abs(Lepton_pdgId[1]) == 11)) *\
              (ttHMVA_SF_Up_2[0]*(abs(Lepton_pdgId[2]) == 13) + (abs(Lepton_pdgId[2]) == 11))',
    'nominalOnly' : True,
    'samples' : mc
}

aliases['ttHMVA_3l_mu_SF_Down'] = {
    'expr' : '(ttHMVA_SF_Down_0[0]*(abs(Lepton_pdgId[0]) == 13) + (abs(Lepton_pdgId[0]) == 11)) *\
              (ttHMVA_SF_Down_1[0]*(abs(Lepton_pdgId[1]) == 13) + (abs(Lepton_pdgId[1]) == 11)) *\
              (ttHMVA_SF_Down_2[0]*(abs(Lepton_pdgId[2]) == 13) + (abs(Lepton_pdgId[2]) == 11))',
    'nominalOnly' : True,
    'samples' : mc
}


# data/MC scale factors
aliases['SFweight'] = {
    'expr': ' * '.join(['SFweight3l', 'ttHMVA_SF_3l', 'LepWPCut', 'btagSF', 'PrefireWeight']),
    'samples': mc
}

aliases['Muon_ttHMVA_SF'] = {
    'expr': '( (abs(Lepton_pdgId[0]) == 13)*(Lepton_tightMuon_cut_Tight_HWWW_tthmva_80_IdIsoSF[0]/Lepton_tightMuon_cut_Tight_HWWW_IdIsoSF[0])+(abs(Lepton_pdgId[0]) == 11) )*( (abs(Lepton_pdgId[1]) == 13)*(Lepton_tightMuon_cut_Tight_HWWW_tthmva_80_IdIsoSF[1]/Lepton_tightMuon_cut_Tight_HWWW_IdIsoSF[1])+ (abs(Lepton_pdgId[1]) == 11) )',
    'samples' : ['Dyemb']
}

'''
aliases['WH3l_dphilllmet_test'] = {
    'linesToAdd': [
        '.L %s/src/PlotsConfigurations/Configurations/WH3l/scripts/WH3l_patch_BDT1718.cc+' % os.getenv('CMSSW_BASE')
    ],
    'class': 'WH3l_patch_BDT1718',
    'args': ("dphilllmet")
}

aliases['WH3l_mOSll_min_test'] = {
    'class': 'WH3l_patch_BDT1718',
    'args': ("mOSllmin")
}

aliases['WH3l_ptOSll_min_test'] = {
    'class': 'WH3l_patch_BDT1718',
    'args': ("ptOSllmin")
}

aliases['WH3l_drOSll_min_test'] = {
    'class': 'WH3l_patch_BDT1718',
    'args': ("drOSllmin")
}

aliases['WH3l_ZVeto_test'] = {
    'class': 'WH3l_patch_BDT1718',
    'args': ("ZVeto")
}

aliases['WH3l_ptlll_test'] = {
    'class': 'WH3l_patch_BDT1718',
    'args': ("ptlll")
}

aliases['WH3l_mtlmet0_test'] = {
    'class': 'WH3l_patch_BDT1718',
    'args': ("mtlmet0")
}

aliases['WH3l_mtlmet1_test'] = {
    'class': 'WH3l_patch_BDT1718',
    'args': ("mtlmet1")
}

aliases['WH3l_mtlmet2_test'] = {
    'class': 'WH3l_patch_BDT1718',
    'args': ("mtlmet2")
}

aliases['WH3l_dphilmet0_test'] = {
    'class': 'WH3l_patch_BDT1718',
    'args': ("dphilmet0")
}

aliases['WH3l_dphilmet1_test'] = {
    'class': 'WH3l_patch_BDT1718',
    'args': ("dphilmet1")
}

aliases['WH3l_dphilmet2_test'] = {
    'class': 'WH3l_patch_BDT1718',
    'args': ("dphilmet2")
}

aliases['WH3l_ptWWW_test'] = {
    'class': 'WH3l_patch_BDT1718',
    'args': ("ptWWW")
}

aliases['WH3l_mtWWW_test'] = {
    'class': 'WH3l_patch_BDT1718',
    'args': ("mtWWW")
}
aliases['WH3l_mlll_test'] = {
    'class': 'WH3l_patch_BDT1718',
    'args': ("mlll")
}

aliases['BDT_SSSF1718'] = {
    'class': 'WH3l_patch_BDT1718',
    'args': ("BDT_SSSF1718")
}
aliases['BDT_SSSFcombin'] = {
    'class': 'WH3l_patch_BDT1718',
    'args': ("BDT_SSSFcombin")
}
aliases['BDT_OSSF1718'] = {
    'class': 'WH3l_patch_BDT1718',
    'args': ("BDT_OSSF1718")
}
aliases['BDT_OSSFcombin'] = {
    'class': 'WH3l_patch_BDT1718',
    'args': ("BDT_OSSFcombin")
}
'''
aliases['Weight2MINLO'] = {
    'linesToAdd': ['.L %s/Differential/weight2MINLO.cc+' % configurations],
    'class': 'Weight2MINLO',
    'args': '%s/src/LatinoAnalysis/Gardener/python/data/powheg2minlo/NNLOPS_reweight.root' % os.getenv('CMSSW_BASE'),
    'samples' : [skey for skey in samples if 'ggH_hww' in skey],
}

# GGHUncertaintyProducer wasn't run for 2017 nAODv5 non-private

thus = [
    'ggH_mu',
    'ggH_res',
    'ggH_mig01',
    'ggH_mig12',
    'ggH_VBF2j',
    'ggH_VBF3j',
    'ggH_pT60',
    'ggH_pT120',
    'ggH_qmtop'
]

for thu in thus:
    aliases[thu+'_2'] = {
        'linesToAdd': ['.L %s/Differential/gghuncertainty.cc+' % configurations],
        'class': 'GGHUncertainty',
        'args': (thu,),
        'samples': mc_ggh,
    }

thusQQ = [
  "qqH_YIELD",
  "qqH_PTH200",
  "qqH_Mjj60",
  "qqH_Mjj120",
  "qqH_Mjj350",
  "qqH_Mjj700",
  "qqH_Mjj1000",
  "qqH_Mjj1500",
  "qqH_PTH25",
  "qqH_JET01",
  "qqH_EWK",
]

for thu in thusQQ:
    aliases[thu] = {
        'linesToAdd': ['.L %s/patches/qqhuncertainty.cc+' % configurations],
        'class': 'QQHUncertainty',
        'args': (thu,),
        'samples': mc_qqh,
        'nominalOnly': True
    }


