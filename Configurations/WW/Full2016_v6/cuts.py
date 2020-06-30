# cuts



_tmp = [ 
     'Lepton_pdgId[0]*Lepton_pdgId[1] == -11*13',
     'mll>12.',
     'Lepton_pt[0]>25.',
     'Lepton_pt[1]>20.',
     '(nLepton>=2 && Alt$(Lepton_pt[2],0)<10.)',
     'PuppiMET_pt > 20.',
     'mpmet > 20.',
     'ptll > 30.',
       ]

supercut = ' && '.join(_tmp)

def addcut(name, exprs):
    cuts[name] = ' && '.join(exprs)


# Jet_btagDeepB

_tmp = [
    'mth > 40.',
    'mll < 76.',
    'drll < 2.5',
    'bVeto',
       ]

addcut('SR_Incl', _tmp)


_tmp = [
    'mth > 40.',
    'mll > 76.',
    'drll < 2.5',
    'bVeto',
       ]

addcut('WWCR_Incl', _tmp)



_tmp = [
    'mth > 40.',
    'mll < 76.',
    'drll < 2.5',
    '(bReq || (!bVeto && zeroJet))',
       ]

addcut('TopCR_Incl', _tmp)

