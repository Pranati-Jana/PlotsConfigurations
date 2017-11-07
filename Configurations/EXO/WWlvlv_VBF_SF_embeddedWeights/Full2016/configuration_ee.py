# example of configuration file

tag = 'HWWhighMass2016_SF_embeddedWeights_ee_VBF_up200'


#ee 1nw, mm 2nw

# used by mkShape to define output directory for root files
outputDir = 'rootFile_SF_ee_VBF_up200'


# file with list of variables
variablesFile = 'variables_SF.py'

# file with list of cuts
cutsFile = 'cuts_ee.py' 

# file with list of samples
samplesFile = 'samples_SF_ee.py' 

# file with list of samples
plotFile = 'plot_SF.py' 



# luminosity to normalize to (in 1/fb)
# lumi = 2.264
#lumi = 6.264

#RIMETTERE baseW in DY e top
#lumi=1.0

lumi = 35.9

# used by mkPlot to define output directory for plots
# different from "outputDir" to do things more tidy
outputDirPlots = 'plotHWWhighMass_SF_ee_VBF_up200'


# used by mkDatacards to define output directory for datacards
outputDirDatacard = 'datacards_SF_ee_VBF'


# structure file for datacard
structureFile = 'structure_SF_merge.py'


# nuisances file for mkDatacards and for mkShape
nuisancesFile = 'nuisances_SF_ee.py'

#nuisancesFile = 'nuisances_VUOTO.py'