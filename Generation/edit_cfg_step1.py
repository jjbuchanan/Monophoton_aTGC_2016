import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--filename")
args = parser.parse_args()
filename = args.filename

pileupinput_filename = '/afs/hep.wisc.edu/cms/jjbuchanan/aTGC_generation/CMSSW_9_4_0_patch1/src/Neutrino_E-10_gun_GSDR_files.txt'
with open(pileupinput_filename, 'r') as file:
  pileupinput_sources = file.read()

pileupinput_sources = str(pileupinput_sources.strip().split())
pileupinput_sources = pileupinput_sources.replace("', '", "',\n'")

with open(filename, 'r') as file:
  cfg_file = file.read()

cfg_file = cfg_file.replace("'$inputFileNames'", "$inputFileNames")
cfg_file = cfg_file.replace("['PILEUPINPUT']", pileupinput_sources)

with open(filename, 'w') as file:
  file.write(cfg_file)
