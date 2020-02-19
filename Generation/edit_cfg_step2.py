import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--filename")
args = parser.parse_args()
filename = args.filename

with open(filename, 'r') as file:
  cfg_file = file.read()

cfg_file = cfg_file.replace("'$inputFileNames'", "$inputFileNames")

with open(filename, 'w') as file:
  file.write(cfg_file)
