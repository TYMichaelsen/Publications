rom BCBio import GFF
import pysam
import sys
import csv
import os
from argparse import ArgumentParser

### INPUT ARGUMENTS ###
parser = ArgumentParser()

parser.add_argument("--bam",
  metavar = "STRING",
  type    = str,
  help    = ".bam file")

parser.add_argument("--gff",
  metavar = "STRING",
  type    = str,
  help    = ".gff file")

parser.add_argument("--out",
  metavar = "STRING",
  type    = str,
  help    = "Output directory",
  default = ".")

args = parser.parse_args()

filBAM = args.bam
name   = os.path.splitext(os.path.basename(filBAM))[0]
filGFF = args.gff

# Open the gff file.
gffHandle = GFF.parse(open(filGFF))
# Open the bam file
bamHandle = pysam.AlignmentFile(filBAM, "rb")

### Functions for read directions ###
def Forward(read):
	if read.is_reverse: return False
	else: return True

def Reverse(read):
	if read.is_reverse: return True
	else: return False

### CALCULATE COVERAGE AND ANNOTATION FOR EVERY LOCUS_TAG ###
coverages=list()

for contig in gffHandle:
	if contig.id in bamHandle.references:

		for feature in contig.features:

			if 'locus_tag' in feature.qualifiers:
				qual = 'locus_tag'
			else:
				if 'ID' in feature.qualifiers:
					qual = 'ID'
				else: 
					continue

			tag      = feature.qualifiers[qual][0]
			f_strand = feature.strand
			f_start  = feature.location.start
			f_stop   = feature.location.end
			f_length = f_stop - f_start

			# Def function: Count reads according to sense and antisense.
			if f_strand == 1:
				sense     = bamHandle.count(contig.id,start = f_start,stop = f_stop,read_callback = Forward)
				antisense = bamHandle.count(contig.id,start = f_start,stop = f_stop,read_callback = Reverse)
			if f_strand == -1:
				sense     = bamHandle.count(contig.id,start = f_start,stop = f_stop,read_callback = Reverse)
				antisense = bamHandle.count(contig.id,start = f_start,stop = f_stop,read_callback = Forward)

			# Store in object.
			coverages.append((name,tag,str(contig.id),f_strand,f_start,f_length,sense,antisense))
			
### PRINT OBJECT AS CSV ###
with open(args.out+'/'+name+'.csv', "w") as f:
  csv_out = csv.writer(f)
  csv_out.writerows(coverages)