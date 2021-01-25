import sys
import json
import csv
from Bio.Seq import Seq

if not len(sys.argv) == 4:
	sys.stdout.write('Usage: python s4_demux_check [Stats.json] [Initial sample-sheet.csv] [Output sample-sheet.csv]' + '\n' + '\n')
	exit()


stats_file = sys.argv[1]
sample_sheet = sys.argv[2]
out_sheet = open(sys.argv[3], "w")


with open(stats_file) as f: #read in the Stats.json
	data = json.load(f)
f.close()

laneID = data['ReadInfosForLanes'][0]['LaneNumber'] #what lane are we in

unknowns_dict = data['UnknownBarcodes'][0]['Barcodes'] #easier to call dict of unknown indices and their read numbers from json file

sample_dict = {} #a dict of the demux results from the json file
for i in range(len(data['ConversionResults'][0]['DemuxResults'])):
	sample_dict[data['ConversionResults'][0]['DemuxResults'][i]['SampleId']] = [data['ConversionResults'][0]['DemuxResults'][i]['IndexMetrics'][0]['IndexSequence'],data['ConversionResults'][0]['DemuxResults'][i]['NumberReads']]

with open(sample_sheet) as sheet1: #add the project name metadata to each sampleID
	sheet1reader = csv.reader(sheet1)
	for row in sheet1reader:
		if not row[0].isdigit():
			continue
		if int(row[0]) == laneID:
			sample_dict[row[1]].append(row[5])
sheet1.close()


s = "+" #we'll use this in join functions below

def i2revcomper(first_index): #this will revcomp i2
	indices = first_index.split("+")
	I2 = Seq(indices[1])
	revcomp_I2 = I2.reverse_complement()
	fix_list = (indices[0], str(revcomp_I2))
	fixed_index = s.join(fix_list)
	return(fixed_index)

def i1revcomper(first_index): # this will revcomp i1
	indices = first_index.split("+")
	I1 = Seq(indices[0])
	revcomp_I1 = I1.reverse_complement()
	fix_list = (str(revcomp_I1), indices[1])
	fixed_index = s.join(fix_list)
	return(fixed_index)

def iswapper(first_index): #this will swap the positions of i1 and i2
	indices = first_index.split("+")
	fix_list = (indices[1], indices[0])
	fixed_index = s.join(fix_list)
	return(fixed_index)

def reporter(x1, x2, x3): #this will tell us if a change is made
	out = "Found an issue with " + x1 + " - changing  " + x2 + "  to  " + x3
	sys.stdout.write(out + '\n')

#write headers to the sample sheet
out_sheet.write('[Data],,,,,' + '\n')
out_sheet.write('Lane,Sample_Name,Sample_ID,Index,Index2,Sample_Project' + '\n')

for key in sample_dict: #we're going to iterate through the demuxed  samples

	index1 = sample_dict[key][0] #this is the original index
	index_split = index1.split("+")
	fixed_sheet = [laneID, key, key, index_split[0], index_split[1], sample_dict[key][2]] #setting up the new sample sheet line
	
	#mutating the indices 
	index_rci2 = i2revcomper(index1)
	index_rci1 = i1revcomper(index1)
	index_iswap = iswapper(index1)
	
	#are the mutations in the unknowns?
	if index_rci2 in unknowns_dict:
		if unknowns_dict[index_rci2] > sample_dict[key][1]: #check to make sure there's more of the unknown that what we had to start with
			index_rci2_split = index_rci2.split("+")
			fixed_sheet[4] = index_rci2_split[1]
			reporter(key, index1, index_rci2)
	elif index_rci1 in unknowns_dict:
		if unknowns_dict[index_rci1] > sample_dict[key][1]:
			index_rci1_split = index_rci1.split("+")
			fixed_sheet[3] = index_rci2_split[0]
			reporter(key, index1, index_rci1)
	elif index_iswap in unknowns_dict:
		if unknowns_dict[index_iswap] > sample_dict[key][1]:
			index_iswap_split = index_iswap.split("+")
			fixed_sheet[3] = index_iswap_split[0]
			fixed_sheet[4] = index_iswap_split[1]
			reporter(key, index1, index_iswap)
	
	row = ",".join([str(item) for item in fixed_sheet])
	out_sheet.write(row + '\n')

out_sheet.close()		





	

