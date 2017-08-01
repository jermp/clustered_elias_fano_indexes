#!/usr/bin/python

import sys
import gzip
import struct

input_file = sys.argv[1] # collection filename
split_file = input_file.split('/')
file = split_file[len(split_file) - 1]

print "Processing file: ", file

# discard all plists whose cumulative contribution to the number of total ints is
# < P percent of total ints, e.g., P = 10 for 10%
P = float(sys.argv[2])
P = P / 100

max_length = int(sys.argv[3]) # max postings list length

input = open(input_file, 'rb')

def extract_pos_len(input):
	lines = []
	pos = 0
	total_ints = 0

	while True:
		bytes = input.read(4)
		if bytes == "":
			break
		 # 'I' stands for unsigned int 4 bytes
		length = struct.unpack('I', bytes)[0]
		total_ints += length
		lines.append([pos, length])
		pos += length + 1
		input.seek(4 * pos)
		
	return [lines, total_ints]


lines_total_ints = extract_pos_len(input)
print "positions and lengths extracted"

lines = lines_total_ints[0]
total_ints = lines_total_ints[1]

print "sorting lengths"
lines.sort(key = lambda x : x[1])

excluded_ints = 0
saved_ints = 0
plists = len(lines)
saved_plists = 0
excluded_plists = 0

saved_plists_positions_list = []
excluded_plists_positions_list = []

print "filtering plists"
for i in xrange(0, plists):
	
	line = lines[i]
	pos = line[0]
	size = line[1]

	if float(excluded_ints) / total_ints > P and size <= max_length:
		saved_ints += size
		saved_plists += 1
		saved_plists_positions_list.append(pos)
	else:
		excluded_ints += size
		excluded_plists += 1
		excluded_plists_positions_list.append(pos)

print "total postings:", total_ints
print "saved postings:", saved_ints
print "percentage of saved ints: " + "{:1.2f}".format(saved_ints * 100.0 / total_ints)

print "total number of lists:", plists
print "number of saved lists:", saved_plists
print "number of excluded lists:", excluded_plists

print "sorting positions and writing to files"

def sort_and_write(list, file):
	list.sort()
	for x in list:
		file.write(str(x) + '\n');
	file.close()

saved_plists_positions = gzip.open(file + ".plists_positions.gz", 'w')
excluded_plists_positions = gzip.open(file + ".excluded_plists_positions.gz", 'w')
sort_and_write(saved_plists_positions_list, saved_plists_positions)
sort_and_write(excluded_plists_positions_list, excluded_plists_positions)
