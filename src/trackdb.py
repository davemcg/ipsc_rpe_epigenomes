# python3

import sys
name = sys.argv[1].split(',')
big_wigs = sys.argv[2].split(',')
parents = sys.argv[3].split(',')
big_wig_colors = sys.argv[4].split(';')
bigBeds = sys.argv[5].split(',')
bigBed_colors = sys.argv[6].split(';')


def make_container(name):
	out = 'track ' + name + '\n'
	out += 'container multiWig\n'
	out += 'shortLabel ' + name + '\n'
	out += 'longLabel ' + name + '\n'
	out += 'type bigWig 0 30000\n\
viewLimits 0:50\n\
visibility full\n\
maxHeightPixels 150:90:11\n\
aggregate transparentOverlay\n\
showSubtrackColorOnUi on\n\
windowingFunction mean\n\
priority 1.4\n\
configurable on\n\
autoScale off\n\n'
	print(out)

for element in name:
	make_container(element)

def track_maker(file_name, parent, color, track_type):
	out = 'track ' + file_name + '\n'
	if parent != "None":
		out += 'parent ' + parent + '\n'
	out += 'color ' + color + '\n'
	out += 'bigDataUrl ' + file_name + '\n'
	out += 'shortLabel ' + file_name.split('.')[0] + '\n'
	out += 'longLabel ' + file_name.split('.')[0] + '\n'
	out += 'type ' + track_type + '\n\n'
	print(out)

for (name, parent, color) in zip(big_wigs, parents, big_wig_colors):
	track_maker(name, parent, color, 'bigWig')

for (name, color) in zip(bigBeds, bigBed_colors):
	track_maker(name, "None", color, 'bigBed')
