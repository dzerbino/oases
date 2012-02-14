#!/usr/bin/env python

import sys
import subprocess
import optparse
import shutil
import os

##########################################
## Options and defaults
##########################################
def getOptions():
	parser = optparse.OptionParser('usage: %prog [options] --data "velveth file descriptors"')
	parser.add_option('-d', '--data',dest='data',help='Velveth file descriptors',metavar='FILES', default='')
	parser.add_option('-p', '--options',dest='oasesOptions',help='Oases options that are passed to the command line, e.g., -cov_cutoff 5 -ins_length 300 ',metavar='OPTIONS', default='')
	parser.add_option('-m', '--kmin',dest='kmin',type="int",help='Minimum k',default=19)
	parser.add_option('-M', '--kmax',dest='kmax',type="int",help='Maximum k',default=31)
	parser.add_option('-s', '--kstep',dest='kstep',type="int",help='Steps in k',default=2)
	parser.add_option('-g', '--merge',dest='kmerge',type="int",help='Merge k',default=27)
	parser.add_option('-o', '--output',dest='directoryRoot',help='Output directory prefix',metavar='NAME',default='oasesPipeline')
	parser.add_option('-r', '--mergeOnly',dest='mergeOnly',help='Only do the merge',action='store_true',default=False)
	parser.add_option('-c', '--clean',dest='clean',help='Clean temp files',action='store_true',default=False)
	options, args = parser.parse_args()
	if not options.mergeOnly and len(options.data) == 0:
		parser.print_help()
		print ''
		print 'You forgot to provide some data files!'
		print 'Current options are:'
		print options
		sys.exit(1)
	return options

##########################################
## Assembly procedure
##########################################
def singleKAssemblies(options):
	p = subprocess.Popen(['velveth', options.directoryRoot, '%i,%i,%i' % (options.kmin, options.kmax+options.kstep, options.kstep)] + options.data.split(), stdout=subprocess.PIPE)
	output = p.communicate()
	assert p.returncode == 0, output[0] + "Hash failed\n"
	for k in range(options.kmin, options.kmax+options.kstep, options.kstep):
	    p = subprocess.Popen(['velvetg','%s_%i' % (options.directoryRoot, k), '-read_trkg', 'yes'], stdout=subprocess.PIPE)
	    output = p.communicate()
	    assert p.returncode == 0, "Velvetg failed at k = %i\n%s" % (k, output[0])
	    p = subprocess.Popen(['oases','%s_%i' % (options.directoryRoot, k)] + options.oasesOptions.split(), stdout=subprocess.PIPE)
	    output = p.communicate()
	    assert p.returncode == 0, "Oases failed at k = %i\n%s" % (k, output[0])
	
def mergeAssemblies(options):
	files = ["%s_%i/transcripts.fa" % (options.directoryRoot, X) for X in range(options.kmin, options.kmax+options.kstep, options.kstep)]
	p = subprocess.Popen(['velveth','%sMerged' % options.directoryRoot, str(options.kmerge), '-long'] + files, stdout=subprocess.PIPE)
	output = p.communicate()
	assert p.returncode == 0, output[0] + "Velveth failed at merge\n"
	p = subprocess.Popen(['velvetg','%sMerged' % options.directoryRoot,'-conserveLong','yes','-read_trkg','yes'], stdout=subprocess.PIPE)
	output = p.communicate()
	assert p.returncode == 0, output[0] + "Velvetg failed at merge\n"
	p = subprocess.Popen(['oases','%sMerged' % options.directoryRoot,'-merge','yes'], stdout=subprocess.PIPE)
	output = p.communicate()
	assert p.returncode == 0, output[0] + "Oases failed merge\n"

##########################################
## Checking dependencies
##########################################
def checkVelvet():
	try:
	    p = subprocess.Popen(['velveth'], stdout=subprocess.PIPE)
	except OSError:
	    print "Could not find Velvet"
	    print "Make sure that it is properly installed on your path"
	    sys.exit(1)
	for line in p.stdout:
		items = line.strip().split(' ')
		if items[0] == 'Version':
			items2 = map(int, items[1].split('.'))
			assert items2 >= [1,1,7], "Velvet must have version 1.1.07 or higher (currently %s)" % items[1]
			return
	assert False

def checkOases():
	try:
	    p = subprocess.Popen(['oases'], stdout=subprocess.PIPE)
	except OSError:
	    print "Could not find Oases"
	    print "Make sure that it is properly installed on your path"
	    sys.exit(1)
	for line in p.stdout:
		items = line.strip().split(' ')
		if items[0] == 'Version':
			items2 = map(int, items[1].split('.'))
			assert items2 >= [0,2,1], "Oases must have version 0.2.01 or higher (currently %s)" % items[1]
			return
	assert False

##########################################
## Clean up
##########################################
def clean(options):
	for k in range(options.kmin, options.kmax, options.kstep):
	    shutil.rmtree("%s_%i" % (options.directoryRoot, k))
	os.remove("%sMerged/Sequences" % options.directoryRoot)
	os.remove("%sMerged/Roadmaps" % options.directoryRoot)
	os.remove("%sMerged/PreGraph" % options.directoryRoot)
	os.remove("%sMerged/Graph2" % options.directoryRoot)
	os.remove("%sMerged/LastGraph" % options.directoryRoot)
	os.remove("%sMerged/contigs.fa" % options.directoryRoot)
	os.remove("%sMerged/Log" % options.directoryRoot)

##########################################
## Master function
##########################################
def main():
	options = getOptions()
	checkVelvet()
	checkOases()
	if not options.mergeOnly:
	    singleKAssemblies(options)
	mergeAssemblies(options)
	if options.clean:
	    clean(options)

if __name__ == "__main__":
	main()
