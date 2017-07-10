#!/usr/local/bin/python3

import os
import sys
import argparse
import numpy as np 

import magic

class ArgumentParserError(Exception):
	pass

class NewArgumentParser(argparse.ArgumentParser):
	def error(self, message):
		print(message)
		sys.exit(0)

def parse_args(args):
	p = NewArgumentParser(description='run MAGIC')
	p.add_argument('filetype',
                   choices=['csv', '10x', 'mtx'],
                   help='what is the file type of your original data?')

	a = p.add_argument_group('required arguments')
	a.add_argument('-d', '--data-file', metavar='D', required=True,
                   help='File path of input data file.')
	a.add_argument('-o', '--output-file', metavar='O', required=True,
				   help='File path of where to save the MAGIC imputed data (in csv format).')
	
	n = p.add_argument_group('normalization parameters')
	n.add_argument('-n', '--no-normalize', default=True, action='store_false',
					help='Do not perform library size normalization on the data')
	n.add_argument('-l', '--log-transform', metavar='L', default=None,
					help='Log-transform data with the specified pseudocount.')
	m = p.add_argument_group('MAGIC parameters')
	m.add_argument('-g', '--gene-name-file', metavar='G', 
					help='Gene name file must be specified when loading mtx data.')
	m.add_argument('-p', '--pca-components', metavar='P', default=20,
				   help='Number of pca components to use when running MAGIC (Default = 20).')
	m.add_argument('--pca-non-random', default=True, action='store_false',
				    help='Do not used randomized solver in PCA computation.')
	m.add_argument('-t', metavar='T', default=6, 
					help='t parameter for running MAGIC (Default = 6).')
	m.add_argument('-k', metavar='K', default=30, 
					help='Number of nearest neighbors to use when running MAGIC (Default = 30).')
	m.add_argument('-ka', metavar='KA', default=10, 
					help='knn-autotune parameter for running MAGIC (Default = 10).')
	m.add_argument('-e', '--epsilon', metavar='E', default=1, 
					help='Epsilon parameter for running MAGIC (Default = 1).')
	m.add_argument('-r', '--rescale', metavar='R', default=99, 
					help='Percentile to rescale data to after running MAGIC (Default = 99).')

	try:
		return p.parse_args(args)
	except ArgumentParserError:
		raise


def main(args: list = None):
	args = parse_args(args)
	try:
		if args.filetype == 'csv':
			scdata = magic.mg.SCData.from_csv(os.path.expanduser(args.data_file), 
											  data_type='sc-seq', normalize=args.no_normalize)
		elif args.filetype == 'mtx':
			scdata = magic.mg.SCData.from_mtx(os.path.expanduser(args.data_file),
											  os.path.expanduser(args.gene_name_file), normalize=args.no_normalize)
		elif args.filetype == '10x':
			scdata = magic.mg.SCData.from_10x(args.data_file, normalize=args.no_normalize)

		if args.log_transform != None:
			scdata.log_transform_scseq_data(pseudocount=args.log_transform)

		scdata.run_magic(n_pca_components=args.pca_components, random_pca=args.pca_non_random, t=args.t,
						 k=args.k, ka=args.ka, epsilon=args.epsilon, rescale_percent=args.rescale)

		scdata.save_magic_to_csv(os.path.expanduser(args.output_file))

	except:
		raise

if __name__ == '__main__':
	main(sys.argv[1:])