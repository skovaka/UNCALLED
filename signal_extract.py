import h5py
import sys
import numpy as np
import argparse
import os

parser = argparse.ArgumentParser(description='Get raw electrical signal from a fast5')
parser.add_argument('--input', '-i', type=str, required=True, help='input location fast5 file')
args=parser.parse_args()

hdf=h5py.File(os.path.abspath(args.input),'r')
readnum=args.input.split('_')[-2]
readname='R'+readnum[1:4]+'_'+readnum[4:]

sig=np.asarray(hdf['Raw/Reads/'+readname+'/Signal'])

np.savetxt(readnum+'_sig.txt', sig, fmt='%i', delimiter='\n')

