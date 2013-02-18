#!/usr/bin/python
from glob import glob
from subprocess import Popen
from sys import stdout, argv
from time import sleep

def main():
    if len(argv) < 4:
        print "Usage:", argv[0], "outfile_dir logfile_dir whitespace [separated cmd files]"
        return

    outfile_dir = argv[1].rstrip('/') + '/'
    logfile_dir = argv[2].rstrip('/') + '/'
    fnames = argv[3:]
    procs = []
    for fname in fnames:
        outfname = ".".join(fname.split('/')[-1].split('.')[:-1])
        rootfname = outfile_dir + outfname + '.root'
        logfname = logfile_dir + outfname + '.log'

        print "./run %s %s >& %s" % (fname, rootfname, logfname)
        stdout.flush()

        procs.append(Popen("./run %s %s >& %s" % (fname, rootfname,
            logfname), shell=True))

        sleep(10)

if __name__=='__main__':
    main()
