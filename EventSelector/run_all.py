from sys import argv, stdout
from ROOT import *
from subprocess import call, Popen
from time import sleep

def main():
    if len(argv) < 4:
        print "Usage:", argv[0], "outfile_dir logfile_dir whitespace [separated infiles]"
        return

    # call("pod-server start", shell=True)
    call("make", shell=True)

    outfile_dir = argv[1].rstrip('/') + '/'
    logfile_dir = argv[2].rstrip('/') + '/'
    fnames = argv[3:]
    procs = []
    for fname in fnames:
        outfname = '.'.join(fname.split("/")[-1].split(".")[:-1])
        rootfname = outfile_dir + outfname + '.out.root'
        logfname = logfile_dir + outfname + '.log'

        print "nice ./run %s %s >& %s" % (fname, rootfname, logfname)
        stdout.flush()

        procs.append(Popen("nice ./run %s %s >& %s" % (fname, rootfname,
            logfname), shell=True))

        sleep(30)

if __name__=='__main__':
    main()
