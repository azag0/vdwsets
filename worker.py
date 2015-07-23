#!/usr/bin/env python
import sys
import os
import glob
import time
import subprocess
import signal

nargs = len(sys.argv[1:])
prefix, myid = sys.argv[1:3]
scratch = sys.argv[3] if len(sys.argv[1:]) == 3 else None
os.chdir(prefix)


def handler_sigint(sig, frame):
        print('Worker {} interrupted, aborting.'.format(myid))
        sys.exit()
signal.signal(signal.SIGINT, handler_sigint)

print('Worker {} alive and ready.'.format(myid))
while True:
    tasks = glob.glob('*.start')
    if not tasks:
        break
    startname = tasks[0]
    basename = os.path.splitext(startname)[0]
    runname = basename + '.running.' + myid
    try:
        os.rename(startname, runname)
    except:
        continue
    print('Worker {} started working on {}...'.format(myid, basename))
    if scratch:
        today = time.strftime('%y-%m-%d')
        rundir = os.path.join(scratch, today, myid, basename)
        os.makedirs(rundir)
        os.symlink(rundir, os.path.join(runname, 'rundir'))
    else:
        rundir = os.path.join(runname, 'rundir')
        os.makedirs(rundir)
    subprocess.call('rsync -a --exclude=rundir ./%s/ %s' % (runname, rundir),
                    shell=True)
    subprocess.call('cd %s && ./run >run.log 2>run.err' % rundir, shell=True)
    os.rename(runname, basename + '.done')
    print('Worker {} finished working on {}.'.format(myid, basename))
print('Worker {} has no more tasks to do, aborting.'.format(myid))
