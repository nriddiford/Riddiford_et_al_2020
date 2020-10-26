import os, sys
from optparse import OptionParser


def submit_jobs():
    purities = [100, 80, 60, 40, 20]

    tumour_id = 1

    for p in purities:
        tid = tumour_id
        nid = tumour_id + 1
        tumour_sample = "visorR%s" % tid
        normal_sample = "visorR%s" % nid

        cmd = "qsub -v VAR1=%s,VAR2=%s,VAR3=%s -o visor/log/%s.runlog -j oe -N %s.visor run_visor.pbs" % (tumour_sample, normal_sample, p, tumour_sample, tumour_sample)

        print(cmd)
        os.system(cmd)

        tumour_id += 2


submit_jobs()