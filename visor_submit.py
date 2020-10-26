import os, sys
from optparse import OptionParser


def submit_jobs():
    purities = [100, 80, 60, 40, 20]

    tumour_id = 1
    tumour_coverage=30
    normal_coverage=50

    for p in purities:
        tid = tumour_id
        nid = tumour_id + 1
        tumour_sample = "visorR%s" % tid
        normal_sample = "visorR%s" % nid

        t_cmd = "qsub -v VAR1=%s,VAR2=%s,VAR3=%s,VAR4=%s -o visor/log/%s.runlog -j oe -N %s.visor run_visor.pbs" % (tumour_sample, p, 'tumour', tumour_coverage, tumour_sample, tumour_sample)
        n_cmd = "qsub -v VAR1=%s,VAR2=%s,VAR3=%s,VAR4=%s -o visor/log/%s.runlog -j oe -N %s.visor run_visor.pbs" % (normal_sample, p, 'normal', normal_coverage, normal_sample, normal_sample)

        print(t_cmd)
        print(n_cmd)
        os.system(t_cmd)
        os.system(n_cmd)

        tumour_id += 2


submit_jobs()
