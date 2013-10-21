# qsub_run_samtools.sh
qsub -N myjobtest -e ./logs -o ./logs -l mem=8G,time=2:: -S /bin/sh -cwd run_samtools.sh