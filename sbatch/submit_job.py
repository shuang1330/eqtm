import os
import shlex, subprocess

def bash_command(cmd):
    subprocess.Popen(['/bin/bash', '-c', cmd])

jobs_path = './jobList'
for job in os.listdir(jobs_path):
    bash_command('sbatch {}'.format(job))
    print('Submitted Job:',job)
print('Check submitted jobs:')
bash_command('squeue -u umcg-sli')
