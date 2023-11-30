## HPC Notes

# vpn login
vpn.arizona.edu

ssh -X dannyjackson@hpc.arizona.edu

# after login, type "shell" and hit return to access the job submission hosts for all environments
# otherwise, you're on a "bastion host" which is a secure portal from "the outside and serves no computational function" (literally why is this involved?)
# this only gets you to a login node, which are not computational nodes.
# you can submit jobs, request interactive sessions, track jobs, and move/edit files from here

# there are 3 clusters: puma, ocelot, and elgato
# switch between them by using their name as a command

# upon opening, it says "dannyjackson@wentletrap" for some reason
# i'm in /home/u15/dannyjackson/ which contains a directory "ondemand" which has a bunch of stuff but none of it makes sense or has any files. the deepest directory is titled, "UAz_vscode" which makes me think this exists because I added the vscode interactive stuff through the HPC website? idk

# use /tmp for working space
# to check usage:
uquota

# only Sabrina and I have full access to the xdisk. I've asked the HPC folks to add a subdirectory that only I have access to (and Sabrina) so that I have a space to work within but that the undergrads can't affect. 

# per month, and refreshed on the first of the month, each group gets:
7,000 CPU hrs on elgato
70,000 on ocelote
100,000 on puma

# puma has a lot of hours but it's chockablock full at the moment and seems to always be. ocelote is also full. i think my workflow should be to test scripts on elgato and run them on one of the other two.

# to check on allocation:
va

# jobs can only run for a max of 10 days

# to submit an interactive, first login to the desired node (elgato works fastest) then use:
interactive -a mcnew

# without the account option, it goes into windfall which takes a lot longer

## Batch jobs 
# header

#!/bin/bash

### define job name
#SBATCH --job-name=XXX

### define output file
### %x is jobname and #j is job ID
#SBATCH --output=%x-%j.out
### REQUIRED specify PI
#SBATCH --account=mcnew
### REQUIRED set partition
#SBATCH --partition=standard
### REQUIRED set number of nodes
#SBATCH --nodes=1
### REQUIRED set number of CPUs
#SBATCH --ntasks=1
### REQUIRED set memory
#SBATCH --mem-per-cpu=5gb
### REQUIRED specify time, hhh:mm:ss
#SBATCH --time=00:01:00


# submit job with
sbatch jobname.slurm

# check on job status with
squeue --job #jobID

