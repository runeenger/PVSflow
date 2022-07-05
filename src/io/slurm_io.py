################################
### Module to write the slurm files for the supercomputer
##################################


def get_slurmtemplate (jobname,templatefile='../src/io/template.slurm'):
    
    
    with open(templatefile, 'r') as f:

        slurm_template=[]
        slurm_template+='#!/bin/bash\n'
        slurm_template+='# Job name:\n'
        slurm_template+='#SBATCH --job-name='+jobname+'\n'
        slurm_template+='#\n'
        slurm_template+='# Project:\n'
        slurm_template+='#SBATCH --account=NN9279K\n'
        slurm_template+='#\n'
        slurm_template+='# Wall clock limit:\n'
        slurm_template+='#SBATCH --time=50:00:00\n'
        slurm_template+='#\n'
        slurm_template+='#SBATCH --partition=bigmem\n'
        slurm_template+='#\n'
        slurm_template+='# Memory per CPU:\n'
        slurm_template+='#SBATCH --mem-per-cpu=10G\n'
        slurm_template+='#\n'
        slurm_template+='# Number of processes:\n'
        slurm_template+='#SBATCH --ntasks=1\n'
        slurm_template+='#\n'

        slurm_template.extend(f.readlines())

        
        return slurm_template



    
def base_commandline(lpvs=200e-4, c0init='constant', c0valuePVS=50, c0valueSAS=0, sigma=1e-4, sasbc='scenarioB', tend=800, toutput=1, dt=1e-3, r=-1, nr=8, nl=100, d=2e-7, refineleft=True, sas=False) :
    """
    Create a base for the command line

    Parameters
    ----------
    lpvs : FLOAT, optional
        length of the vessel. The default is 200e-4.
    c0init : STRING,optional
        kind of initialisation for the concentration . 
        The default is 'null' which means the initial concentration is zero everywhere in the PVS. 
        If set to 'constant' the concentration value in the PVS will be c0value.
    c0value : FLOAT, optional
        Initial value for the concentration in the SAS and in the PVS depending on c0init. The default is 50.
    sasbc : STRING, optional
        Kind of boundary condition in the SAS. 
        The default is 'scenarioB' : mass conservation, no outflow of CSF out of the SAS.
        'scenarioA' : zero concentration in the SAS - the tracer is flushed out.
        'scenarioB': massconservation no outflow
        'scenarioC' : constant CSF outflow.
        'scenarioD' : pressure dependent outflow.
    tend : FLOAT, optional
        end time of the simulation. The default is 800.
    toutput : FLOAT, optional
        period for the outputs. The default is 1.
    dt : FLOAT, optional
        time step. The default is 1e-3.
    r : FLOAT, optional
        Flow resistance at the right boundary. 
        The default is -1 which leads to a zero flow condition.
    nr : INT, optional
        Number of cells in the radial direction. The default is 8.
    nl : INT, optional
        number of cells in the longitudinal direction. The default is 100.
    d : FLOAT, optional
        Modelecular diffusion coefficient. The default is 2e-7.
    refineleft : BOOL, optional
        State if the mesh must have a geometrical progression of the left side. 
        The default is True.

    Returns
    -------
    STRING : base for the command line

    """
    

    
    jobcommand='srun -n 1 python3 PVS_simulation.py'

    jobcommand+=' -lpvs '+str(lpvs)
    jobcommand+=' -c0init '+c0init
    jobcommand+=' -c0valueSAS '+str(c0valueSAS)
    jobcommand+=' -c0valuePVS '+str(c0valuePVS)
    jobcommand+=' -sasbc '+sasbc
    jobcommand+=' -tend '+str(tend)
    jobcommand+=' -toutput '+str(toutput)
    jobcommand+=' -dt '+str(dt)
    jobcommand+=' -r '+str(r)
    jobcommand+=' -nr '+str(nr)
    jobcommand+=' -nl '+str(nl)
    jobcommand+=' -d '+str(d)
    if refineleft:
        jobcommand+=' -refineleft '+str(refineleft)
    jobcommand+=' -s '+str(sigma)
    if sas :
        jobcommand+=' -issas '+str(sas)
    
    if c0init=='gaussian':
         jobcommand+=' -xi '+str(lpvs/2)
    
    return jobcommand


def base_PVSBraincommandline(lpvs=200e-4, mesh='regular', Emem=10e4, E=10e4, K=1e-13, Kmem=1e-15, hmem=1e-4, tend=20, toutput=1e-1, toutputcycle=1e-2, dt=5e-3, r=-1, nrpvs=8, nrmem=8,nl=400, d=2e-7,rbrain=100e-4) :
    """
    Create a base for the command line

    Parameters
    ----------
    lpvs : FLOAT, optional
        length of the vessel. The default is 200e-4.
    c0init : STRING,optional
        kind of initialisation for the concentration . 
        The default is 'null' which means the initial concentration is zero everywhere in the PVS. 
        If set to 'constant' the concentration value in the PVS will be c0value.
    c0value : FLOAT, optional
        Initial value for the concentration in the SAS and in the PVS depending on c0init. The default is 50.
    sasbc : STRING, optional
        Kind of boundary condition in the SAS. 
        The default is 'scenarioB' : mass conservation, no outflow of CSF out of the SAS.
        'scenarioA' : zero concentration in the SAS - the tracer is flushed out.
        'scenarioB': massconservation no outflow
        'scenarioC' : constant CSF outflow.
        'scenarioD' : pressure dependent outflow.
    tend : FLOAT, optional
        end time of the simulation. The default is 800.
    toutput : FLOAT, optional
        period for the outputs. The default is 1.
    dt : FLOAT, optional
        time step. The default is 1e-3.
    r : FLOAT, optional
        Flow resistance at the right boundary. 
        The default is -1 which leads to a zero flow condition.
    nr : INT, optional
        Number of cells in the radial direction. The default is 8.
    nl : INT, optional
        number of cells in the longitudinal direction. The default is 100.
    d : FLOAT, optional
        Modelecular diffusion coefficient. The default is 2e-7.
    refineleft : BOOL, optional
        State if the mesh must have a geometrical progression of the left side. 
        The default is True.

    Returns
    -------
    STRING : base for the command line

    """
    

    
    jobcommand='srun -n 1 python3 PVSBrainMembrane_nitsche.py'

    jobcommand+=' -lpvs '+str(lpvs)
    jobcommand+=' -tend '+str(tend)
    jobcommand+=' -toutput '+str(toutput)
    jobcommand+=' -toutputcycle '+str(toutputcycle)
    jobcommand+=' -dt '+str(dt)
    jobcommand+=' -rbrain '+str(rbrain)
    jobcommand+=' -nrpvs '+str(nrpvs)
    jobcommand+=' -nrmem '+str(nrmem)
    jobcommand+=' -mesh '+str(mesh)
    jobcommand+=' -nl '+str(nl)
    jobcommand+=' -d '+str(d)
    jobcommand+=' -perm '+str(K)
    jobcommand+=' -permmem '+str(Kmem)
    jobcommand+=' -E '+str(E)
    jobcommand+=' -Emem '+str(Emem)
    jobcommand+=' -hmem '+str(hmem)
    

    return jobcommand



def write_state_slurm(jobname, fi, ai, rv, h0, outputfolder, file, **kargs) :
    """
    Function to write a slurm file with the proper command line to launch a simulation of a given oscillatory state.


    Parameters
    ----------
    jobname : STRING
        name of the job : will be used as label for the outputfiles.
    fi: LIST of FLOAT
        provides the frequencies (Hz) of the oscillations to be superposed
    ai: LIST of FLOAT
        provides the amplitudes (ratio) of the oscillations to be superposed
    Rv: FLOAT
        vessel mean radius (cm)
    h0: FLOAT
        PVS mean thickness (cm)
    file : STRING
        name of the slurm file to create.
    

    Returns
    -------
    None. The function creates a file.

    """
    
    
    jobcommand=base_commandline(**kargs)
    jobcommand+=' -j '+jobname
    jobcommand+=' -ai '
    
    for a in ai :
        jobcommand+=str(a)+' ' 
        
    jobcommand+=' -fi '
    for f in fi :
        jobcommand+=str(f)+' ' 
        
    jobcommand+=' -rv '+str(rv)
    jobcommand+=' -rpvs '+str(rv+h0)
    
    jobcommand+=' -o '+outputfolder
    
    slurm_template=get_slurmtemplate(jobname)
    
    with open(file, 'w') as f:
        for line in slurm_template :
            f.write(line)
        f.write('\n')
        f.write(jobcommand)
        
    f.close()
    
