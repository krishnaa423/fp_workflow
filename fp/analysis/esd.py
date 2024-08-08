#region: Modules.
import warnings
from typing import IO, Optional, Union

import numpy as np
from numpy.linalg import eigh

from ase import Atoms
from ase.optimize.optimize import UnitCellFilter

import sys
import pickle
import time
from math import sqrt
from os.path import isfile

import numpy as np

from ase.parallel import rank, barrier
from ase.io.trajectory import PickleTrajectory
from ase.io import read, write 
import os 

from fp.flows import *
from fp.io import *
from fp.analysis.esf import *
#endregion

#region: Variables.
#endregion

#region: Functions.
#endregion

#region: Classes.
class Dynamics:
    """Base-class for all MD and structure optimization classes.

    Dynamics(atoms, logfile)

    atoms: Atoms object
        The Atoms object to operate on
    logfile: file object or str
        If *logfile* is a string, a file with that name will be opened.
        Use '-' for stdout.
    trajectory: Trajectory object or str
        Attach trajectory object.  If *trajectory* is a string a
        PickleTrajectory will be constructed.  Use *None* for no
        trajectory.
    """
    def __init__(self, atoms, logfile, trajectory):
        self.atoms = atoms

        if rank != 0:
            logfile = None
        elif isinstance(logfile, str):
            if logfile == '-':
                logfile = sys.stdout
            else:
                logfile = open(logfile, 'a')
        self.logfile = logfile
        
        self.observers = []
        self.nsteps = 0

        if trajectory is not None:
            if isinstance(trajectory, str):
                trajectory = PickleTrajectory(trajectory, 'w', atoms)
            self.attach(trajectory)

    def get_number_of_steps(self):
        return self.nsteps

    def insert_observer(self, function, position=0, interval=1, 
                        *args, **kwargs):
        """Insert an observer."""
        if not callable(function):
            function = function.write
        self.observers.insert(position, (function, interval, args, kwargs))

    def attach(self, function, interval=1, *args, **kwargs):
        """Attach callback function.

        At every *interval* steps, call *function* with arguments
        *args* and keyword arguments *kwargs*."""

        if not hasattr(function, '__call__'):
            function = function.write
        self.observers.append((function, interval, args, kwargs))

    def call_observers(self):
        for function, interval, args, kwargs in self.observers:
            if self.nsteps % interval == 0:
                function(*args, **kwargs)

class Optimizer(Dynamics):
    """Base-class for all structure optimization classes."""
    def __init__(self, atoms, restart, logfile, trajectory):
        """Structure optimizer object.

        atoms: Atoms object
            The Atoms object to relax.
        restart: str
            Filename for restart file.  Default value is *None*.
        logfile: file object or str
            If *logfile* is a string, a file with that name will be opened.
            Use '-' for stdout.
        trajectory: Trajectory object or str
            Attach trajectory object.  If *trajectory* is a string a
            PickleTrajectory will be constructed.  Use *None* for no
            trajectory.
        """
        Dynamics.__init__(self, atoms, logfile, trajectory)
        self.restart = restart

        if restart is None or not isfile(restart):
            self.initialize()
        else:
            self.read()
            barrier()
    def initialize(self):
        pass

    def run(self, fmax=0.05, steps=100000000):
        """Run structure optimization algorithm.

        This method will return when the forces on all individual
        atoms are less than *fmax* or when the number of steps exceeds
        *steps*."""

        self.fmax = fmax
        step = 0
        while step < steps:
            f = self.atoms.get_forces()
            self.log(f)
            self.call_observers()
            if self.converged(f):
                return
            self.step(f)
            self.nsteps += 1
            step += 1

    def converged(self, forces=None):
        """Did the optimization converge?"""
        if forces is None:
            forces = self.atoms.get_forces()
        if hasattr(self.atoms, 'get_curvature'):
            return (forces**2).sum(axis=1).max() < self.fmax**2 and \
                   self.atoms.get_curvature() < 0.0
        return (forces**2).sum(axis=1).max() < self.fmax**2

    def log(self, forces):
        fmax = sqrt((forces**2).sum(axis=1).max())
        e = self.atoms.get_potential_energy()
        T = time.localtime()
        if self.logfile is not None:
            name = self.__class__.__name__
            self.logfile.write('%s: %3d  %02d:%02d:%02d %15.6f %12.4f\n' %
                               (name, self.nsteps, T[3], T[4], T[5], e, fmax))
            self.logfile.flush()
        
    def dump(self, data):
        if rank == 0 and self.restart is not None:
            pickle.dump(data, open(self.restart, 'wb'), protocol=2)

    def load(self):
        return pickle.load(open(self.restart))

class NDPoly:
    def __init__(self, ndims=1, order=3):
        """Multivariate polynomium.

        ndims: int
            Number of dimensions.
        order: int
            Order of polynomium."""
        
        if ndims == 0:
            exponents = [()]
        else:
            exponents = []
            for i in range(order + 1):
                E = NDPoly(ndims - 1, order - i).exponents
                exponents += [(i,) + tuple(e) for e in E]
        self.exponents = np.array(exponents)
        self.c = None
        
    def __call__(self, *x):
        """Evaluate polynomial at x."""
        return np.dot(self.c, (x**self.exponents).prod(1))

    def fit(self, x, y):
        """Fit polynomium at points in x to values in y."""
        A = (x**self.exponents[:, np.newaxis]).prod(2)
        self.c = np.linalg.solve(np.inner(A, A), np.dot(A, y))

def polyfit(x, y, order=3):
    """Fit polynomium at points in x to values in y.

    With D dimensions and N points, x must have shape (N, D) and y
    must have length N."""
    
    p = NDPoly(len(x[0]), order)
    p.fit(x, y)
    return p

class BFGS(Optimizer):
    # default parameters
    defaults = {**Optimizer.defaults, 'alpha': 70.0}

    def __init__(
        self,
        atoms: Atoms,
        restart: Optional[str] = None,
        logfile: Optional[Union[IO, str]] = '-',
        trajectory: Optional[str] = None,
        append_trajectory: bool = False,
        maxstep: Optional[float] = None,
        alpha: Optional[float] = None,
        **kwargs,
    ):
        """BFGS optimizer.

        Parameters
        ----------
        atoms: :class:`~ase.Atoms`
            The Atoms object to relax.

        restart: str
            JSON file used to store hessian matrix. If set, file with
            such a name will be searched and hessian matrix stored will
            be used, if the file exists.

        trajectory: str
            Trajectory file used to store optimisation path.

        logfile: file object or str
            If *logfile* is a string, a file with that name will be opened.
            Use '-' for stdout.

        maxstep: float
            Used to set the maximum distance an atom can move per
            iteration (default value is 0.2 Å).

        alpha: float
            Initial guess for the Hessian (curvature of energy surface). A
            conservative value of 70.0 is the default, but number of needed
            steps to converge might be less if a lower value is used. However,
            a lower value also means risk of instability.

        kwargs : dict, optional
            Extra arguments passed to
            :class:`~ase.optimize.optimize.Optimizer`.

        """
        if maxstep is None:
            self.maxstep = self.defaults['maxstep']
        else:
            self.maxstep = maxstep

        if self.maxstep > 1.0:
            warnings.warn('You are using a *very* large value for '
                          'the maximum step size: %.1f Å' % self.maxstep)

        self.alpha = alpha
        if self.alpha is None:
            self.alpha = self.defaults['alpha']
        Optimizer.__init__(self, atoms=atoms, restart=restart,
                           logfile=logfile, trajectory=trajectory,
                           append_trajectory=append_trajectory,
                           **kwargs)

    def initialize(self):
        # initial hessian
        self.H0 = np.eye(3 * len(self.optimizable)) * self.alpha

        self.H = None
        self.pos0 = None
        self.forces0 = None

    def read(self):
        file = self.load()
        if len(file) == 5:
            (self.H, self.pos0, self.forces0, self.maxstep,
             self.atoms.orig_cell) = file
        else:
            self.H, self.pos0, self.forces0, self.maxstep = file

    def step(self, forces=None):
        optimizable = self.optimizable

        if forces is None:
            forces = optimizable.get_forces()

        pos = optimizable.get_positions()
        dpos, steplengths = self.prepare_step(pos, forces)
        dpos = self.determine_step(dpos, steplengths)
        optimizable.set_positions(pos + dpos)
        if isinstance(self.atoms, UnitCellFilter):
            self.dump((self.H, self.pos0, self.forces0, self.maxstep,
                       self.atoms.orig_cell))
        else:
            self.dump((self.H, self.pos0, self.forces0, self.maxstep))

    def prepare_step(self, pos, forces):
        forces = forces.reshape(-1)
        self.update(pos.flat, forces, self.pos0, self.forces0)
        omega, V = eigh(self.H)

        # FUTURE: Log this properly
        # # check for negative eigenvalues of the hessian
        # if any(omega < 0):
        #     n_negative = len(omega[omega < 0])
        #     msg = '\n** BFGS Hessian has {} negative eigenvalues.'.format(
        #         n_negative
        #     )
        #     print(msg, flush=True)
        #     if self.logfile is not None:
        #         self.logfile.write(msg)
        #         self.logfile.flush()

        dpos = np.dot(V, np.dot(forces, V) / np.fabs(omega)).reshape((-1, 3))
        steplengths = (dpos**2).sum(1)**0.5
        self.pos0 = pos.flat.copy()
        self.forces0 = forces.copy()
        return dpos, steplengths

    def determine_step(self, dpos, steplengths):
        """Determine step to take according to maxstep

        Normalize all steps as the largest step. This way
        we still move along the direction.
        """
        maxsteplength = np.max(steplengths)
        if maxsteplength >= self.maxstep:
            scale = self.maxstep / maxsteplength
            # FUTURE: Log this properly
            # msg = '\n** scale step by {:.3f} to be shorter than {}'.format(
            #     scale, self.maxstep
            # )
            # print(msg, flush=True)

            dpos *= scale
        return dpos

    def update(self, pos, forces, pos0, forces0):
        if self.H is None:
            self.H = self.H0
            return
        dpos = pos - pos0

        if np.abs(dpos).max() < 1e-7:
            # Same configuration again (maybe a restart):
            return

        dforces = forces - forces0
        a = np.dot(dpos, dforces)
        dg = np.dot(self.H, dpos)
        b = np.dot(dpos, dg)
        self.H -= np.outer(dforces, dforces) / a + np.outer(dg, dg) / b

    def replay_trajectory(self, traj):
        """Initialize hessian from old trajectory."""
        if isinstance(traj, str):
            from ase.io.trajectory import Trajectory
            traj = Trajectory(traj, 'r')
        self.H = None
        atoms = traj[0]
        pos0 = atoms.get_positions().ravel()
        forces0 = atoms.get_forces().ravel()
        for atoms in traj:
            pos = atoms.get_positions().ravel()
            forces = atoms.get_forces().ravel()
            self.update(pos, forces, pos0, forces0)
            pos0 = pos
            forces0 = forces

        self.pos0 = pos0
        self.forces0 = forces0

class Esd(BFGS):
    def __init__(self, fmax=0.05, max_steps=10, *args, **kwargs):
        self.fmax = fmax 
        self.max_steps = max_steps 
        self.pos_new = None 
        self.is_converged = False
        self.is_max_steps = False 
        super().__init__(*args, **kwargs)
        self.nsteps = -1 

    def dump_hdf5(self, data):
        H, pos0, pos_new, forces0, nsteps, maxstep, is_converged, is_max_steps = data 

        with h5py.File('esd.h5', 'w') as w:
            ds_H = w.create_dataset('H', shape=(maxstep, H.shape[0], H.shape[1]), dtype=H.dtype); ds_H[:] = np.zeros((maxstep, H.shape[0], H.shape[1])); ds_H[nsteps, ...] = H 
            ds_pos0 = w.create_dataset('pos0', shape=(maxstep, pos0.shape[0]), dtype=pos0.dtype); ds_pos0[:] = np.zeros((maxstep, pos0.shape[0])); ds_pos0[nsteps, ...] = pos0 
            ds_pos_new = w.create_dataset('pos_new', shape=(maxstep, pos_new.shape[0]), dtype=pos_new.dtype); ds_pos_new[:] = np.zeros((maxstep, pos_new.shape[0])); ds_pos_new[nsteps, ...] = pos_new 
            ds_forces0 = w.create_dataset('forces0', shape=(maxstep, forces0.shape[0]), dtype=forces0.dtype); ds_forces0[:] = np.zeros((maxstep, forces0.shape[0])); ds_forces0[nsteps, ...] = forces0 
            w.create_dataset('nsteps', data=nsteps)
            w.create_dataset('maxstep', data=maxstep)
            w.create_dataset('is_converged', data=is_converged)
            w.create_dataset('is_max_steps', data=is_max_steps)

    def load_hdf5(self, step_idx=-1):
        with h5py.File('esd.h5', 'r') as r:
            nsteps = r['nsteps'][()] if step_idx==-1 else step_idx 
            H = r['H'][nsteps, ...]
            pos0 = r['pos0'][nsteps, ...]
            pos_new = r['pos_new'][nsteps, ...]
            forces0 = r['forces0'][nsteps, ...]
            maxstep = r['maxstep'][()]
            is_converged = r['is_converged'][()]
            is_max_steps = r['is_max_steps'][()]

        return (H, pos0, pos_new, forces0, nsteps, maxstep, is_converged, is_max_steps)

    def read_hdf5(self, step_idx=-1):
        file = self.load_hdf5(step_idx)
        self.H, self.pos0, self.pos_new, self.forces0, self.nsteps, self.maxstep, self.is_converged, self.is_max_steps = file

    # def dump(self, data):
    #     # # TODO: 
    #     # if rank == 0 and self.restart is not None:
    #     #     pickle.dump(data, open(self.restart, 'wb'), protocol=2)
    #     self.dump_hdf5(data)

    # def load(self, step_idx=-1):
    #     # # TODO: 
    #     # return pickle.load(open(self.restart))
    #     return self.load_hdf5(step_idx)

    # def read(self, step_idx=-1):
    #     # # TODO: 
    #     # file = self.load()  
    #     # self.H, self.pos0, self.pos_new, self.forces0, self.nsteps, self.maxstep, self.is_converged, self.is_max_steps = file
    #     self.read_hdf5(step_idx)

    def converged(self, forces=None):
        """Did the optimization converge?"""
        # if forces is None:
        #     forces = self.atoms.get_forces()
        if hasattr(self.atoms, 'get_curvature'):
            return (forces**2).sum(axis=1).max() < self.fmax**2 and \
                   self.atoms.get_curvature() < 0.0
        return (forces**2).sum(axis=1).max() < self.fmax**2

    @staticmethod
    def static_step(fmax, max_steps, ymlfile):
        # Read atoms.
        atoms = read('./esd_atoms.xsf') if os.path.exists('./esd_atoms.xsf') else read('./sc_atoms.xsf')
        atoms.calc = EsfCalculator()
        
        # Read esd if present.
        esd = Esd(fmax=fmax, max_steps=max_steps, atoms=atoms)
        if os.path.exists('./esd.h5'): esd.read_hdf5()
        if esd.is_converged:
            print('ESD has converged')
            return 
        if esd.is_max_steps:
            print(f'ESD has reached max steps: {esd.max_steps}')
            return
        
        # Create the flow.
        fullgridflow = FullGridFlow.from_yml(ymlfile)
        flowmanage = fullgridflow.get_flowmanage([
                # Relax,
                Scf,
                Dfpt,
                # Phbands,
                # Phdos,
                # Phmodes,
                # Dos,
                # Pdos,
                Dftelbands,
                # Kpdos,
                # Wannier,
                Wfn,
                Epw,
                Wfnq,
                # Wfnfi,
                # Wfnqfi,
                Epsilon,
                Sigma,
                Inteqp,
                Kernel,
                Absorption,
                # Plotxct,
                # Bseq,
                Esf,
                # Esd(input=input),
                # XctPh(input=input),
                # Pol(input=input),
                # Xctpol(input=input),
            ], 
            save_pkl=True
        )
        flowmanage.create()
        flowmanage.create_job_all_script(
            filename='job_esd_step_flow.sh',
            start_job='job_scf.sh',
            stop_job='job_esf.sh',
            flowfile_to_read='flowmanage.pkl',
        )

        total_time = run_and_wait_command('./job_esd_step_flow.sh', total_time=0.0) 
        print(f'Total time for step flow: {esd.nsteps} is: {total_time}')
        esd.step()

    def step(self, forces=None):

        atoms = self.atoms

        if forces is None:
            forces = atoms.get_forces()

        # Return if converged.
        if self.converged(forces):
            self.is_converged = True 
            self.dump_hdf5((self.H, self.pos0, self.pos_new, self.forces0, self.nsteps, self.maxstep, self.is_converged, self.is_max_steps))
            return 
        
        # Return if is max_steps.
        if not self.nsteps < self.max_steps:
            self.is_max_steps = True 
            self.dump_hdf5((self.H, self.pos0, self.pos_new, self.forces0, self.nsteps, self.maxstep, self.is_converged, self.is_max_steps))
            return

        pos = atoms.get_positions()
        dpos, steplengths = self.prepare_step(pos, forces)
        dpos = self.determine_step(dpos, steplengths)
        self.pos_new = pos + dpos
        
        # Update new position. 
        atoms.set_positions(self.pos_new)
        # Increment the number of steps.
        self.nsteps += 1
        
        # Save and move all to save folder.
        self.dump_hdf5((self.H, self.pos0, self.pos_new, self.forces0, self.nsteps, self.maxstep, self.is_converged, self.is_max_steps))
        save_folder = f'./esd/iter_{self.nsteps}'
        os.system('mkdir -p ./esd')
        os.system(f'mkdir -p {save_folder}')
        flow: FlowManage = load_obj('./flowmanage.pkl')
        flow.save_job_results(save_folder)

        # Write the new positions for the next step.
        write('./esd_atoms.xsf', atoms)
#endregion
