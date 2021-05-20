import numpy as np
import tkinter as tk
from tkinter import ttk
import math
from scipy import stats, optimize, sparse
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg)
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import LogNorm
import matplotlib.cm
import copy
from matplotlib.figure import Figure


class Muon:
    """Cosmic ray muon object. Coordinate positions are taken relative to an origin in bottom left corner of detector.
    Attributes:
    azim: Azimuthal angle, random value distributed evenly between 0 and 2pi radians
    zen: Zenith angle, random value distributed proportional to cos^2 between 3pi/2 and pi/2
    energy: Random log-normal value in GeV of mean 6.55 and stdev 1.8
    charge: Randomly generated charge of muon takes values +1 (with probability 0.53) and -1 (with probability 0.47)
    x_coord: x coordinate of muon on plane in cm (direction into the page in relation to diagram)
    y_coord: y coordinate of muon on plane in cm (horizontal direction in relation to diagram)
    height: height of plane above detector in cm (vertical direction in relation to diagram). Must be >= 30 (height of detector)
    """

    def __init__(self, width=100, length=80, height=60):
        """
        Initialises random muon attributes and takes horizontal dimensions and height of plane above detector
        as optional kwargs
        """
        self.energy = np.random.lognormal(6.55, 1.8)
        if np.random.random() > 0.43:
            self.charge = 1
        else:
            self.charge = -1
        # Accept/reject method to get random zenith angle, distributed proportional to cos squared
        x, y = np.pi * np.random.random() + 0.5 * np.pi, np.random.random()
        while y >= (np.cos(x)) ** 2:
            x, y = np.pi * np.random.random() + 0.5 * np.pi, np.random.random()
        self.zen = x
        self.azi = 2 * np.pi * np.random.random()
        if height < 30:
            raise Exception(f'{height} arg is too small! Must be > 30 as muons must start above detector.')
        self.height = height
        self.x_coord = width * np.random.random() - 0.5 * width
        self.y_coord = length * np.random.random() - 0.5 * length + 25

    def __str__(self):
        return f'''Cosmic ray muon of energy {round(self.energy)} GeV, with charge {self.charge} e, 
    zenith angle {round(self.zen, 2)} radians and azimuthal angle {round(self.azi, 2)} radians, 
y coord {round(self.y_coord, 2)}'''

    def starting_coords(self):
        """
        Return coordinates of where muon enters detector.
        Raises MissedDetector if muon path never intersects detector.
        """
        is_edge = lambda t: 0 <= self.height + t * np.cos(self.zen) <= 30 and 0 <= self.y_coord + t * np.sin(self.zen) * np.sin(self.azi) <= 50
        intersections = list(filter(is_edge, ((30 - self.height) / np.cos(self.zen),
                                              (-self.y_coord) / (np.sin(self.zen) * np.sin(self.azi)), (50 - self.y_coord) /
                                              (np.sin(self.azi) * np.sin(self.zen)))))
        if not intersections:
            raise MissedDetector(self)
        t0 = min(intersections)
        return self.y_coord + t0 * np.sin(self.zen) * np.sin(self.azi), self.height + t0 * np.cos(self.zen)


class MissedDetector(Exception):
    """
    Exception raised when Muon object attributes result in the muon missing the chamber entirely.
    """

    def __init__(self, muon):
        """
        Takes Muon object that raises exception as positional arg.
        Constructor stores muon object in class scope as well as message string.
        """
        self.muon = muon
        self.message = 'Muon has missed detector!'
        super().__init__(self.message)

    def __str__(self):
        """
        String repr of muon object to get some info as to why exception has occured from muon attributes.
        """
        return f'{self.muon} -> {self.message}'


class DriftChamber:
    """
    Drift Chamber object
    Attributes:
    grid: numpy 2d array representing charge distribution in drift chamber
    spacing: grid spacing (floating point)
    lu: sparse lu factorisation of mat, implemented with a getter @property decorator to ensure that any changes
    in the electric field or diffusivity are reflected in the matrix
    E: numerical electric field
    D: numerical diffusivity
    tau: positive numerical timestep for drift diffusion matrix
    """

    def __init__(self, spacing=1., E=1e5, D=0.1, tau=1e-6, mat=True):
        """
        Contructor initialises class attributes as well as charge grid and drift diffusion matrix (also lu factorises).
        Takes positive numerical spacing, numerical E (electric field), positive numerical D (diffusivity)
        and positive numerical tau (timestep) as optional kwargs.
        Also takes mat kwarg, True by default, this should be set to False if need to avoid constructing matrix
        (e.g for investigating small spacings in the initial ionisation track that would otherwise generate very large matricies).
        """
        self.grid = np.zeros((int(30 / spacing), (int(50 / spacing))))
        self.spacing = spacing
        self.tau = tau
        self.E = E
        self.D = D
        if mat:
            # Construct matrix
            self.alpha = D * tau / ((spacing * 1e-2) ** 2)
            self.beta = tau * E * 0.02 / (2 * spacing * 1e-2)
            nz, ny = self.grid.shape
            self._mat = sparse.diags([1 + 4 * self.alpha + self.beta, -self.alpha, -self.alpha - self.beta,
                                      -self.alpha, -self.alpha], [0, 1, -1, ny, -ny], shape=(nz * ny, nz * ny),
                                     format='csc')
            # Lu factorise for faster solving
            self._lu = sparse.linalg.spilu(self._mat)

    def __str__(self):
        return f'''Drift Chamber object, grid spacing of {self.spacing}cm, 
        diffusivity of {self.D}, E-field of {self.E}.'''

    def initial_ion(self, muon):
        """Updates charge grid to reflect initial ionisation pattern from incident muon.
        Takes muon (Muon object) as a required positional arg.
        Here the energy loss is modelled using the moyal distribution which is an accurate description 
        of the energy loss of high energy radiation in matter.
        """
        #  Get starting coordinates in units of spacing
        y, z = muon.starting_coords()[0] / self.spacing, muon.starting_coords()[1] / self.spacing
        # Boolean variable right, True if muon moving left to right
        right = np.sin(muon.zen) * np.sin(muon.azi) > 0
        while True:
            #  Get values of next intersections in y and z directions
            if right:
                y_next = math.floor(y + 1)
            else:
                y_next = math.ceil(y - 1)
            z_next = math.ceil(z - 1)
            t = min((z_next - z) / np.cos(muon.zen), (y_next - y) / (np.sin(muon.zen) * np.sin(muon.azi)))
            # Update coordinates
            y += t * np.sin(muon.zen) * np.sin(muon.azi)
            z += t * np.cos(muon.zen)
            #  Check if new coordinates outside detector
            if z < 0 or y < 0 or y > 50 / self.spacing:
                return
            # Get number of ionised electrons, assuming all energy loss is through ionisation
            # Thus energy loss is quantised and so rounded to the nearest smallest integer by wrapping in int
            d = t * self.spacing
            n_electrons = int(stats.moyal.rvs(94 * d, (94 * d) ** 0.5))
            # Convert coordinates to grid indices
            if right:
                j = math.ceil(y - 1)
            else:
                j = math.floor(y)
            i = math.floor(z)
            # Only take into account energy losses so take negative losses as 0 ionisation
            self.grid[i, j] = max(n_electrons / (self.spacing ** 2), 0)

    @property
    def lu(self):
        """
        Lu factorisation of drift diffusion matrix. 
        @property used to allow flexibility in changing E field mid simulation
        """
        # Check if matrix corresponds to other attributes
        if self._mat[1, 0] == -self.alpha - self.beta and self._mat[0, 1] == -self.alpha:
            return self._lu
        # Recontruct matrix if there is a discrepancy 
        nz, ny = self.grid.shape
        self._mat = sparse.diags([1 + 4 * self.alpha + self.beta, -self.alpha, -self.alpha - self.beta,
                                  -self.alpha, -self.alpha], [0, 1, -1, ny, -ny], shape=(nz * ny, nz * ny),
                                 format='csc')
        self._lu = sparse.linalg.spilu(self._mat)
        return self._lu

    def drift_diff(self):
        """
        Updates charge grid by solving matrix equation for single iteration.
        Dampens periodicity by setting all grid cells within 2cm from boundary to 0.
        """
        # Flatten 2d array of charge
        q = self.grid.flatten()
        # Solve matrix equation to update charge values
        q = self.lu.solve(q)
        # Update 2d array of charge by reshaping q
        self.grid = q.reshape(self.grid.shape)
        # Set all cells two 2cm in from each boundary to zero to dampen periodicity
        boundary = int(2 / self.spacing)
        self.grid[:, -boundary:] = 0
        self.grid[:, :boundary] = 0
        self.grid[-boundary:, :] = 0
        self.grid[:boundary, :] = 0


class Widget:
    """Tkinter window allows user to select grid spacing, generate random muon tracks and start/stop animation"""

    def __init__(self):
        """Initialises window with buttons and matplotlib figure"""
        # Set up matplotlib figure
        self.fig = Figure(figsize=(6, 4), dpi=150)
        self.sp = self.fig.add_subplot(111)
        self.sp.set_xlabel('y coordinate (cm)')
        self.sp.set_ylabel('z coordinate (cm)')
        self.sp.set_ylim(0, 30)
        self.sp.set_xlim(0, 50)
        # Set up Tkinter window
        self.root = tk.Tk()
        self.root.title('Drift Chamber')
        # Integrate matplotlib figure into window
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.root)
        self.canvas.get_tk_widget().grid(column=0, row=3)
        spacings_label = tk.Label(self.root, text='Select a grid spacing')
        spacings_label.grid(row=0, column=1)
        # Combobox to get grid spacing from user
        self.spacings = ttk.Combobox(self.root, values=['1', '0.5', '0.1'], state='readonly')
        self.spacings.grid(row=1, column=1)
        spacing_button = tk.Button(self.root, text='Press to confirm spacing', command=self.spacing_press)
        spacing_button.grid(row=2, column=1)
        # Button to generate muon track
        muon_button = tk.Button(self.root, text='Press to generate random muon track', command=self.gen_muon)
        muon_button.grid(row=3, column=1)
        # Animation start button
        ani_button = tk.Button(self.root, text='Press to start animation', command=self.ani)
        ani_button.grid(row=4, column=1)
        # Animation stop button
        stop_button = tk.Button(self.root, text='Press to stop animation', command=self.stop_press)
        stop_button.grid(row=5, column=1)
        clear_button = tk.Button(self.root, text='Press to clear axis', command=self.clear_press)
        clear_button.grid(row=6, column=1)
        tk.mainloop()

    def stop_press(self):
        """Button command stops animation by setting stop variable to True"""
        self.stop = True

    def clear_press(self):
        """Button command clears plot and rebuilds blank axis"""
        self.chamber.grid *= 0.
        self.fig.clf()
        self.sp = self.fig.add_subplot(111)
        self.sp.set_xlabel('y coordinate (cm)')
        self.sp.set_ylabel('z coordinate (cm)')
        self.sp.set_ylim(0, 30)
        self.sp.set_xlim(0, 50)
        self.canvas.draw()

    def spacing_press(self):
        """Button sets grid spacing to user input from spacings combobox"""
        self.chamber = DriftChamber(spacing=float(self.spacings.get()))

    def gen_muon(self):
        """Generates a random muon track, handling MissedDetector exceptions"""
        # Need custom cmap to handle 0 value pixels in lognorm normalisation
        # Copy the jet cmap from matplotlib
        my_cmap = copy.copy(matplotlib.cm.get_cmap('jet'))
        # Set bad (zero value) pixels to dark blue
        my_cmap.set_bad((0, 0, 0.5))
        while True:
            try:
                m = Muon()
                self.chamber.initial_ion(m)
                self.im = self.sp.imshow(self.chamber.grid, cmap=my_cmap, origin='lower', extent=[0, 50, 0, 30],
                                         norm=LogNorm(vmin=1e-2, vmax=100))
                self.fig.colorbar(self.im, label=r"Charge density ($em^{-2}$)")
                self.canvas.draw()
                return
            except MissedDetector:
                continue

    def __animate(self, i):
        """Animation function called at each frame. Updates charge grid by iteratively solving matrix equation by 
        called drift_diff method from DriftChamber instance"""
        self.chamber.drift_diff()
        self.im.set_array(self.chamber.grid)
        self.fig.suptitle(f't = {round(i * self.chamber.tau, 6)} seconds')
        return [self.im]

    def __gen(self):
        """Stops animation when stop attribute set to False by stop button"""
        i = 0
        while not self.stop:
            i += 1
            yield i

    def ani(self):
        """Button commmand begins animation. Stores animation object in class scope to avoid garbage collection
        which is neccessary for animation to work"""
        self.stop = False
        self.an = animation.FuncAnimation(self.fig, self.__animate, interval=15, frames=self.__gen, repeat=False)
        self.canvas.draw()
        
        
if __name__ == '__main__':
    w = Widget()
