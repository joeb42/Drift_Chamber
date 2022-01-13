import numpy as np
from scipy import stats, sparse
import math
import copy
import matplotlib

matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import LogNorm
import matplotlib.cm
from muon import Muon, MissedDetector


class DriftChamber:
    """
    Drift Chamber object

    Attributes:
    grid: numpy 2d array representing charge distribution in drift chamber
    spacing: grid spacing (floating point)
    _mat: (private) sparse matrix for iterative computation of drift diffusion of charge
    lu: sparse lu factorisation of mat, implemented with a getter @property decorator to ensure that any changes
    in the electric field or diffusivity are reflected in the matrix
    fig: Matplotlib figure for animating drift diffusion
    ax: axes for animation
    im: image representing grid of charge for animation
    anim: Matplotlib FuncAnimation instance, needs to be stored in class scope
    so not garbage collected in order for animation to work
    E: numerical electric field
    D: numerical diffusivity
    """

    def __init__(self, spacing=1.0, E=1e5, D=0.1, tau=1e-6, mat=True):
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
            self._mat = sparse.diags(
                [
                    1 + 4 * self.alpha + self.beta,
                    -self.alpha,
                    -self.alpha - self.beta,
                    -self.alpha,
                    -self.alpha,
                ],
                [0, 1, -1, ny, -ny],
                shape=(nz * ny, nz * ny),
                format="csc",
            )
            # Lu factorise for faster solving
            self._lu = sparse.linalg.spilu(self._mat)

    def __str__(self):
        """String representation of DriftChamber object for ease of debugging"""
        return f"""Drift Chamber object, grid spacing of {self.spacing}cm, 
        diffusivity of {self.D}, E-field of {self.E}."""

    def initial_ion(self, muon):
        """Updates charge grid to reflect initial ionisation pattern from incident muon.
        Takes muon (Muon object) as a required positional arg.
        """
        #  Get starting coordinates in units of spacing
        y, z = (
            muon.starting_coords()[0] / self.spacing,
            muon.starting_coords()[1] / self.spacing,
        )
        # Boolean variable right, True if muon moving left to right
        right = np.sin(muon.zen) * np.sin(muon.azi) > 0
        while True:
            #  Get values of next intersections in y and z directions
            if right:
                y_next = math.floor(y + 1)
            else:
                y_next = math.ceil(y - 1)
            z_next = math.ceil(z - 1)
            t = min(
                (z_next - z) / np.cos(muon.zen),
                (y_next - y) / (np.sin(muon.zen) * np.sin(muon.azi)),
            )
            # Update coordinates
            y += t * np.sin(muon.zen) * np.sin(muon.azi)
            z += t * np.cos(muon.zen)
            #  Check if new coordinates outside detector
            if z < 0 or y < 0 or y > 50 / self.spacing:
                return
            # Get number of ionised electrons, assuming all energy loss is through ionisation
            # Thus energy loss is quantised and so rounded to the nearest smallest integer by casting to int
            d = t * self.spacing
            n_electrons = int(stats.moyal.rvs(94 * d, (94 * d) ** 0.5))
            #  Convert coordinates to grid indices
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
        Using the property decorator allows flexibility to easily change the matrix
        mid simulation (e.g by turning off the electric field) by simply changing the values of the alpha or beta
        attribute.  First checks if alpha or beta attributes have changed, then reconstructs sparse matrix
        and lu factorises if either have changed.
        """
        if (
            self._mat[1, 0] == -self.alpha - self.beta
            and self._mat[0, 1] == -self.alpha
        ):
            return self._lu
        nz, ny = self.grid.shape
        self._mat = sparse.diags(
            [
                1 + 4 * self.alpha + self.beta,
                -self.alpha,
                -self.alpha - self.beta,
                -self.alpha,
                -self.alpha,
            ],
            [0, 1, -1, ny, -ny],
            shape=(nz * ny, nz * ny),
            format="csc",
        )
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

    def __init(self):
        """
        Private method only used for animation.
        """
        self.im.set_data(np.zeros(self.grid.shape))
        return [self.im]

    def __ani(self, i):
        """
        Private method updates animation at each frame. Calls drift_diff method at each frame.
        """
        self.drift_diff()
        self.im.set_array(self.grid)
        self.fig.suptitle(f"t = {round(i * self.tau, 6)} seconds")
        return [self.im]

    def animate(self, vmin=1e-4, vmax=1):
        """
        Plots a Matplotlib funcanimation of charge drift-diffusion.
        Takes positive numerical vmin and vmax arg as min, max values for log scale colour normalisation.
        """
        self.fig = plt.figure()
        self.ax = plt.axes(xlim=(0, 50), ylim=(0, 30))
        # Copy the jet cmap from matplotlib
        my_cmap = copy.copy(matplotlib.cm.get_cmap("jet"))
        # Set bad (zero value) pixels to dark blue
        my_cmap.set_bad((0, 0, 0.5))
        self.im = plt.imshow(
            self.grid,
            origin="lower",
            cmap=my_cmap,
            extent=[0, 50, 0, 30],
            norm=LogNorm(vmin=vmin, vmax=vmax),
        )
        self.ax.set_xlabel("y coordinate (cm)")
        self.ax.set_ylabel("z coordinate (cm)")
        self.anim = animation.FuncAnimation(
            self.fig,
            self.__ani,
            init_func=self.__init,
            blit=True,
            interval=1,
            repeat=False,
        )
        plt.colorbar(label=r"Charge density ($em^{-2}$)")
        plt.show()
