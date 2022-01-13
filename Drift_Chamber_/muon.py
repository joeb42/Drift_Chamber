import numpy as np


class Muon:
    """Cosmic ray muon object. Coordinate positions are taken relative to an origin in bottom left corner of detector.

    Attributes:
    azim: Azimuthal angle, random value distributed evenly between 0 and 2pi radians
    zen: Zenith angle, random value distributed proportional to cos^2 between 3pi/2 and pi/2
    energy: Random log-normal value in GeV of mean 6.55 and stdev 1.8
    charge: Randomly generated charge of muon takes values +1 (with probability 0.53) and -1 (with probability 0.47)
    x_coord: x coordinate of muon on plane in cm (direction into the page in relation to diagram above)
    y_coord: y coordinate of muon on plane in cm (horizontal direction in relation to diagram above)
    height: height of plane above detector in cm (vertical direction in relation to diagram above). Must be >= 30 (height of detector)
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
            raise Exception(
                f"{height} arg is too small! Must be > 30 as muons must start above detector."
            )
        self.height = height
        self.x_coord = width * np.random.random() - 0.5 * width
        self.y_coord = length * np.random.random() - 0.5 * length + 25

    def __str__(self):
        """
        String representation of Muon object for ease of debugging and for use in MissedDetector exception
        """
        return f"""Cosmic ray muon of energy {round(self.energy)} GeV, with charge {self.charge} e, 
    zenith angle {round(self.zen, 2)} radians and azimuthal angle {round(self.azi, 2)} radians, 
y coord {round(self.y_coord, 2)}"""

    def starting_coords(self):
        """
        Return coordinates of where muon enters detector.
        Raises MissedDetector if muon path never intersects detector.
        """
        is_edge = (
            lambda t: 0 <= self.height + t * np.cos(self.zen) <= 30
            and 0 <= self.y_coord + t * np.sin(self.zen) * np.sin(self.azi) <= 50
        )
        intersections = list(
            filter(
                is_edge,
                (
                    (30 - self.height) / np.cos(self.zen),
                    (-self.y_coord) / (np.sin(self.zen) * np.sin(self.azi)),
                    (50 - self.y_coord) / (np.sin(self.azi) * np.sin(self.zen)),
                ),
            )
        )
        if not intersections:
            raise MissedDetector(self)
        t0 = min(intersections)
        return self.y_coord + t0 * np.sin(self.zen) * np.sin(
            self.azi
        ), self.height + t0 * np.cos(self.zen)


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
        self.message = "Muon has missed detector!"
        super().__init__(self.message)

    def __str__(self):
        """
        Print string repr of muon object to get some info as to why exception has occured from muon attributes.
        """
        return f"{self.muon} -> {self.message}"
