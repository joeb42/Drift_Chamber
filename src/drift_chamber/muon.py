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

    POSITIVE_CHARGE_PROBABILITY = 0.53

    def __init__(self, width: int = 100, length: int = 80, height: int = 60) -> None:
        """
        Initialises random muon attributes and takes horizontal dimensions and height of plane above detector
        as optional kwargs
        """
        self.energy: float = np.random.lognormal(6.55, 1.8)
        r = np.random.random()
        self.charge: int = 1 if r > (1 - Muon.POSITIVE_CHARGE_PROBABILITY) else -1
        # Accept/reject method to get random zenith angle, distributed proportional to cos squared
        x, y = np.pi * np.random.random() + 0.5 * np.pi, np.random.random()
        while y >= (np.cos(x)) ** 2:
            x, y = np.pi * np.random.random() + 0.5 * np.pi, np.random.random()
        self.zen: float = x
        self.azi: float = 2 * np.pi * np.random.random()
        self.height = height
        self.x_coord: float = width * np.random.random() - 0.5 * width
        self.y_coord: float = length * np.random.random() - 0.5 * length + 25

    def __repr__(self) -> str:
        """
        String representation of Muon object for ease of debugging and for use in MissedDetector exception
        """

        return f"Muon: energy: {self.energy:.0f} GeV, charge: {self.charge}e, zenith: {self.zen:.3f}, azimuthal: {self.azi:.3f}, x: {self.x_coord:.3f}, y: {self.y_coord:.3f}"

    def starting_coords(self) -> tuple[float, float]:
        """
        Return coordinates of where muon enters detector.
        Raises MissedDetector if muon path never intersects detector.
        """

        def is_edge(t) -> bool:
            return (0 <= self.height + t * np.cos(self.zen) <= 30) and (
                0 <= self.y_coord + t * np.sin(self.zen) * np.sin(self.azi) <= 50
            )

        edges = [
            (30 - self.height) / np.cos(self.zen),
            (-self.y_coord) / (np.sin(self.zen) * np.sin(self.azi)),
            (50 - self.y_coord) / (np.sin(self.azi) * np.sin(self.zen)),
        ]
        intersections = [edge for edge in edges if is_edge(edge)]
        if not intersections:
            raise MissedDetector(self)
        t0 = min(intersections)
        x = self.y_coord + t0 * np.sin(self.zen) * np.sin(self.azi)
        y = self.height + t0 * np.cos(self.zen)
        return x, y


class MissedDetector(Exception):
    """
    Exception raised when Muon object attributes result in the muon missing the chamber entirely.
    """

    def __init__(self, muon: Muon):
        """
        Takes Muon object that raises exception as positional arg.
        Constructor stores muon object in class scope as well as message string.
        """
        self.muon = muon
        self.message = "Muon has missed detector!"
        super().__init__(self.message)

    def __str__(self) -> str:
        """
        Print string repr of muon object to get some info as to why exception has occured from muon attributes.
        """
        return f"{self.muon} -> {self.message}"
