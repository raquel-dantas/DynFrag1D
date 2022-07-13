from operator import index
from textwrap import indent
import akantu as aka
import numpy as np
import DFMesh2d



# Apply BCs

class FixedVelocity (aka.DirichletFunctor):
    """Fixed velocity at the boundaries."""

    def __init__(self, axis, vel):
        super().__init__(axis)
        self.axis = axis
        self.time = 0
        self.vel = vel
    
    def set_time(self, t):
        self.time = t
    
    def __call__(self, node, flags, disp, coord):
        flags[int(self.axis)] = True
        disp[int(self.axis)] = self.vel*self.time


