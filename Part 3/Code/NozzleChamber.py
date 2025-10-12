# BRUKER IKKE KODEMAL!!!!
# Skrevet av Bendik Thune

# Imports
import  numpy        as     np
import  scipy.stats  as     st
import  math         as     mt
# AST imports
import ast2000tools.constants as const

from FuelChamber import FuelChamber

class NozzleChamber(FuelChamber):
    def __init__(self,Length,Temp,NumParticles, length_nozzle): #Nozzle is square with sides equal to length
        super().__init__(Length,Temp,NumParticles,10**(-12))  #Using Class from fuelChamber
        
        self.length_nozzle = length_nozzle 
        self.Force = []
    
    def NozzleStep(self): #Finding particles that leave the nozzle
        
        Nozzle_length = self.length_nozzle/2 
        
        #Making new arrays with each dim
        x_cords = abs(self.Positions[:,0])
        y_cords = abs(self.Positions[:,1])
        z_cords = (self.Positions[:,2])    

        #Finding indexes where each dimmensjon fills their conditions
        x_indexes = np.where(x_cords < Nozzle_length)
        y_indexes = np.where(y_cords < Nozzle_length)
        z_indexes = np.where(z_cords < -self.Length/2)

        
        all_indexes = np.intersect1d(np.intersect1d(x_indexes,y_indexes),z_indexes) #finds the common indexes
        Leaving_part_vel = self.Velocities[all_indexes] #saves their velocities for later use

        #Make new particles thaat enter from the top of the gas tank
        self.Velocities[all_indexes] = np.random.normal(loc = 0, scale = self.sigma,size = (len(all_indexes),3)) 
        self.Positions[all_indexes] = (0,0, self.Length/2 *0.95)

        # Calculating and returning leaving momentum (Derived from velocity in negative z direction) and n leaving particles
        Momentum = sum(abs(Leaving_part_vel[:,2] * self.ParticleMass))
        LeavingParticles = len(Leaving_part_vel[:,2])
        return (Momentum,LeavingParticles)

    def TimeStep(self):
        self.EulerStep()
        Momentum,LeavingParticles = self.NozzleStep()
        self.CollisionStep()

        return (Momentum,LeavingParticles)

    def TimeLoop(self):

        while self.t < self.t_max:

            self.EulerStep()
            self.NozzleStep()
            self.CollisionStep()

            self.t += self.dt