# BRUKER IKKE KODEMAL!!!!
# Skrevet av Bastian Eggum Huuse

# Imports
import  numpy        as     np
import  math         as     mt
import  matplotlib.pyplot as plt
import matplotlib.animation as animation
# AST imports
import ast2000tools.constants as const

class FuelChamber:

    def __init__(self,Length,Temp,NumParticles,dt = 10**(-12)):

        # Time parameters
        self.t     = 0
        self.dt    = dt
        self.t_max = 10**(-9)

        # Chamber Parameters
        self.Length       = Length
        self.Temp         = Temp
        self.NumParticles = NumParticles
        self.ParticleMass = const.m_H2

        # Pressure, which we track for stats
        self.TotalPressure = 0
        self.counter = 0

        # Points, to animate
        self.SelectPositions = np.zeros((1000,500,3))
        self.NumPoints = 500
        self.i = 0

        # Maxwell-boltzmann deviation, used for generating the particle velocities later
        self.sigma = ((self.Temp*const.k_B)/self.ParticleMass)**(1/2)
        # Method for generating a random velocity, according to the 
        self.MaxwellBoltzmann = lambda : np.random.normal(loc = 0, scale = self.sigma) 

        # Arrays we will store positions and velocities in later
        self.Positions = np.zeros((self.NumParticles,3))
        self.Velocities = np.zeros((self.NumParticles,3))

        # There isn't any better place to put this so it goes here.
        # The coordinate system inside the chamber goes from -self.Length/2 to self.Length/2.
        # This means that the origin is in the center of the box, which is nice

        self.SetPositionsAndVelocities()

    def SetPositionsAndVelocities(self):

        # We essentially want to tile a 3D grid (a box)
        # We tile within a slightly smaller box than the actual chamber, so no particles begin within the walls
        SmallerHalf = (self.Length * 0.95) / 2
        NumParticlesPerLine = mt.ceil((self.NumParticles)**(1/3)) # Since we want all the particles to fit within a cube, we find one side length of said cube
        SpacedPoints = np.linspace(-SmallerHalf,SmallerHalf,NumParticlesPerLine) # Linearly spacing points along previously mentioned side length

        # In the case where the total number of particles isn't a cube number, we assume that it is the closest cube number (rounded up)
        # We then fill the cube up until we have generated all positions, meaning that the cube goes partially unfinished.

        # Looping through the spaced points, nested twice, to create our cube of particles
        # This is not particularily fast, but we only need to do this once, so it's ok :)
        n = 0 # Iterator
        for i in SpacedPoints:
            for j in SpacedPoints:
                for k in SpacedPoints:
                    self.Positions[n] = np.array([i,j,k]) # Setting the position of the current particle
                    # In the flowchart, velocities and positions were generated in different loops, but its more efficient to put them both in the same loop :)
                    self.Velocities[n] = np.array([self.MaxwellBoltzmann(),self.MaxwellBoltzmann(),self.MaxwellBoltzmann()])
                    
                    n += 1 # We keep track of how many particles we have generated, so we can break out of the loop when finished
                    if(n == self.NumParticles): # Because of the syntax we see in the lines below, we only need to break out of this loop
                        break

                else: # This is an incredibly strange way of breaking nested loops that i found on stackoverflow
                    continue
                break # Breaking out of the second loop
            else: # It seems that python reads the for-loop as an if statement? The else only kicks in if the loop ends on it's own, not if it's broken. Pretty cool :)
                continue
            break # Breaking out of the first loop

        self.FirstPositions = self.Positions.copy()

    def EulerStep(self):
        
        # We integrate to the next timestep using the Euler method
        # Numpy lets us do this with the entire array at once!!
        self.Positions += (self.Velocities * self.dt)

    def CalculatePressure(self):

        # This method calculates pressure along one wall, which we use to check if our simulation is sound.
        Indexes = np.where(self.Positions[:,0] > self.Length/2) # Finding all indexes where particles are colliding with the wall
        Velocities = self.Velocities[:,0][Indexes] # Finding related velocities
        Forces = (Velocities * self.ParticleMass * 2) / (self.dt) # Finding the forces these particles apply 
        self.TotalPressure += sum(Forces) / ((self.Length)**2) # Summing them together
        

        if(len(Forces) > 0):
            self.counter += 1

    def CollisionStep(self):

        # Finding all indexes where the particles are outside of the box
        Indexes = np.where(abs(self.Positions) > self.Length/2)
        # Reversing all velocities where this is the case :)
        self.Velocities[Indexes] *= -1

    def TimeLoop(self):

        while self.t < self.t_max:

            self.SelectPositions[self.i] = self.Positions[0:self.NumParticles-1:int(self.NumParticles/self.NumPoints)]
            self.i += 1

            self.EulerStep()
            self.CalculatePressure()
            self.CollisionStep()

            self.t += self.dt

# Runtime code
if __name__ == "__main__":
    
    N = 10**5

    TestChamber = FuelChamber(Length = 10**(-6),Temp = 3*10**3, NumParticles = N)
    TestChamber.TimeLoop()

    # Calculating Velocity
    V = TestChamber.Velocities
    MeanV = sum((V[:,0]**2 + V[:,1]**2 + V[:,2]**2)**(1/2)) / N
    AnalyticalV = 4*((const.k_B*TestChamber.Temp)/(2 * const.pi * TestChamber.ParticleMass))**(1/2)

    # Calculating Pressure
    TotalP = TestChamber.TotalPressure
    P  = (TotalP / TestChamber.counter)
    AnalyticalP = (N * const.k_B * TestChamber.Temp) / (TestChamber.Length**3)

    # Calculating Energy
    MeanE =  ((1/2)*TestChamber.ParticleMass*(sum(V[:,0]**2 + V[:,1]**2 + V[:,2]**2)))/N
    AnalyticalE = (3/2)*const.k_B*TestChamber.Temp

    # First looking at Velocity
    print(f"Mean velocity derived from simulation          :{MeanV:.5e}")
    print(f"Mean velocity derived analyticaly              :{AnalyticalV:.5e}")
    print(f"Ratio between simulated and analytical answers :{MeanV/AnalyticalV}\n")

    # Second looking at Pressure
    print(f"Pressure calculated from simulation            :{P:.5e}")
    print(f"Pressure calculated analyticaly                :{AnalyticalP:.5e}")
    print(f"Ratio between simulated and analytical answers :{P/AnalyticalP}\n")

    # Third looking at Energy
    print(f"Mean energy calculated from simulation         :{MeanE:.5e}")
    print(f"Mean energy calculated analyticaly             :{AnalyticalE:.5e}")
    print(f"Ratio between simulated and analytical answers :{MeanE/AnalyticalE}\n")

    # Plotting some stuff

    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    line = ax.plot(TestChamber.SelectPositions[0][:,0],TestChamber.SelectPositions[0][:,1],TestChamber.SelectPositions[0][:,2],".")[0]

    def animate(i):
        line.set_data_3d(TestChamber.SelectPositions[i][:,0],TestChamber.SelectPositions[i][:,1],TestChamber.SelectPositions[i][:,2])

    ani = animation.FuncAnimation(
        fig, animate, 300, interval=100)


    plt.show()
