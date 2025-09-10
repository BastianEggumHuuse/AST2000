# BRUKER IKKE KODEMAL!!!!
# Skrevet av Bastian Eggum Huuse og Bendik Thune

# imports
import matplotlib.pyplot as plt
# AST imports
import ast2000tools.constants as const
import ast2000tools.utils as utils
# From imports
from NozzleChamber import NozzleChamber
from ast2000tools.solar_system import SolarSystem
from ast2000tools.space_mission import SpaceMission

# Ast init
seed = utils.get_seed('bmthune')
system = SolarSystem(seed)
mission = SpaceMission(seed)

class FuelRocket:

    def __init__(self,FuelMass,SpeedBoost,NumMotors,NumParticles = 10**5,dt = 10**(-3)):

        # Parameters
        self.FuelMass = FuelMass
        self.RocketMass = mission.spacecraft_mass
        self.TotalMass = FuelMass + self.RocketMass
        self.SpeedBoost = SpeedBoost
        self.NumMotors = NumMotors
        self.dt = dt
        self.t = 0

        # We only care about velocity
        self.Velocity = 0
        self.counter = 0

        # Motor Parameters
        self.MotorLength = 10**(-6)
        self.NozzleLength = self.MotorLength * 0.25
        self.Temp = 3000
        self.NumParticles = NumParticles

        self.Thrust, self.FuelConsumption = self.SimulateEngine()

        # Lists used to save data for plotting later
        self.VelocityList = []
        self.FuelList = []

    def SimulateEngine(self):

        # Initializing motor
        print("Initializing motor...")
        self.Motor = NozzleChamber(self.MotorLength,self.Temp,self.NumParticles,self.NozzleLength)
        
        # looping through the motor for a small period of time to calculate
        # Fuel consumption and the force of the motor.
        t = 0
        MomentumSum = 0
        ParticleSum = 0
        while t < self.Motor.t_max:

            Momentum,LeavingParticles = self.Motor.TimeStep()
            MomentumSum += Momentum
            ParticleSum += LeavingParticles

            t += self.Motor.dt
            if(int((t/self.Motor.t_max)* 1000) % 100 == 0):
                print(f"Simulating motor: {int((t/self.Motor.t_max)*100):4}%")

        # Calculating the force and fuelconsumption of one motor
        Force = (MomentumSum/self.Motor.t_max) 
        FuelConsumption = ((ParticleSum/self.Motor.t_max) * self.Motor.ParticleMass)

        # We have N motors, so to simulate this, we multiply by NumMotors
        # We assume that all the motors operate the same.
        TotalForce = Force * self.NumMotors
        TotalFuelConsumption = FuelConsumption * self.NumMotors

        print(f"Finished Simulating {self.NumMotors} motors.")
        print(f"Calculated Force per motor : {Force:.5e}, Calculated Fuel Consumption per motor : {FuelConsumption:.5e}")
        print(f"Calculated Force           : {TotalForce:.5e}, Calculated Fuel Consumption      : {TotalFuelConsumption:.5e}")
        return(TotalForce,TotalFuelConsumption)

    def TimeStep(self):

        # Eulering Velocity, Fuel, and time
        self.Velocity += ((self.Thrust/self.TotalMass)) * self.dt
        self.FuelMass -= self.FuelConsumption * self.dt
        self.TotalMass = self.FuelMass + self.RocketMass
        self.counter += 1
        self.t += self.dt

        # Adding to lists (Comment out to increase performance)
        self.VelocityList.append(self.Velocity)
        self.FuelList.append(self.FuelMass)

        # Printing info every 1000th timestep (uncomment to see info)
        if self.counter % 1000 == 0:
            pass
            #print(f"Current Velocity : {self.Velocity:.3f}, Current Fuel Mass : {self.FuelMass:.3f},Current time in seconds : {self.t:.1f}, Current time in minutes : {self.t/60:.1f}")

    def TimeLoop(self):

        while(self.Velocity < self.SpeedBoost):
            self.TimeStep()
            #break

            if(self.FuelMass <= 0):
                print("HOUSTON WE HAVE A PROBLEM.... \nBAAANG")
                break
        
if __name__ == "__main__":
    

    # Creating rocket instance
    NumMotors = int((1000000**3)/60) # 1 cube meter grid / 60 :)
    Fuel = 200000
    EscapeVelocity = ((2 *(system.masses[0]*const.m_sun)*const.G) / (system.radii[0] * 1000))**(1/2)
    Particles = 10**5

    TestRocket = FuelRocket(FuelMass=Fuel,SpeedBoost=EscapeVelocity,NumMotors=NumMotors,NumParticles=Particles)
    TestRocket.TimeLoop()

    print("\nThe rocket has reached escape velocity!!! (or crashed)")
    print(f"Total fuel consumed = {(Fuel - TestRocket.FuelMass):.2f} kg")

    # Plotting
    fig, ax = plt.subplots()
    ax.plot(TestRocket.VelocityList,TestRocket.FuelList,color = "crimson")
    ax.axvline(EscapeVelocity,linestyle = "dashed",color = "k")
    ax.text(EscapeVelocity - 2500,100000,f"{EscapeVelocity:.1f} m/s")

    plt.title("Drivstoff mot hastighet")
    plt.xlabel("Hastighet [m/s]")
    plt.ylabel("Resterende drivstoff [kg]")

    # 189333.25 kg
    plt.show()