# UCN-Transport-Simulations
Upgrade and use of Adam Holley's neutron transport code for use in the NCSU UCN project.

To Do:
  -Edit Transport Simulator.c to take command line arguments for (necessary so that we do not have to move and copy the sim file everywhere):
        -location of regionfile
        -location of connexfile
        -location to output files
        -location of batch file (optional, if it exists then run a batch, if it doesn't only run a single run)
     -all other parameters can be changed in the Regionfile of Connexfile or changed temporarily in the Transport Simulator and then changed back (preferrably make a copy and then delete it afterwards)
     
  -Edit Transport simulator so that poof() can be called with a random offset from the cutplane to enable particle generation anywhere in a region: use for Leung Next Generation Source
  
  -Edit Transport simulator to enable 108neV boost in random direction on contact with a cutplane: use for simulation of frost on surface of deuterium
