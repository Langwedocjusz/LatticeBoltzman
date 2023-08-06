# LatticeBoltzman

Implementation of a Lattice Boltzman (BKG) 2D9Q scheme for simulating flows.
Core simulation is written in C++17.

### Building
Regardless of the platform make sure cmake is installed and added to your path.

#### On Linux:
You can use the provided script `BuildProjects.sh`. It will ask you to choose either 'Debug' or 'Release' configuration and build the program.
Resulting executable can be found at:

	./build/bin/LatticeBoltzman
  
#### On Windows:
The provided batchfile `WIN_GenerateProjects.bat` will generate a Visual Studio solution.
After running it you can open `build/LatticeBoltzman.sln` to select configuration and build the program.

### Running
Simulation executable requires one argument: a path to a config json file:

	./build/bin/LatticeBoltzman <PATH TO A SIMULATION CONFIG FILE (JSON)>

Examples of such files can be found in:

     ./simconfigs

The program will produce a number, defined in the simulation config, of space-separated txt files containing states of the lattice (macroscopic velocities and densities) at different time steps.
They can be visualized by using the provided Python script:

    ./DrawGraphs.py

