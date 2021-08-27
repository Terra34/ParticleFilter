# ParticleFIlter
Particle filter as part of a master thesis project

Part of a particle filter implementation on ARM Cortex-M0 core. 
The C++ version was run on Cortex-A9 core (Zybo Z7-20) and performed better(faster) than the C version.

Several pseudorandom gaussian generators were implemented (with python API):
  1. Box-Muller transformation
  2. Marsaglia's polar method
  3. Inverse method
  4. Ziggurat method

The C version is modified to be used on ARM Cortex-M1/M0 core.
Includes a fixed point implementation. 
This implementation runs about 10 times faster than the fastest version of the floating point filter while also having a reduced spatial complexity.
