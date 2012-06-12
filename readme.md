
Demo Fluids Simulation
==============================



*-- Description --*

This is a basic Eulerian fluid simulator.
The method is based on the first chapters of Fluid Simulation for Computer Graphics by Robert Bridson.
Programmation is by Nadine Dommanget (ndommanget@gmail.com).



*-- Compilation --*

This code now works on MacOSX Lion.
A few modifications are necessary on Linux and Windows.

It needs :
- OpenGL 3.2 core profile,
- SFML,
- ppm2tiff,
- ffmpeg.

It additionnally uses Glew, Eigen and the Mersenne-Twister random generator.

'cd ./src'

'make'



*-- Execution --*

'cd ./src'

'./Demo [-c <Cells number on width>] [-p <particles density [1-3]>] [-r]'

Examples : 
- './Demo -c 20 -p 1' generates a grid of [20x10x6], partially filled with 9600 particles.
- './Demo -p 2 -c 40 -r' generates a grid of [40x20x13], partially filled with 665600 particles, and records each frame in ./images/.



*-- Animation creation  --*

Last execution mandatory option : '-r'.

'sh video [name]'

Example :
- 'sh video testDemo' creates testDemo.mp4 in ./videos/.

