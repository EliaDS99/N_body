# N-body system time evolution

Project Status: Completed (Academic Project)
Development Period: May 2024- Jul 2024
Language: C

## Description

This is a standalone C implementation of a numerical solver for the time evolution of a three-dimensional N-body system under a central gravitational interaction. The code integrates the equations of motion using a fourth-order Runge–Kutta scheme and allows the study of the dynamical evolution, stability, and energy behavior of the system.

The solver supports different initial spatial configurations and includes the option to add a massive central object, enabling the modeling of clustered systems with different geometries and centrally dominated gravitational potentials.

## How to use

The code is contained in a single file for ease of use.
Compile with any standard C compiler:

$ gcc nbody.c -o nbody.x -lm
