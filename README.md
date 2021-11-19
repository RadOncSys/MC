# MC

Monte-Carlo project for dosimetry modelling in the field of radiation oncology.

History
=======

This project has started from the 
[EGS-Nova: An Adaptation of EGS in C/C++](http://rcwww.kek.jp/research/egs/epub/aap/js3nov98.html)
many years ago.

Thanks a lot to James C. Satterthwaite who has ported EGS4 to well structured C language application.

***MC*** package was evolving in the process of many tasks solving in the field of radiotherapy dosimetry.
Among them simulations of pencil beam dose distributions, radiotherapy unit radiation beam forming systems simulations, radiation source modeling, etc.

Among them are simulations of pencil beam dose distributions, radiotherapy unit radiation beam forming systems simulations, radiation source modeling, etc.

Project motivation
==================

***EGS4*** is simple and attractive radiation transport simulation tool. Unfortunately, it is almost impossible to manage this system in modern software business world. Other Monte Carlo solutions (***NRCEGS*** is still ***EGS*** with ***Mortran*** in core) are less benchmarked and friendly to radiotherapy community.

***MC*** was not designed to simply port ***EGS4***. Porting was perfectly done in ***EGS-Nova***. For unclear reasons NOVA disappeared from the internet some years ago. ***MC*** development was originally driven by the necessity of radiation source modelling and simplifying complex geometry scenes programming. So, it is not as flexible as ***EGS4***.

Important reason for the creating new ***MC*** engine was requirement of programming new complex geometries.

The last but not least reason is that ***C++*** classical object oriented MC system is a good educational and self-educational tool for studying radiation transport field and improving programming skills.

Important decisions
===================

The biggest philosophical difference compare to ***EGS4*** is assumption, that scene objects form a linear chain. Particle leaving any object can hit only next object if moved in positive Z direction or previous object if moved in negative direction. This is strong restriction. There are two mechanisms to work around. One is programming complex objects as a separate object with its own logic like conventional ***EGS4***. The other work around is supporting nested objects.

Transport logic (managing events and their processing) keeps physics nature from ***EGS4***, but has completely different implementation. It may be difficult to understand. Advantage is all this logic is implemented in compact form in a single class ***mcTransport***.

Physics
=======

***MC*** package implements ***EGS4*** physics models without modifications.

Feature directions
==================

- Replace PEGS4 data file with cross-sections generation on the fly.
- Verify and improve radiation transport physics.
- Documentation.
- More samples.

Installation
===================

Project assume to build under ***Visual Studio 2017 Community edition***.

To build it just clone or download it from the ***GitHub*** and open solution.

Solution consists of four projects.

| Project | Description |
| :---- | :---- |
| ***Geometry***  |Small internal use geometry base classes, like vectors, matrixes, planes, etc |
| ***MC***  | Core Monte Carlo library |
| ***MC.Tests***  | Some collection of tests. Valuable learning resource |
| ***MCSimulator***  | Complex executable, which uses full power of MC library and can be used standalone. Most practical tasks can be completely described by XML input files |

Running simulations
===================

***MC*** library can be used to program custom tasks. Source code for ***MCSimulator*** fully demonstrates how to do that and can be considered as reference implementation. It is very like ***EGS4*** workflow.

***MCSimulator*** can be used as ready to use *Monte Carlo* simulator. It takes two *XML* input files to describe the task. First file describes simulation conditions. Second file describes simulation geometry.

***Data/XML*** folder contains some input examples.

To run simulations you need to create a folder for the task and copy there ***MCSimulator.exe*** and input files. Then you need to have side by side to this folder data folder with ***PEGS4*** input file (***AcceleratorSimulator.pegs4dat*** for now).

Here are some examples to run simulations from the samples data:

	MCSimulator S_C60_sphere.xml G_C60_sphere.xml
	MCSimulator S_CK_6X.xml G_CyberKnife.xml
	MCSimulator S_Ross.xml G_Ross.xml
	MCSimulator S_Target.xml G_Target.xml

Output and visualization
========================

***MCSimulator*** produces two output files. Their names indicated in simulation conditions ***XML*** file. One file (**.dat*) contains numerical results of simulations. The other one is VRML graphical dump (**.wrl*). It can be viewed many ***VRML*** viewers. The recommended old style perfect one can be unzipped from this project ***.\Data\vrmlview.zip*** file.

Documentation
=============

Project documentation can be found here
[here...](https://radoncsys.github.io/MC)

References
==========

 1. [MC code presentaion at meetup](https://radoncsys.github.io/MC/files/2017_05_23_Meetup_MC.pdf) (Russian)
 1. [MC validation presentaion](https://radoncsys.github.io/MC/files/20201021_MCValid.pdf) (Russian)