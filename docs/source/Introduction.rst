Introduction
============

Purpose of the program
^^^^^^^^^^^^^^^^^^^^^^

**MC** is an ionizing radiation Monte-Carlo transport library and reference application.

Primarily **MC** is Monte Carlo project for dosimetry modelling in the field of radiation oncology.

History
^^^^^^^

This project has started from the 
`EGS-Nova: An Adaptation of EGS in C/C++ <http://rcwww.kek.jp/research/egs/epub/aap/js3nov98.html>`_
many years ago.

Thanks a lot to James C. Satterthwaite who has ported EGS4 to well structured C language application.

**MC** package was evolving in the process of many tasks solving in the field of radiotherapy dosimetry.
Among them simulations of pencil beam dose distributions, radiotherapy unit 
radiation beam forming systems simulations, radiation source modeling, etc.

Among them are simulations of pencil beam dose distributions, 
radiotherapy unit radiation beam forming systems simulations, radiation source modeling, etc.

Project motivation
^^^^^^^^^^^^^^^^^^

**EGS4** is simple and attractive radiation transport simulation tool. 
Unfortunately, it is almost impossible to manage this system in modern software business world. 
Other Monte Carlo solutions (**NRCEGS** is still **EGS** with **Mortran** in core) are
less benchmarked and friendly to radiotherapy community.

**MC** was not designed to simply port **EGS4**. Porting was perfectly done in **EGS-Nova**. 
For unclear reasons NOVA disappeared from the internet some years ago. **MC** development
was originally driven by the necessity of radiation source modelling and simplifying complex geometry scenes programming. 
So, it is not as flexible as **EGS4**.

Important reason for the creating new **MC** engine was requirement of programming new complex geometries.

The last but not least reason is that **C++** classical object oriented MC system is a 
good educational and self-educational tool for studying radiation transport field and improving programming skills.

Important decisions
^^^^^^^^^^^^^^^^^^^

The biggest philosophical difference compare to **EGS4** is assumption, 
that scene objects form a linear chain. 
Particle leaving any object can hit only next object if moved in 
positive Z direction or previous object if moved in negative direction. 
This is strong restriction. There are two mechanisms to work around. 
One is programming complex objects as a separate object with its own logic like
conventional **EGS4**. The other work around is supporting nested objects.

Transport logic (managing events and their processing) keeps physics nature 
from **EGS4**, but has completely different implementation. 
It may be difficult to understand. 
Advantage is all this logic is implemented in compact form in a single class **mcTransport**.

Physics
^^^^^^^

**MC** package implements **EGS4** physics models without modifications.

Feature directions
^^^^^^^^^^^^^^^^^^

- Replace PEGS4 data file with cross-sections generation on the fly.
- Verify and improve radiation transport physics.
- Documentation.
- More samples.
