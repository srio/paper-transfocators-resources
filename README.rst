=============================
paper-transfocators-resources
=============================

Calculations for f2 at 177m 
Slit size in the following calculations: 37 microns (H) x 205 microns (V).

Select slit opening
===================

Select the aperture to match the desired coherent fraction:

.. image:: https://github.com/srio/paper-transfocators-resources/blob/main/Figures/Figure_1cf.png
.. image:: https://github.com/srio/paper-transfocators-resources/blob/main/Figures/Figure_1cf_und.png



(f1,f2) trajectories
====================

For a given f1, calculate f2 in order to have the focus at the sample. 

.. image:: https://github.com/srio/paper-transfocators-resources/blob/main/Figures/Figure_1traj_comparison.png

Focal sizes
===========

This are the data from the two methods, each method use its own (f1,f2) trajectories, that are not the same. 

.. image:: https://github.com/srio/paper-transfocators-resources/blob/main/Figures/Figure_1size_comparison1.png

To compare the methods I calculated the sizes using the same (f1,f2) trajectory for the two codes (the one from WOFRY)

.. image:: https://github.com/srio/paper-transfocators-resources/blob/main/Figures/Figure_1size_comparison2.png

and with the same (f1,f2) trajectory for the two codes (the one from GSMM)

.. image:: https://github.com/srio/paper-transfocators-resources/blob/main/Figures/Figure_1size_comparison3.png



TODO
====
- Recalculate using undulator coherent mode decomposition, rectangular slit (Gaussian used) and real lenses (ideal used)
- calculate f2 position at 190
- Calculate all wanted energies
- compute transmission (intensity at center and integrated intensity)
- Select parameters for 2D simulations and compare

DONE
====
- Obtain GSM trajectories (f1,f2) and corresponding size for better comparison
- Double-check the papameters used in GSMM and WOFRY
