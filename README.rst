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

Standard system 1
=================

Paremeters at: https://confluence.esrf.fr/display/PROJEBSL1/Parameters+for+hybrid+ray+tracing+and+wave+propagation

.. list-table:: Wofry Results
   :widths: 30 25 25 25 25
   :header-rows: 1

   * - How?
     - H FWHM [um]
     - V FWHM [um]
     - Coherent Fraction H
     - Coherent Fraction V
   * - Wofry 1D GSM Gausian slit
     - 19.8
     - 19.9
     - 0.83
     - 0.84
   * - Wofry 1D Und Rect slit
     - 16.4
     - 24.6
     - 0.90
     - 0.92
     
Note that coherent fraction refers to the final beam after all optics (not at the source).

Images at sample position for Und with Rect slit H and V:

.. image:: https://github.com/srio/paper-transfocators-resources/blob/main/Figures/Figure_System1_H.png

.. image:: https://github.com/srio/paper-transfocators-resources/blob/main/Figures/Figure_System1_V.png


The differences between GSM and Und+RectSlit come from both the source and the slit, but it seems that theeffect of the slit is more important. This is a Wofry 2D calculation for the fist GSM mode using a Rectangular slit:

.. image:: https://github.com/srio/paper-transfocators-resources/blob/main/Figures/Figure_System1_2D_GSM.png

TODO
====
- compute transmission (intensity at center and integrated intensity)
- calculate f2 position at 190
- Calculate all wanted energies

DONE
====
- Obtain GSM trajectories (f1,f2) and corresponding size for better comparison
- Double-check the papameters used in GSMM and WOFRY
- Recalculate using undulator coherent mode decomposition, rectangular slit (Gaussian used) and real lenses (ideal used)
- Select parameters for 2D simulations and compare
