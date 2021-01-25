# INCA-analyzer
Matlab script for analyzing INCA PV loops data

This is a self-made MATLAB script for analyzing pressure and volume data from the INCA system.
This script calculates the usual parameters derived from the PV analysis: Ees (also single-beat Ees with a 5th order polynomial), Ea, EDV, ESV, EF, PRSW, SW, VAC, LV dPdtmax, etc... It also analyzes diastolic funciton by the non-linear fit of EDPVR (beta index) and different tau calculations . It also provides estimation of the atrial function and contribution to the LV filling.

Please, consider that I am not an engineer, I'm a dedicated physician to critical care patients with a deep pasion for physiology.

Load the script and run it in Matlab, select Example.csv for testing and check the results. Look inside the code for more variables and calculations.

Feel free to contribute to the code. Email me at: ignaciomonge@gmail.com.

Some examples of the graphical output.

[Waveforms](docs/figures/waveforms.png)

[PV loops](docs/figures/PVloops.png)

[Single PV loop with characteristic points](docs/figures/PVloop.png)

[Single beat EES calculation](docs/figures/EESsb.png)

[LV filling analysis](docs/figures/LVfilling.png)

[Diastolic function: tau calculations (four methods)](docs/figures/taudiastolic.png)

