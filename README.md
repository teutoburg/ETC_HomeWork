# Build a HAWKI ETC

### User Story:

*Kieran wants to study the initial mass function of an open cluster in the Large Magellanic Cloud.
He plans to use the HAWK-I instrument on the VLT and would like to know what the faintest stars are that will be detectable (S/N > 5) with 1 hour of observing time in the K filter.
He assumes average (50% percentile) observing conditions will apply for the observation*

### Deliverable:
The goal of this exercise is to build a **very** simple exposure time calculator (ETC) module in Python for the [HAWK-I imager instrument on the VLT](https://www.eso.org/sci/facilities/paranal/instruments/hawki.html).
The only conditions for how to implement the ETC are:
- written in Python >=3.7
- the ETC function or class must be importable by other scripts (i.e. a function or class)
- a test suite is included in a separate file (e.g. `test_etc.py`)

There are no design constraints for the implementation of the ETC, however the code should be written in Python and include a test suite.

If you make any simplifying assumptions during the coding process, please document these in comments so that the reason for the assumption isn't lost to the sands of time. 

### Hints:
The aim of this exercise is to demonstrate coding style, not to showcase radiometric talents. 
Therefore, it' is perfectly fine to make use of the radiometric-gymnastics package: [How Many Photons](https://pypi.org/project/HowManyPhotons/)

For comparing results, you may want to use the official [ESO HAWK-I ETC](https://www.eso.org/observing/etc/bin/gen/form?INS.NAME=HAWK-I+INS.MODE=imaging).
