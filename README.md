# Accretion-efficiency
The required code for Accretion efficiency for high redshift central SMBHs in quasars using the QUOTAS and QuasarNET datasets. The flux and volume limited (corrected) data were imported using the [Flux and Volume limited](https://github.com/artakh10/Flux-and-Volume-limit) code with a bin width of 2.5.

## Getting Started

These instructions will guide you in setting up the project on your local machine for development and testing. Please refer to the **Code Description** section for instructions on how to deploy it on a live system.

## Code Description
Here is a description of the different programs and files used in this repository:


* ```On the Accretion efficiency of High redshift Quasars.ipynb:``` runs the code for the original QUOTAS and QuasarNET data and afterward the corrected or flux and volume limited dataset using the Friends-of-Friends (FOF) Algorithm; and compares plots of accretion efficiency Vs. parameters such as redshift and BH mass. 

* ```constants_Acc.py:``` Relevant cosmological and physical constants.

* ```function_Acc.py:``` Functions used in the ```On the Accretion efficiency of High redshift Quasars.ipynb``` program based on bolometric luminosity, redshift, BH mass, optical luminosity, Eddington luminosity, and accretion efficiency.


## Requisites
This code makes use of several Python libraries:

* ```NumPy```
* ```Matplotlib```
* ```Pandas```
* ```SciPy```
*  ```Seaborn```
## Authors

* **Arta Khosravi** - [artakh10](https://github.com/artakh10)


## Citation
If you use the code, please link this repository, and cite [the related arXiv article](http://arxiv.org/abs/2405.03240).

## Contact
For comments, questions, etc. you can reach me at artakh10@gmail.com.
