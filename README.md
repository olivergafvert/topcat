TOPCAT
===
Topcat is Java library for computing invariants on multidimensional persistence modules. Topcat has been developed during my work on my master's thesis about multidimensional persistence at KTH Royal Institute of Technology in Stockholm, Sweden. The theoretical framework is based on the paper "Multidimensional Persistence and Noise" by Chachólski et al. (arXiv:1505.06929). A detailed description of the algorithms implemented in Topcat can be found in my <a href="http://kth.diva-portal.org/smash/record.jsf?pid=diva2%3A939842&dswid=-4958">thesis</a>.

Topcat currently only supports computations over the field Z/2Z but all algorithms that are used in Topcat extend to an arbitrary field. Also, Topcat currently only supports computation of the homology of a multifiltered simplicial complex, although the framework may be extended to computing the homology of any type of multifiltered space for which there are algorithms to do so.

The input that may be given to Topcat is in the form of a list of distance matrices together with one or more lists of filtration values. Then a multifiltered Vietoris-Rips complex will be constructed using this data and passed to the next step in the pipeline. Given this input Topcat computes the complete persistence module, including all maps, for each dimension up to a specified max dimension. 

When the persistence modules have been computed it is straightforward to compute any invariant. Topcat currently supports computation of the rank invariant (see "The theory of multidimensional persistence" by Carlsson and Zomorodian) and the feature counting function (previously called basic barcode) (see "Multidimensional persistence and noise" by Chachólski et al.) for a special type of domain noise and a lower bound of the feature counting function for standard noise in the direction of a ray.


Requirements
---

* Maven 2
* Java 7 or higher
* MATLAB (optional)


Installation
---

In the project root folder, compile the java code with Maven using the following:

~~~
mvn clean install
~~~

Usage
---

The library can be accessed either within java or via the MATLAB-interface. To access the library via MATLAB, add the topcat/matlab-folder to the matlab path with the command:

~~~
addpath(genpath('path to/topcat/matlab'))
~~~ 

Then execute the script load_topcal.m.

