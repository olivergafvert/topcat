TOPCAT
===
Topcat is Java library for computing and analyzing multiparameter persistence modules. It comes with a python interface and a jupyter notebook with examples of how to use it. 

Input to Topcat can be given as either:
* a list of distance matrices and a list of filtration values for each parameter
* a multifiltered simplicial complex

The current implementation of Topcat has been inspired by Ulrich Bauers work on Ripser (http://ripser.org), where he uses the following optimizations (among others):
* computes the cohomology instead of homology,
* uses the combinatorial number system to index simplices,
* performs an implicit matrix reduction of the coboundary matrix .

We make heavy use of these techniques in Topcat to be able to handle millions of simplices. 

---
A detailed description of the algorithms implemented in Topcat to compute the multiparameter persistence modules can be found in my <a href="http://kth.diva-portal.org/smash/record.jsf?pid=diva2%3A939842&dswid=-4958">master's thesis</a>. A description of the theory of persistence contours, noise systems and the stable rank can be found in the following papers: 
* "Multidimensional Persistence and Noise" by Scolamiero et al. (arXiv:1505.06929)
* "Stable Invariants for Multidimensional Persistence" by G and Chachólski (arXiv:1703.03632)

---
Future work:
* Make use of the Chunk reduction described in the paper ''Chunk Reduction for Multi-Parameter Persistent Homology'' by Fugacci and Kerber (arXiv:1812.08580),
* Implement multiparameter persistence landscapes as described in the paper ''Multiparameter Persistence Landscapes'' by Oliver Vipond (arXiv:1812.09935).

Requirements
---

* Maven 2
* Java 7 or higher
* Python 2.7 (optional)
* py4j (optional)


Installation
---

In the project root folder, compile the java code with Maven using the following:

~~~
mvn clean install
~~~

Usage
---

The library can be accessed either within java or via the Python-interface.


