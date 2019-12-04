from py4j.java_gateway import (JavaGateway, GatewayParameters, launch_gateway)
import numpy as np
import os, sys

def absolutejarpath():
	filepath = os.path.abspath(__file__)
	projectdir = os.path.split(os.path.split(filepath)[0])[0]
	jarpath = os.path.join(os.path.join(projectdir, "target"), "topcat-1.0-SNAPSHOT.jar")
	return jarpath

port = launch_gateway(classpath=absolutejarpath(), die_on_exit=True, redirect_stdout=sys.stdout)
gateway = JavaGateway(gateway_parameters=GatewayParameters(port=port, auto_convert=True))
topcat = gateway.jvm.topcat.mains.PythonInterface

def persistenceModules_dist(distanceMatrices, filtrationValues, maxdim, contour=None):
	if contour == None:
		return list(map(PersistenceModule, topcat.computePersistenceModules(distanceMatrices, filtrationValues, maxdim)))
	return list(map(PersistenceModule, topcat.computePersistenceModules(distanceMatrices, filtrationValues, maxdim, contour)))

'''
	Computes the multiparameter persistence modules from a list of points up to dimension 'maxdim'.

	@param points - points a numpy array of points
	@param distances - a list of strings of the distances to be used ('euclidean', 'euclidean_codensity')
	@param filtrationValues - a numpy array of filtration values for each distance
	@param maxdim - the max dimension of the homology to be computed
	
	Returns a list of python PersistenceModule objects.
'''
def persistenceModules(points, distances, filtrationValues, maxdim):
	return list(map(PersistenceModule, topcat.computePersistenceModules(points, distances, filtrationValues, maxdim)))
	

def stableRank_dist(distanceMatrices, filtrationValues, maxdim, contour=None):
	if contour == None:
		return np.asarray(list(topcat.computeStableRank(distanceMatrices, filtrationValues, maxdim)))
	return np.asarray(list(topcat.computeStableRank(distanceMatrices, filtrationValues, maxdim, contour)))

'''
	Computes the stable rank of the multiparameter persistence modules computed from a list 
	of points up to dimension 'maxdim'.

	@param points - points a numpy array of points
	@param distances - a list of strings of the distances to be used ('euclidean', 'euclidean_codensity')
	@param filtrationValues - a numpy array of filtration values for each distance
	@param maxdim - the max dimension of the homology to be computed
	@param contour - a numpy array of values for the step contours to be used for the stable rank (one for each distance).

	Returns a list of python PersistenceModule objects.
'''
def stableRank(points, distances, filtrationValues, maxdim, contour=None):
	if contour == None:
		return np.asarray(list(topcat.computeStableRank(points, distances, filtrationValues, maxdim)))
	return np.asarray(list(topcat.computeStableRank(points, distances, filtrationValues, maxdim, contour)))


class PersistenceModule(object):
	def __init__(self, module):
		self.module = module

	'''
		Computes the stable rank of the persistence module for shift values 'values'.
		@param values - a list of floats specifying the shift values
		returns a StableRankFunction object containing the stable rank for the shift values.
	'''
	def stableRank(self, values):
		return self.module.computeStableRank(values)

	'''
		Computes the stable rank for a specified contour
		@param contour - a PersistenceContour java object
	'''
	def stableRank(self, values, contour):
		return self.module.computeStableRank(values, contour)

	'''
		Computes the rank of the map from index u to index v.
		@param u - a list of integers
		@param v - a list of integers
		returns the rank of the map.
	'''
	def rank(self, u, v):
		return self.module.rank(u, v)

	def __str__(self):
		return self.module.getFunctor().toString()


