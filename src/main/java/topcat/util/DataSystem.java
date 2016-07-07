/*
Topcat - a multidimensional persistent homology library.
Copyright (C) 2016 Oliver Gäfvert

This file is part of Topcat.

Topcat is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Topcat is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

package topcat.util;

import topcat.matrix.distancematrix.ArrayDistanceMatrix;
import topcat.matrix.distancematrix.DistanceMatrix;

import java.util.ArrayList;
import java.util.List;

/**
 * A data system is a collection of the information about a data set.
 * It contains a set of distance matrices and the filtration values.
 */
public class DataSystem {
    private List<DistanceMatrix> distanceMatrices = new ArrayList<>();
    private List<List<Double>> filtrationValues = new ArrayList<>();

    public DataSystem(){}

    public DataSystem(List<DistanceMatrix> distanceMatrices, List<List<Double>> filtrationValues){
        this.distanceMatrices = distanceMatrices;
        this.filtrationValues = filtrationValues;
    }

    public void addDistanceMatrix(DistanceMatrix distanceMatrix){
        this.distanceMatrices.add(distanceMatrix);
    }

    public void addDistanceMatrix(double[][] distanceMatrix){
        this.distanceMatrices.add(new ArrayDistanceMatrix(distanceMatrix));
    }

    public void addFiltrationValues(List<Double> filtrationValues){
        this.filtrationValues.add(filtrationValues);
    }

    public void addFiltrationValues(double[] filtrationValues){
        List<Double> filtrationValueList = new ArrayList<>();
        for(int i=0;i<filtrationValues.length;i++){
            filtrationValueList.add(filtrationValues[i]);
        }
        this.filtrationValues.add(filtrationValueList);
    }

    public List<DistanceMatrix> getDistanceMatrices(){
        return distanceMatrices;
    }

    public List<List<Double>> getFiltrationValues(){
        return filtrationValues;
    }

    public static DataSystem create(){
        return new DataSystem();
    }
}
