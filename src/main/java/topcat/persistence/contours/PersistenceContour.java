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

package topcat.persistence.contours;

import topcat.util.IntTuple;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Abstract noise.
 */
public abstract class PersistenceContour {
    protected List<List<Double>> filtrationValues;

    public PersistenceContour(List<List<Double>> filtrationValues){
        this.filtrationValues = filtrationValues;
    }

    public abstract IntTuple shift(IntTuple position, double epsilon);


    public List<Double> filtrationValue(IntTuple index){
        List<Double> values = new ArrayList<Double>();
        for(int i=0;i<index.length();i++){
            values.add(filtrationValues.get(i).get(index.get(i)));
        }
        return values;
    }

    public IntTuple filtrationIndex(List<Double> filtrationValue){
        IntTuple index = IntTuple.zeros(filtrationValue.size());
        for(int i=0;i<index.length();i++){
            int findex = Collections.binarySearch(filtrationValues.get(i), filtrationValue.get(i));
            if (findex < 0){
                findex = -findex >= filtrationValues.get(i).size() ? filtrationValues.get(i).size()-1 : -(findex+2);
            }
            index.set(i, findex);
        }
        return index;
    }

    public IntTuple filtrationIndex(Double... vals){
        return filtrationIndex(Arrays.asList(vals));
    }
}
