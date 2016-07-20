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

package topcat.persistence.noise;

import topcat.persistence.fcf.FeatureCountingFunction;
import topcat.persistence.functor.Functor;
import topcat.util.IntTuple;
import topcat.util.Pair;
import topcat.util.Tuple;

import java.util.ArrayList;
import java.util.List;

/**
 * Abstract noise.
 */
public abstract class Noise {

    public static List<IntTuple> getIndexSequence(List<List<Double>> filtrationValues){
        List<Double> weightList = new ArrayList<>();
        for(int i=0;i<filtrationValues.size(); i++){
            weightList.add(1.0);
        }
        Tuple<Double> weight = new Tuple<>(weightList);
        return getIndexSequence(filtrationValues, weight);
    }

    /**
     * Returns the filtration indices in N^r of a ray in the direction given by 'weight' in Q^r.
     * @param filtrationValues
     * @param weight
     * @return
     */
    public static List<IntTuple> getIndexSequence(List<List<Double>> filtrationValues, Tuple<Double> weight){
        List<IntTuple> indices = new ArrayList<>();
        IntTuple size = new IntTuple(filtrationValues.size());
        for(int i=0;i<filtrationValues.size();i++){
            size.set(i, filtrationValues.get(i).size());
        }
        IntTuple current = new IntTuple(filtrationValues.size());
        indices.add(new IntTuple(current));
        while(true){
            List<Pair<Integer, Double>> pairs = new ArrayList<>();
            double min = Double.MAX_VALUE;
            for(int v=0; v<current.length(); v++){
                if(current.get(v) < size.get(v)-1) {
                    Pair<Integer, Double> pair = new Pair<>(v, weight.get(v)*filtrationValues.get(v).get(current.get(v)+1));
                    if(pair._2() < min){
                        min = pair._2();
                    }
                    pairs.add(pair);
                }
            }
            if(pairs.size()== 0){
                break;
            }
            for(Pair<Integer, Double> v : pairs){
                if(v._2() <= min){
                    current.set(v._1(), current.get(v._1())+1);
                }
            }
            indices.add(new IntTuple(current));
        }
        return indices;
    }

    public List<IntTuple> getDiagonalIndexSequence(List<List<Double>> filtrationValues){
        List<IntTuple> indices = new ArrayList<>();

        IntTuple current = new IntTuple(filtrationValues.size());
        IntTuple unit = new IntTuple(filtrationValues.size());
        IntTuple size = new IntTuple(filtrationValues.size());
        for(int i=0;i<filtrationValues.size();i++){
            current.set(i, 0);
            unit.set(i, 1);
            size.set(i, filtrationValues.get(i).size());
        }

        while(current.lt(size)){
            indices.add(current);
            current = current.plus(unit);
        }
        return indices;
    }


    public FeatureCountingFunction computeFCF(Functor F, List<List<Double>> filtrationValues) {
        List<Double> weightList = new ArrayList<>();
        for(int i=0;i<filtrationValues.size(); i++){
            weightList.add(1.0);
        }
        Tuple<Double> weight = new Tuple<>(weightList);
        return computeFCF(F, filtrationValues, weight);
    }

    public abstract FeatureCountingFunction computeFCF(Functor F, List<List<Double>> filtrationValues, Tuple<Double> weight);
}
