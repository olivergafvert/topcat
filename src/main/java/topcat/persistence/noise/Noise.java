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

import topcat.persistence.barcode.BasicBarcode;
import topcat.persistence.functor.Functor;
import topcat.util.IntPair;
import topcat.util.IntTuple;
import topcat.util.Pair;
import topcat.util.Tuple;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Abstract noise.
 */
public abstract class Noise {

    public static List<IntTuple> getIndexSequence(List<List<Double>> filtrationValues){
        Tuple<Double> weight = new Tuple<>(
                IntStream.range(0, filtrationValues.size()).mapToObj(v->1.0).collect(Collectors.toList())
        );
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
            List<Pair<Integer, Double>> pairs = IntStream.range(0, current.length())
                    .mapToObj(v -> {
                        if(current.get(v) < size.get(v)-1) {
                            return new Pair<>(v, weight.get(v)*filtrationValues.get(v).get(current.get(v)+1));
                        }
                        return null;
                    })
                    .filter(v -> v != null)
                    .collect(Collectors.toList());
            if(pairs.size()== 0) break;
            double min = pairs.stream()
                    .min((v, w) -> v._2().compareTo(w._2()))
                    .get()
                    ._2();
            pairs.stream()
                    .filter(v -> v._2() <= min)
                    .forEach(v -> current.set(v._1(), current.get(v._1())+1));
            indices.add(new IntTuple(current));
        }
        return indices;
    }

    public abstract BasicBarcode computeBasicBarcode(Functor F, List<List<Double>> filtrationValues);

}
