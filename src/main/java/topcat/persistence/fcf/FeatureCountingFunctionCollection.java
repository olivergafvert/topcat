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

package topcat.persistence.fcf;

import topcat.util.Pair;

import java.util.HashMap;
import java.util.List;
import java.util.Set;

/**
 * A collection that holds the basic barcodes for each dimension.
 */
public class FeatureCountingFunctionCollection {

    private HashMap<Integer, FeatureCountingFunction> basicBarcodes = new HashMap<>();

    public void set(Integer dim, FeatureCountingFunction featureCountingFunction){
        basicBarcodes.put(dim, featureCountingFunction);
    }

    public FeatureCountingFunction get(Integer dim){
        if(!basicBarcodes.containsKey(dim)){
            FeatureCountingFunction ret = new FeatureCountingFunction();
            ret.add(new Pair<>(0.0, 0));
            return ret;
        }
        return basicBarcodes.get(dim);
    }

    public Set<Integer> getPresentDimensions(){
        return basicBarcodes.keySet();
    }

    @Override
    public String toString(){
        StringBuilder sb = new StringBuilder();
        for(Integer dim : basicBarcodes.keySet()) {
            sb.append("Dim: ").append(dim).append(" [");
            List<Pair<Double, Integer>> basicBarcode = basicBarcodes.get(dim);
            for (Pair<Double, Integer> pair : basicBarcode){
                sb.append("(").append(pair._1()).append(", ").append(pair._2()).append("), ");
            }
            sb.append("]\n");
        }
        return sb.toString();
    }

}
