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

import topcat.persistence.contours.kernels.KernelFunction;
import topcat.persistence.contours.kernels.StepKernelFunction;
import topcat.util.IntTuple;
import java.util.ArrayList;
import java.util.List;

/**
 * Represents standard noise in the direction of a ray (which is equivalent to
 * standard noise in the direction of a cone when the dimension is less than 3).
 * For details see the paper Stable Invariants for Multidimensional Persistence by
 * Gäfvert and Chachólski (arXiv:1703.03632).
 */
public class StandardContour extends ProductContour {

    public StandardContour(List<List<Double>> filtrationIndices){
        this(null, filtrationIndices);
    }

    public StandardContour(IntTuple direction, List<List<Double>> filtrationIndices){
        super(filtrationIndices, null);
        if(direction == null){
            direction = IntTuple.ones(filtrationIndices.size());
        }
        this.kernels = new ArrayList<>();
        for(int i=0;i<direction.length();i++){
            KernelFunction kernel = new StepKernelFunction(new double[]{direction.get(i)>0 ? 1.0/direction.get(i) : 0}, new double[]{0});
            this.kernels.add(kernel);
        }
    }
}
