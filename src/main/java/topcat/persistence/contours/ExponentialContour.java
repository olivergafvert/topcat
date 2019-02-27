/*
Topcat - a multidimensional persistent homology library.
Copyright (C) 2019 Oliver Gäfvert

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

import topcat.persistence.contours.kernels.ExponentialKernel;
import topcat.persistence.contours.kernels.KernelFunction;
import topcat.util.IntTuple;

import java.util.ArrayList;
import java.util.List;

public class ExponentialContour extends ProductContour {

    public ExponentialContour(List<List<Double>> filtrationIndices){
        this(null, filtrationIndices);
    }

    public ExponentialContour(List<Double> factors, List<List<Double>> filtrationIndices) {
        super(filtrationIndices, null);
        if(factors==null){
            factors = new ArrayList<>();
            for(int i=0;i<filtrationIndices.size();i++) factors.add(1.0);
        }
        this.kernels = new ArrayList<>();
        for(int i=0;i<factors.size();i++){
            this.kernels.add(new ExponentialKernel(factors.get(i)));
        }
    }
}
