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

import topcat.persistence.contours.kernels.KernelFunction;
import topcat.util.IntTuple;

import java.util.ArrayList;
import java.util.List;

public class ProductContour extends PersistenceContour{
    List<KernelFunction> kernels;

    public ProductContour(List<List<Double>> filtrationValues, List<KernelFunction> kernels) {
        super(filtrationValues);
        this.kernels = kernels;
    }

    @Override
    public IntTuple shift(IntTuple position, double epsilon) {
        List<Double> v = new ArrayList<>();
        List<Double> fposition = filtrationValue(position);
        for(int i=0;i<fposition.size();i++){
            v.add(kernels.get(i).integrate(fposition.get(i), epsilon));
        }
        return filtrationIndex(v);
    }
}
