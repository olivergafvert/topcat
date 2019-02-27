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

import org.junit.Assert;
import org.junit.Test;
import topcat.persistence.contours.kernels.KernelFunction;
import topcat.persistence.contours.kernels.StepKernelFunction;
import topcat.util.IntTuple;
import topcat.util.Tuple;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Created by oliver on 2016-03-20.
 */
public class PersistenceContourTest {

    @Test
    public void filtrationIndexLookupTest(){
        List<List<Double>> filtrationValues = new ArrayList<>();
        List<Double> f1 = new ArrayList<>();
        List<Double> f2 = new ArrayList<>();
        for(int i=0;i<10;i++){
            f1.add(i*1.0);
            f2.add(i*2.0);
        }
        filtrationValues.add(f1);filtrationValues.add(f2);
        PersistenceContour contour = new ProductContour(filtrationValues, null);
        Assert.assertEquals(contour.filtrationIndex(1.5, 0.5), new IntTuple(1, 0));
        Assert.assertEquals(contour.filtrationIndex(1.0, 2.0), new IntTuple(1, 1));
        Assert.assertEquals(contour.filtrationIndex(5.1, 4.9), new IntTuple(5, 2));
        Assert.assertEquals(contour.filtrationIndex(11.0, 21.0), new IntTuple(9, 9));
        Assert.assertEquals(contour.filtrationIndex(10.0, 0.0), new IntTuple(9, 0));
        Assert.assertEquals(contour.filtrationIndex(10.0, 1.0), new IntTuple(9, 0));
        Assert.assertEquals(contour.filtrationIndex(10.0, 1.9), new IntTuple(9, 0));
    }

    @Test
    public void StandardContourTest(){
        List<List<Double>> filtrationValues = new ArrayList<>();
        List<Double> f1 = new ArrayList<>();
        List<Double> f2 = new ArrayList<>();
        for(int i=0;i<10;i++){
            f1.add(i*1.0);
            f2.add(i*1.0);
        }

        //Positions

        IntTuple p0 = new IntTuple(0, 0);
        IntTuple p1 = new IntTuple(1, 0);
        IntTuple p2 = new IntTuple(0, 1);
        IntTuple p3 = new IntTuple(1, 1);


        filtrationValues.add(f1);filtrationValues.add(f2);
        PersistenceContour contour = new StandardContour(filtrationValues);
        Assert.assertEquals(contour.shift(p0, 0), p0);
        Assert.assertEquals(contour.shift(p0, 1), p3);
        Assert.assertEquals(contour.shift(p1, 0), p1);
        Assert.assertEquals(contour.shift(p1, 0.5), p1);
        Assert.assertEquals(contour.shift(p1, 0.9), p1);
        Assert.assertEquals(contour.shift(p1, 1), new IntTuple(2, 1));
        Assert.assertEquals(contour.shift(p3, 0), p3);
        Assert.assertEquals(contour.shift(p3, 4), new IntTuple(5, 5));


        IntTuple direction = new IntTuple(1, 2);
        contour = new StandardContour(direction, filtrationValues);
        Assert.assertEquals(contour.shift(p0, 0), p0);
        Assert.assertEquals(contour.shift(p0, 1), new IntTuple(1, 2));
        Assert.assertEquals(contour.shift(p1, 0), p1);
        Assert.assertEquals(contour.shift(p1, 0.5), p3);
        Assert.assertEquals(contour.shift(p1, 0.9), p3);
        Assert.assertEquals(contour.shift(p1, 1), new IntTuple(2, 2));
        Assert.assertEquals(contour.shift(p3, 0), p3);
        Assert.assertEquals(contour.shift(p3, 4), new IntTuple(5, 9));
    }


    @Test
    public void StepContourTest(){
        double[] theta = new double[]{1, 2, 3, 4, 5};
        double[] intervals = new double[]{1, 2, 3, 4, 5};

        List<List<Double>> filtrationValues = new ArrayList<>();
        List<Double> f = new ArrayList<>();
        for(int i=0;i<10;i++){
            f.add((double)i);
        }
        filtrationValues.add(f);

        KernelFunction kernel = new StepKernelFunction(theta, intervals);
        List<KernelFunction> kernels = new ArrayList<>();
        kernels.add(kernel);

        PersistenceContour contour = new ProductContour(filtrationValues, kernels);

        Assert.assertEquals(contour.shift(new IntTuple(0), 1), new IntTuple(1));
        Assert.assertEquals(contour.shift(new IntTuple(0), 2), new IntTuple(1));
        Assert.assertEquals(contour.shift(new IntTuple(0), 3), new IntTuple(2));
        Assert.assertEquals(contour.shift(new IntTuple(1), 2), new IntTuple(2));
        Assert.assertEquals(contour.shift(new IntTuple(2), 3), new IntTuple(3));
        Assert.assertEquals(contour.shift(new IntTuple(3), 9), new IntTuple(5));
        Assert.assertEquals(contour.shift(new IntTuple(5), 4.9), new IntTuple(5));
        Assert.assertEquals(contour.shift(new IntTuple(5), 5), new IntTuple(6));
        Assert.assertEquals(contour.shift(new IntTuple(5), 5.1), new IntTuple(6));
    }
}
