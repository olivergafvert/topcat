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

import org.junit.Assert;
import org.junit.Test;
import topcat.util.IntTuple;
import topcat.util.Tuple;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Created by oliver on 2016-03-20.
 */
public class NoiseTest {

    @Test
    public void indexSequenceTest(){
        List<List<Double>> filtrationValues = new ArrayList<>();
        filtrationValues.add(Arrays.asList(new Double[]{
                0.0, 1.0, 5.0, 10.0
        }));
        filtrationValues.add(Arrays.asList(new Double[]{
                0.0, 0.5, 1.0, 1.5
        }));
        List<IntTuple> indices = new ArrayList<>();
        indices.add(new IntTuple(0, 0));
        indices.add(new IntTuple(1, 1));
        indices.add(new IntTuple(1, 2));
        indices.add(new IntTuple(1, 3));
        indices.add(new IntTuple(2, 3));
        indices.add(new IntTuple(3, 3));
        Assert.assertEquals(Noise.getIndexSequence(
                filtrationValues, new Tuple<>(0.5, 1.0)), indices);
    }
}
