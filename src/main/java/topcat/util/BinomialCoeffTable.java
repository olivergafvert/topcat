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

package topcat.util;

import it.unimi.dsi.fastutil.longs.LongArrayList;
import it.unimi.dsi.fastutil.longs.LongList;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * A table used to implement the combinatorial number system to index the simplices. Inspired by
 * its use in the implementation of Ripser [1].
 *
 * [1] - http://ripser.org
 */
public class BinomialCoeffTable {
    List<LongList> B;

    public BinomialCoeffTable(int n, int k){
        this.B = new ArrayList<>();
        for(int i=0;i<=n;i++){
            LongList l = new LongArrayList();
            for(int j=0;j<=i && j<=k+1; j++){
                if(j==0 || j==i){
                    l.add(1);
                }else{
                    l.add(B.get(i-1).getLong(j-1)+B.get(i-1).getLong(j));
                }
            }
            B.add(l);
        }
    }

    public long get(int n, int k){
        if(k>n) {
            return 0;
        }
        return B.get(n).getLong(k);
    }

    /**
     * Computes the index of a k-simplex in the combinatorial number system.
     * @param vertices
     * @return
     */
    public long computeIndex(List<Integer> vertices){
        Collections.sort(vertices);
        long index = 0;
        for (int i = 0; i < vertices.size(); i++) {
            index += get(vertices.get(i), i+1);
        }
        return index;
    }

    public long computeIndex(Integer... v){
        return computeIndex(Arrays.asList(v));
    }
}
