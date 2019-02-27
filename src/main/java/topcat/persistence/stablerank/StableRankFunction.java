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

package topcat.persistence.stablerank;

import topcat.util.Pair;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Represents the stable rank function of a persistence module for some persistence contour.
 * See the paper Stable Invariants for Multidimensional Persistence by Gäfvet and Chachólski
 * (arXiv:1703.03632) for details.
 */
public class StableRankFunction extends ArrayList<Pair<Double, Integer>> {

    /**
     * Returns the critical epsilon values of the feature counting function.
     * @return
     */
    public double[] getEpsilons(){
        double[] epsilons = new double[size()];
        for(int i=0;i<size();i++){
            epsilons[i] = get(i)._1();
        }
        return epsilons;
    }

    /**
     * Returns the values of the feature counting function at the critical epsilon
     * values.
     * @return
     */
    public int[] getValues(){
        int[] values = new int[size()];
        for(int i=0;i<size();i++){
            values[i] = get(i)._2();
        }
        return values;
    }

    @Override
    public boolean equals(Object o){
        if(o==null || !o.getClass().equals(getClass())){
            return false;
        }
        StableRankFunction other = (StableRankFunction) o;
        if(size() != other.size()){
            return false;
        }
        for(int i=0;i<size();i++){
            if(Math.abs(get(i)._1()-other.get(i)._1()) > 1E-15 ||
                    !get(i)._2().equals(other.get(i)._2())){
                return false;
            }
        }
        return true;
    }

    public List<List<Double>> toList(){
        List<Double> epsilons = new ArrayList<>();
        List<Double> values = new ArrayList<>();
        for(int i=0;i<size();i++){
            epsilons.add(get(i)._1());
            values.add(get(i)._2().doubleValue());
        }
        List<List<Double>> out = new ArrayList<>();
        out.add(epsilons);out.add(values);
        return out;
    }

    /**
     * Computes the shift to interleave f into g.
     * @param f
     * @param g
     * @return
     */
    private static double computeInterleaving(StableRankFunction f, StableRankFunction g){
        double f_e = -1;
        for(Pair<Double, Integer> f_pair : f){
            double min = Double.POSITIVE_INFINITY;
            for(Pair<Double, Integer> g_pair : g){
                if(f_pair._2() >= g_pair._2()){
                    min = g_pair._1()-f_pair._1();
                    break;
                }
            }
            f_e = min > f_e ? min : f_e;
        }
        return f_e;
    }

    /**
     * Computes the interleaving distance of two basic barcodes.
     * @param f
     * @param g
     * @return
     */
    public static double interleavingDistance(StableRankFunction f, StableRankFunction g){
        double f_e = computeInterleaving(f, g);
        double g_e = computeInterleaving(g, f);
        return f_e < g_e ? g_e : f_e;
    }

    /**
     * Computes the maximum interleaving distance over each dimension of two basic barcode collections.
     * @param f
     * @param g
     * @return
     */
    public static double interleavingDistance(StableRankFunctionCollection f, StableRankFunctionCollection g){
        if(f == null || g == null || f.getPresentDimensions() == null || g.getPresentDimensions() == null) return Double.MAX_VALUE;

        Set<Integer> dimensions = new HashSet<Integer>();
        try {
            dimensions.addAll(f.getPresentDimensions());
            dimensions.addAll(g.getPresentDimensions());
        }
        catch(Exception e){
            return Double.MAX_VALUE;
        }

        double epsilon = -1;
        for(Integer dim : dimensions){
            double t_epsilon = interleavingDistance(f.get(dim), g.get(dim));
            if(t_epsilon > epsilon){
                epsilon = t_epsilon;
            }
        }
        return epsilon;
    }
}
