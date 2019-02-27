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

package topcat.persistence.contours.kernels;

import java.util.List;

public class StepKernelFunction extends KernelFunction {
    double[] theta;
    double[] intervals;

    public StepKernelFunction(int n, double[] intervals){
        this.intervals = intervals;
        this.theta = new double[n];
        for(int i=0;i<n;i++){
            this.theta[i] = 1;
        }
    }

    public StepKernelFunction(double[] theta, double[] intervals){
        this.theta = theta;
        this.intervals = intervals;
    }

    public StepKernelFunction(List<Double> theta, List<Double> intervals){
        this.theta = new double[theta.size()];
        this.intervals = new double[intervals.size()];
        for(int i=0;i<theta.size();i++){
            this.theta[i] = theta.get(i);
            this.intervals[i] = intervals.get(i);
        }
    }

    @Override
    public double integrate(double v, double t){
        if (t == 0) {
            return v;
        }
        int start = -1;
        for (int i = 0; i < intervals.length && start < 0; i++) {
            if (v < intervals[i]) {
                start = i;
            }
        }
        if (start < 0) {
            return t / theta[theta.length - 1] + v;
        }

        double rho_int = theta[start] * (intervals[start] - v);
        while (rho_int < t && start < intervals.length - 1) {
            start++;
            rho_int += theta[start] * (intervals[start] - intervals[start - 1]);
        }
        if (rho_int < t) {
            return (t - rho_int) / theta[theta.length - 1] + intervals[intervals.length - 1];
        }
        return intervals[start] - (rho_int-t) / theta[start];
    }
}
