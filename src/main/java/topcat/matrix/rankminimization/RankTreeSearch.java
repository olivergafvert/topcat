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

package topcat.matrix.rankminimization;

import topcat.matrix.BMatrix;
import topcat.matrix.BVector;
import topcat.util.Pair;

import java.util.*;

/**
 * Implements depth first search for rank minimization.
 */
public class RankTreeSearch {
    List<AffineVectorSpaceIterator> vectorSpaces;

    public RankTreeSearch(List<AffineVectorSpaceIterator> vectorSpaces){
        this.vectorSpaces = vectorSpaces;
    }

    public Pair<Integer, BMatrix> findMinRank(){
        if(vectorSpaces.size() == 0){
            return new Pair<>(0, new BMatrix(0, 0));
        }
        Collections.sort(vectorSpaces, (o1, o2) -> o1.getDimension() - o2.getDimension());
        vectorSpaces.forEach(f -> System.out.print(f.getDimension()+", "));
        System.out.print("\n");

        int k=0;
        BMatrix A = new BMatrix(vectorSpaces.size(), vectorSpaces.get(0).getAmbientDimension());
        Stack<Pair<Integer, AffineVectorSpaceIterator>> chosen = new Stack<>(); //Subspaces which have a vector in A
        Stack<Pair<Integer, AffineVectorSpaceIterator>> remaining = new Stack<>(); //Subspaces that don't have a vector in A
        int minRank = Integer.MAX_VALUE;
        BMatrix Amin = null;
        for(int i=0;i<A.rows;i++){
            A.setRow(i, vectorSpaces.get(i).next());
            chosen.push(new Pair<>(i, vectorSpaces.get(i)));
        }
        while(!chosen.isEmpty()){

            //System.out.println("A: \n"+A);

            //Compute the rank of A
            int r = BMatrix.rank(A);
            if(r < minRank){
                minRank = r;
                Amin = new BMatrix(A);
            }

            if(minRank == 1) break;

            //System.out.println("Rank(A) = " + r + " (minrank = " + minRank + ")");

            //Change the inner loop vector
            if(chosen.peek()._2().hasNext()){
                A.setRow(chosen.peek()._1(), chosen.peek()._2().next());
            }else{ //If no more remaining reset and push the iterator to remaining
                Pair<Integer, AffineVectorSpaceIterator> vs = chosen.pop();
                vs._2().reset();
                A.setRow(vs._1(), new BVector(vs._2().getAmbientDimension()));
                remaining.push(vs);
            }

            //Move up the stack and replace vectors whose iterators are exceeded
            while(!remaining.empty()){
                if(chosen.empty()) break;
                if(chosen.peek()._2().hasNext()){
                    A.setRow(chosen.peek()._1(), chosen.peek()._2().next());
                    int rankA = BMatrix.rank(A);
                    while(rankA >= minRank && chosen.peek()._2().hasNext()){ //Only move down the tree if rank can be made lower than minrank
                        A.setRow(chosen.peek()._1(), chosen.peek()._2().next());
                        rankA = BMatrix.rank(A);
                    }
                    if(rankA < minRank) {
                        while (!remaining.isEmpty()) {
                            chosen.push(remaining.pop());
                            A.setRow(chosen.peek()._1(), chosen.peek()._2().next());
                            rankA = BMatrix.rank(A);
                            while(rankA >= minRank && chosen.peek()._2().hasNext()){ //Only move down the tree if rank can be made lower than minrank
                                A.setRow(chosen.peek()._1(), chosen.peek()._2().next());
                                rankA = BMatrix.rank(A);
                            }
                            if(rankA >= minRank){
                                Pair<Integer, AffineVectorSpaceIterator> vs = chosen.pop();
                                vs._2().reset();
                                A.setRow(vs._1(), new BVector(vs._2().getAmbientDimension()));
                                remaining.push(vs);
                                break;
                            }
                        }
                    }else{
                        Pair<Integer, AffineVectorSpaceIterator> vs = chosen.pop();
                        vs._2().reset();
                        A.setRow(vs._1(), new BVector(vs._2().getAmbientDimension()));
                        remaining.push(vs);
                    }
                }else{
                    Pair<Integer, AffineVectorSpaceIterator> vs = chosen.pop();
                    vs._2().reset();
                    A.setRow(vs._1(), new BVector(vs._2().getAmbientDimension()));
                    remaining.push(vs);
                }
            }
        }
        return new Pair<>(minRank, Amin);
    }

    public static void main(String[] args){
        BVector v = new BVector(5, new int[]{0, 2, 4});
        List<BVector> vectors = Arrays.asList(
                new BVector[]{new BVector(5, new int[]{1}),
                        new BVector(5, new int[]{3})});

        List<BVector> vectors2 = Arrays.asList(
                new BVector[]{new BVector(5, new int[]{1})});

        AffineVectorSpaceIterator it1 = new AffineVectorSpaceIterator(new BMatrix(vectors), v, 4);
        AffineVectorSpaceIterator it2 = new AffineVectorSpaceIterator(new BMatrix(vectors2), v, 4);
        List<AffineVectorSpaceIterator> iterators = new ArrayList<>();
        iterators.add(it1);
        iterators.add(it2);
        RankTreeSearch rt = new RankTreeSearch(iterators);
        int minrank = rt.findMinRank()._1();

        System.out.println("Minimal rank: "+minrank);
    }
}
