package topcat.matrix;

import topcat.persistence.simplex.Simplex;

public class GradedColumn<T> extends Column<T> implements Comparable<GradedColumn<T>>{
    private Simplex grade;

    public GradedColumn(Simplex grade){
        super();
        this.grade = grade;
    }

    public Simplex getGrade(){
        return grade;
    }

    public void setGrade(Simplex grade){
        this.grade = grade;
    }

    @Override
    public int compareTo(GradedColumn<T> o) {
        return this.grade.compareTo(o.getGrade());
    }

    @Override
    public boolean equals(Object o){
        if(!o.getClass().equals(GradedColumn.class)) return false;
        GradedColumn<T> column = (GradedColumn<T>) o;
        return column.grade.equals(this.getGrade());
    }

    @Override
    public String toString(){
        return "Grade: "+grade+", "+super.toString();
    }
}
