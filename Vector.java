import java.util.*; 
import java.math.*; 

public class Vector {
	
	public Vector() {return; }
	
}

class PVector extends Vector {
	
	//fields and defaults***********************************************
	private double x, y;
	
	public PVector() {
		
		super(); 
		
		this.x = 0; 
		this.y = 0; 
	}
	
	public PVector(double x, double y) {
		
		super(); 
		
		this.x = x; 
		this.y = y; 
		
	}
	
	@Override 
	public String toString() {
		
		return "("+Double.toString(this.x)+", "+Double.toString(this.y)+")"; 
		
	}
	
	//*******************************************************************	
	//sets gets**********************************************************
	
	public void set_x(double x) {this.x = x; }
	public void set_y(double y) {this.y = y; }
	
	public double get_x() {return this.x; }
	public double get_y() {return this.y; }
	
	public double get_mag() {return Math.sqrt(Math.pow(x,  2)+Math.pow(y, 2)); }
	public double get_theta() {
		
		if(this.x!=0) {return Math.atan((this.y/this.x)); }
		if (this.y>0) {return Math.PI/2;}
		if (this.y<0) {return 3*Math.PI/2;}
		else {return 0; }
			
	}
	
	//*******************************************************************
	//behaviorals********************************************************
	
	public void plus(PVector v) {
		
		this.x = this.x+v.get_x(); 
		this.y = this.y+v.get_y(); 
		
	}
	
	public void minus(PVector v) {
		
		this.x = this.x-v.get_x(); 
		this.y = this.y-v.get_y(); 
		
	}
	
	//*******************************************************************
	//static behaviorals********************************************************
		
	public static PVector add(PVector v1, PVector v2) {
		
		return new PVector(v1.get_x()+v2.get_x(), v1.get_y()+v2.get_y()); 
		
	}
	
	public static PVector subtract(PVector v1, PVector v2) {
		
		PVector res = new PVector(v1.get_x()-v2.get_x(), v1.get_y()-v2.get_y()); 
		return res; 
		
	}
	
	public static double dot(PVector v1, PVector v2) {
		
		return v1.get_x()*v2.get_x()+v1.get_y()*v2.get_y(); 
		
	}
	
	public static double cross(PVector v1, PVector v2) {
		
		return v1.get_x()*v2.get_y()-v1.get_y()*v2.get_x(); 
		
	}
	
	public static PVector multiply(PVector v, double scale) {
		
		return new PVector(v.get_x()*scale, v.get_y()*scale); 
		
	}
	
}
