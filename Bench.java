import java.awt.*;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.Random;
import java.util.Stack;

import javax.swing.*;

public class Bench extends JPanel {

	public static int BEAKER_UNIT = 10; 
	public static int BEAKER_X = 10; 
	public static int BEAKER_Y = 10; 
	public static int BEAKER_WIDTH = 100; 
	public static int BEAKER_DEPTH = 80; 
	
	public static int STANDARD_MOLES = 500; 
	
	public ArrayList<Compound> listNO2 = new ArrayList<Compound>(); 
	public ArrayList<Compound> listN2O4 = new ArrayList<Compound>(); 
	
	public Bench() {
		
		setLayout(null); 
		repaint(); 
		
		Technician.create_NO2_Gas(this, STANDARD_MOLES, 1, 1, Bench.BEAKER_WIDTH-1, Bench.BEAKER_DEPTH-1); 
		
		Timer clock = new Timer(50, new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				frame(); 
			}
		}); 
		
		JButton unfreezebutton = new JButton("Unfreeze");		
		unfreezebutton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				clock.start(); 
			}
		});
		unfreezebutton.setBounds(1200, 700, 120, 30);
		add(unfreezebutton); 
		
		JButton freezebutton = new JButton("freeze");		
		freezebutton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				clock.stop(); 
			}
		});
		freezebutton.setBounds(1500, 700, 120, 30);
		add(freezebutton); 
	}
	
	public void paintComponent(Graphics g) {
		
		super.paintComponent(g);
		setBackground(Color.LIGHT_GRAY);

		//drawing out the beaker
		g.setColor(Color.WHITE);
		g.fillRect(this.BEAKER_UNIT*this.BEAKER_X, this.BEAKER_UNIT*this.BEAKER_Y, this.BEAKER_UNIT*this.BEAKER_WIDTH, this.BEAKER_UNIT*this.BEAKER_DEPTH);
		
		//gauge panel
		g.setColor(Color.BLACK);
		g.setFont(new Font("Times New Roman", 0, 30));
		g.drawString("Count NO2: "+Integer.toString(this.listNO2.size()), 1200,100); 
		g.drawString("Count N2O4: "+Integer.toString(this.listN2O4.size()), 1200,200); 
		g.drawString("Q value: "+Double.toString((int)((this.listN2O4.size()/Math.pow(this.listNO2.size(), 2))*10000)/10000.0), 1200, 300); 
		
		Technician.image_NO2_Gas(g, this.listNO2); 
		Technician.image_N2O4_Gas(g, this.listN2O4);
		
	}
	
	public void frame() {
		
		Gas.Brownian(this.listNO2); 
		Synthesis.NO2toN2O4(this, 0.9);
		Gas.Brownian(this.listN2O4);
		SpontaneousDecomposition.N2O4toNO2(this, 0.001); 
		repaint(); 
		
	}
	
}

class Technician extends Bench {
	
	public static void create_NO2_Gas (Bench b, int moles, int minX, int minY, int maxX, int maxY) {
		
		
		double NOFixingDis = 1.197; 
		double alpha = 2.339; 
		Stack<double[]> initialState = getRandStateSet(Bench.STANDARD_MOLES, minX, minY, maxX, maxY); ; 
		
		
		for (int i = 0; i<moles; i++) {
			
			Compound thisNO2 = new Compound(3); 
			
			double[] randState = initialState.pop(); 
			
			double randX = randState[0]; 
			double randY = randState[1]; 
			double randTheta = randState[2]; 
			
			Atom N = new Atom(14.0067, 0.5, new PVector(randX, randY), 7, 5, 0.0, 3.04, 1, thisNO2); 
			Atom O1 = new Atom(15.999, 0.5, new PVector(randX+NOFixingDis*Math.cos(randTheta), randY+NOFixingDis*Math.sin(randTheta)), 8, 6, 0.0, 3.44, -1, thisNO2); 
			Atom O2 = new Atom(15.999, 0.5, new PVector(randX+NOFixingDis*Math.cos(randTheta+alpha), randY+NOFixingDis*Math.sin(randTheta+alpha)), 8, 6, 0.0, 3.44, -1, thisNO2); 
			
			Atom[] componentsNO2 = {N, O1, O2}; 
			thisNO2.set_components(componentsNO2);
			thisNO2.set_mass(46.0047);
			thisNO2.set_pos(N.get_pos());
			
			Intramolecular N_O_1 = new Intramolecular(N, O1, -1, -1, NOFixingDis, 0); 
			Intramolecular N_O_2 = new Intramolecular(N, O2, -1, -1, NOFixingDis, 0); 
			
			ArrayList<Intramolecular> bondsNO2 = new ArrayList<Intramolecular>(); 
			bondsNO2.add(N_O_1); 
			bondsNO2.add(N_O_2); 
			thisNO2.set_bonds(bondsNO2);
			thisNO2.set_dipole(new PVector(0,0));
	
			b.listNO2.add(thisNO2); 
			
			}
		
	}//end create_NO2_gas
	
	public static void image_NO2_Gas(Graphics g, ArrayList<Compound> listNO2) {
		
		for (int i =0; i<listNO2.size(); i++) {
			
			listNO2.get(i).image(g);
		}
		
	}
	
	public static void image_N2O4_Gas(Graphics g, ArrayList<Compound> listN2O4) {
		
		for (int i = 0; i<listN2O4.size(); i++) {
			
			listN2O4.get(i).image(g);
		}
		
	}
	
	public static Stack<double[]> getRandStateSet(int moles, int minX, int minY, int maxX, int maxY) {
		
		Random rand = new Random(); 
		Stack<double[]> sysState = new Stack<double[]>(); 
		
		for (int i = 0; i<moles; i++) {
			
			double randX=minX+(maxX-minX)*rand.nextDouble();
			double randY=minY+(maxY-minY)*rand.nextDouble(); 
			double randTheta=2*Math.PI*rand.nextDouble(); 
			double[] partState = {randX, randY, randTheta}; 
			sysState.push(partState); 
			
		}
		
		return sysState; 
		
	}
	
}//end Technician
