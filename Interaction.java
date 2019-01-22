import java.util.*;
import java.awt.*;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Random;

import javax.swing.*;
import java.math.*; 

public class Interaction {

	//fields and defaults***********************************************

	private ArrayList<Material> subjects; 
	private double energy; 
	
	public Interaction() {
		
		this.subjects = new ArrayList<Material>(); 
		this.energy = 0; 
		
	}
	
	public Interaction(ArrayList<Material> subjects, double energy) {
		
		this.subjects = subjects; 
		this.energy = energy; 
		
	}
	
	//*******************************************************************	
	//sets gets**********************************************************
	
	public void set_subjects(ArrayList<Material> subjects) {this.subjects = subjects; }
	public void set_energy(double energy) {this.energy = energy; }
	
	public ArrayList<Material> get_subjects() {return this.subjects; }
	public double get_energy() {return this.energy; }
	
	//*******************************************************************
	//behaviorals********************************************************
	
	
	
	//*******************************************************************
}

class Intermolecular extends Interaction {
	
	//fields and defaults***********************************************

	private Molecule molecule1; 
	private Molecule molecule2; 
	private double effectiveDis; 
	
	public Intermolecular() {
		
		super(); 
		
		this.molecule1 = null; 
		this.molecule2 = null; 
		this.effectiveDis = -1; 
		
	}
	
	public Intermolecular(Molecule molecule1, Molecule molecule2, double effectiveDis) {
		
		super(); 
		
		this.molecule1 = molecule1; 
		this.molecule2 = molecule2; 
		this.effectiveDis = effectiveDis; 
		
	}
	
	public boolean equals(Intermolecular other) {
		
		return (this.molecule1.equals(other.molecule1) && this.molecule2.equals(other.molecule2))
			|| (this.molecule1.equals(other.molecule2) && this.molecule2.equals(other.molecule1)); 
	}
	
	//*******************************************************************	
	//sets gets**********************************************************

	public void set_molecule1(Molecule molecule1) {this.molecule1 = molecule1; }
	public void set_molecule2(Molecule molecule2) {this.molecule2 = molecule2; }
	
	public Molecule get_molecule1() {return this.molecule1; }
	public Molecule get_molecule2() {return this.molecule2; }
	
	//*******************************************************************
	//behaviorals********************************************************



	//*******************************************************************
}

class HydrogenBond extends Intermolecular {
	
	//fields and defaults***********************************************
	
	public static double EFFECTIVE_DIS = 4.0; //Distances measured in Angstrom!!
	
	private Atom donor; 
	private Atom receiver; 
	private Compound donorMolecule; 
	private Compound receiverMolecule; 
	private double dis; 
	
	public HydrogenBond() {
		
		super(null, null, HydrogenBond.EFFECTIVE_DIS); 
		
		this.donor = null; 
		this.receiver = null; 
		this.donorMolecule = null; 
		this.receiverMolecule = null; 
		this.dis = -1; 
				
	}
	
	public HydrogenBond(Atom donor, Atom receiver, Compound donorMolecule, Compound receiverMolecule, double initialDis) {
		
		super(donorMolecule, receiverMolecule, HydrogenBond.EFFECTIVE_DIS); 
		
		this.donor = donor; 
		this.receiver = receiver; 
		this.donorMolecule = donorMolecule; 
		this.receiverMolecule = receiverMolecule; 
		this.dis = initialDis; 
		
	}
	
	//*******************************************************************	
	//sets gets**********************************************************

	public void set_donor(Atom donor){
		
		if (donor.get_number()==1 && donor.get_partCharge()>0) {
			this.donor = donor; 
			this.donorMolecule = donor.get_inCompound(); 
		}
		
		else {System.out.println("PHYS ERROR: "+donor.toString()+" cannot be a hydrogen bond donor.");}
		
	}
	
	public void set_receiver(Atom receiver) {
		
		if (donor.get_countLonePair()>1 && donor.get_partCharge()<0) {
			this.receiver = receiver; 
			this.receiverMolecule = receiver.get_inCompound(); 
		}
		
		else {System.out.println("PHYS ERROR: "+receiver.toString()+" cannot receive a hydrogen bond.");}
			
	}
	
		//Missing donor&receiver Molecule setters
	
	public void set_dis(double dis) {
		
		if (dis>=0) {this.dis = dis; }
		else {System.out.println("PHYS ERROR: Distance cannot be negative. ");}
		
	}
	
	public Atom get_donor() {return this.donor; }
	public Atom get_receiver() {return this.receiver; }
	public Compound get_donorMolecule() {return this.donorMolecule; }
	public Compound get_receiverMolecule() {return this.receiverMolecule; }

	//*******************************************************************
	//behaviorals********************************************************



	//*******************************************************************
}

class Coulombic extends Intermolecular {
	
	//fields and defaults***********************************************

	
	
	//*******************************************************************	
	//sets gets**********************************************************



	//*******************************************************************
	//behaviorals********************************************************



	//*******************************************************************
}

class London extends Intermolecular {
	
	//fields and defaults***********************************************

	
	
	//*******************************************************************	
	//sets gets**********************************************************



	//*******************************************************************
	//behaviorals********************************************************



	//*******************************************************************
}


class IonicBond extends Intermolecular {
	
	//fields and defaults***********************************************

	
	
	//*******************************************************************	
	//sets gets**********************************************************



	//*******************************************************************
	//behaviorals********************************************************

	

	//*******************************************************************
}

class Intramolecular extends Interaction {

	//fields and defaults***********************************************
	
	//physical properties
	private Atom atom1; 
	private Atom atom2; 
	private PVector dipole; 
	private double electrodiff; 
	
	//simulation distances
	private double bondingDis; 
	private double separatingDis; 
	private double fixingDis; 
	
	//simulation probs 
	private double separatingProb; 
	
	public Intramolecular() {
		
		super(); 
		
		this.atom1 = null; 
		this.atom2 = null; 
		this.dipole = new PVector(); 
		this.electrodiff = 0; 
		
		this.bondingDis = -1; 
		this.separatingDis = -1; 
		this.fixingDis = -1; 
		
		this.separatingProb = -1; 
		
	}
	
	public Intramolecular(Atom atom1, Atom atom2, double bondingDis, double separatingDis, double fixingDis, double separatingProb) {
		
		super(); 
		
		this.atom1 = atom1; 
		this.atom2 = atom2; 
		
		this.electrodiff = this.atom1.get_electronegativity()-this.atom2.get_electronegativity(); 
		
		this.bondingDis = bondingDis; 
		this.separatingDis = separatingDis; 
		this.fixingDis = fixingDis; 
		
		this.separatingProb = separatingProb; 
		
		if (this.electrodiff<=-0.4) {
	
			this.atom1.set_partCharge(0.5);
			this.atom2.set_partCharge(-0.5);
			
		}
		
		else if (this.electrodiff>=0.4) {
			
			this.atom1.set_partCharge(-0.5);
			this.atom2.set_partCharge(0.5);
			
		}
		
		else {
			
			this.atom1.set_partCharge(0);
			this.atom2.set_partCharge(0);
		}
		
	}
	
	//*******************************************************************	
	//sets gets**********************************************************

	public PVector get_dipole() {return this.dipole; }

	//*******************************************************************
	//behaviorals********************************************************

	public void image(Graphics g) {
		
		g.setColor(Color.DARK_GRAY);
		g.drawLine((int)(Bench.BEAKER_UNIT*Bench.BEAKER_X+Bench.BEAKER_UNIT*this.atom1.get_pos().get_x()), 
				(int)(Bench.BEAKER_UNIT*Bench.BEAKER_Y+Bench.BEAKER_UNIT*this.atom1.get_pos().get_y()), 
				(int)(Bench.BEAKER_UNIT*Bench.BEAKER_X+Bench.BEAKER_UNIT*this.atom2.get_pos().get_x()), 
				(int)(Bench.BEAKER_UNIT*Bench.BEAKER_Y+Bench.BEAKER_UNIT*this.atom2.get_pos().get_y()));
		
	}

	//*******************************************************************
}

class Synthesis extends Interaction {
	
	public static Compound NO2toN2O4(Compound molecule1, Compound molecule2){
		
		Random rand = new Random(); 
		double NOFixingDis = 1.21; 
		double NNFixingDis = 1.75; 
		double angleONO = 2.356; 
		double NNSeparatingProb = 0; 
		
		Compound thisN2O4 = new Compound(6); 
		
		PVector posN2O4 = PVector.multiply(PVector.add(molecule1.get_pos(), molecule2.get_pos()), 0.5); 
		
		Atom N1 = molecule1.get_components()[0]; 
		N1.set_pos(PVector.add(posN2O4, new PVector(-0.5*NNFixingDis, 0)));
		Atom N2 = molecule2.get_components()[0]; 
		N2.set_pos(PVector.add(posN2O4, new PVector(0.5*NNFixingDis, 0)));
		Atom O1 = molecule1.get_components()[1]; 
		O1.set_pos(PVector.add(posN2O4, new PVector(-0.5*NNFixingDis-Math.cos(angleONO/2)*NOFixingDis, Math.sin(angleONO/2)*NOFixingDis)));
		Atom O2 = molecule2.get_components()[1]; 
		O2.set_pos(PVector.add(posN2O4, new PVector(0.5*NNFixingDis+Math.cos(angleONO/2)*NOFixingDis, Math.sin(angleONO/2)*NOFixingDis)));
		Atom O3 = molecule1.get_components()[2]; 
		O3.set_pos(PVector.add(posN2O4, new PVector(-0.5*NNFixingDis-Math.cos(angleONO/2)*NOFixingDis, -Math.sin(angleONO/2)*NOFixingDis))); 
		Atom O4 = molecule2.get_components()[2]; 
		O4.set_pos(PVector.add(posN2O4, new PVector(0.5*NNFixingDis+Math.cos(angleONO/2)*NOFixingDis, -Math.sin(angleONO/2)*NOFixingDis)));
		
		Atom[] componentsN2O4 = {N1, N2, O1, O2, O3, O4}; 
		thisN2O4.set_components(componentsN2O4);
		thisN2O4.set_mass(molecule1.get_mass()+molecule2.get_mass());; 
		thisN2O4.set_pos(posN2O4);
		
		Intramolecular N_N = new Intramolecular(N1, N2, -1, -1, NNFixingDis, NNSeparatingProb); 
		Intramolecular N_O_1 = new Intramolecular(N1, O1, -1, -1, NOFixingDis, 0); 
		Intramolecular N_O_2 = new Intramolecular(N2, O2, -1, -1, NOFixingDis, 0); 
		Intramolecular N_O_3 = new Intramolecular(N1, O3, -1, -1, NOFixingDis, 0); 
		Intramolecular N_O_4 = new Intramolecular(N2, O4, -1, -1, NOFixingDis, 0); 
		
		ArrayList<Intramolecular> bondsN2O4 = new ArrayList<Intramolecular>(); 
		bondsN2O4.add(N_N); 
		bondsN2O4.add(N_O_1); 
		bondsN2O4.add(N_O_2); 
		bondsN2O4.add(N_O_3); 
		bondsN2O4.add(N_O_4);
		thisN2O4.set_bonds(bondsN2O4);
		thisN2O4.set_dipole(new PVector(0,0));
		
		return thisN2O4; 
		
	}
	
	public static void NO2toN2O4(Bench b, double bondingProb) {
		
		Random rand = new Random(); 
		
		double bondingDis = 1.85; 
		
		int i = 0; 
		
		while (i<b.listNO2.size()) {
			
			int j = i+1; 
			while(j<b.listNO2.size()) {
				
				Compound molecule1 = b.listNO2.get(i); 
				Compound molecule2 = b.listNO2.get(j); 
				
				PVector displacement = PVector.subtract(molecule1.get_pos(), molecule2.get_pos()); 
				double distance = displacement.get_mag(); 
				
				if(distance<=bondingDis && rand.nextDouble()<bondingProb) {
					
					b.listNO2.remove(molecule1); 
					b.listNO2.remove(molecule2); 
					
					b.listN2O4.add(NO2toN2O4(molecule1, molecule2)); 
				}
				
				j++; 
			}
			i++; 
		}
		
	}
	
}//end Synthesis

class SpontaneousDecomposition extends Interaction {
	
	public static Stack<Compound> N2O4toNO2(Compound N2O4) {
		
		double separatingDis = 2; 
		
		//find the location of the parent N2O4
		PVector center = N2O4.get_pos(); 
		
		//get the six atoms in N2O4
		Atom N_in_1 = N2O4.get_components()[0]; 
		Atom N_in_2 = N2O4.get_components()[1]; 
		Atom O1_in_1 = N2O4.get_components()[2]; 
		Atom O1_in_2 = N2O4.get_components()[3]; 
		Atom O2_in_1 = N2O4.get_components()[4]; 
		Atom O2_in_2 = N2O4.get_components()[5]; 
		
		//get the five bonds in N2O4
		Intramolecular N_O_1_in_1 = N2O4.get_bonds().get(1); 
		Intramolecular N_O_1_in_2 = N2O4.get_bonds().get(2);
		Intramolecular N_O_2_in_1 = N2O4.get_bonds().get(3); 
		Intramolecular N_O_2_in_2 = N2O4.get_bonds().get(4); 
		
		//put things that will end up in NO2_1 into an array/lists
		Atom[] components_in_1 = {N_in_1, O1_in_1, O2_in_1}; 
		ArrayList<Intramolecular> bonds_in_1 = new ArrayList<Intramolecular>(
				Arrays.asList(N_O_1_in_1, N_O_2_in_1)
				);
		
		//put things that will end up in NO2_2 into an array/list
		Atom[] components_in_2 = {N_in_2, O1_in_2, O2_in_2}; 
		ArrayList<Intramolecular> bonds_in_2 = new ArrayList<Intramolecular>(
				Arrays.asList(N_O_1_in_2, N_O_2_in_2)
				);
		
		//create the two NO2 molecules
		Compound NO2_1 = new Compound(components_in_1, bonds_in_1); 
		NO2_1.set_mass(46.0047);
		NO2_1.set_pos(N_in_1.get_pos());
		NO2_1.set_dipole(new PVector(0,0));
		
		Compound NO2_2 = new Compound(components_in_2, bonds_in_2); 
		NO2_2.set_mass(46.0047);
		NO2_2.set_pos(N_in_2.get_pos());
		NO2_2.set_dipole(new PVector(0,0));
		
		//justify the position of the two molecules
		NO2_1.move(PVector.subtract(center, NO2_1.get_pos()),0);
		NO2_2.move(PVector.subtract(center, NO2_2.get_pos()),0);
		//push the two molecules away to prevent rebond
		NO2_1.move(new PVector(-separatingDis/2, 0), 0);
		NO2_2.move(new PVector(separatingDis/2, 0), 0);
		
		
		//adjust the bond lengths by moving the O's
		PVector O1_in_1_arm = PVector.subtract(O1_in_1.get_pos(),NO2_1.get_pos()); 
		PVector O1_in_2_arm = PVector.subtract(O2_in_1.get_pos(), NO2_2.get_pos()); 
		PVector O2_in_1_arm = PVector.subtract(O2_in_1.get_pos(),NO2_1.get_pos()); 
		PVector O2_in_2_arm = PVector.subtract(O2_in_2.get_pos(),NO2_2.get_pos()); 
		
		O1_in_1.move(PVector.multiply(O1_in_1_arm, -0.013/1.21));
		O1_in_2.move(PVector.multiply(O1_in_2_arm, -0.013/1.21));
		O2_in_1.move(PVector.multiply(O2_in_1_arm, -0.013/1.21));
		O2_in_2.move(PVector.multiply(O2_in_2_arm, -0.013/1.21));
		
		
		//adjust the bond angles by moving the O's with rot-matrix
		//note we pick the O's that should be rotated anti-clockwise
		//in order to keep the signs of angles consistent
		double diffangle = -0.017; 
		O1_in_1.move(new PVector(
				O1_in_1_arm.get_x()*Math.cos(diffangle)-O1_in_1_arm.get_y()*Math.sin(diffangle)-O1_in_1_arm.get_x(),				
				O1_in_1_arm.get_x()*Math.sin(diffangle)+O1_in_1_arm.get_y()*Math.cos(diffangle)-O1_in_1_arm.get_y()
				));	
		O2_in_2.move(new PVector(
				O2_in_2_arm.get_x()*Math.cos(diffangle)-O2_in_2_arm.get_y()*Math.sin(diffangle)-O2_in_2_arm.get_x(),				
				O2_in_2_arm.get_x()*Math.sin(diffangle)+O2_in_2_arm.get_y()*Math.cos(diffangle)-O2_in_2_arm.get_y()
				));	
		
		//put the two NO2 molecules into an array and return it
		Stack<Compound> NO2s = new Stack<Compound>();
		NO2s.push(NO2_1); 
		NO2s.push(NO2_2); 
		return NO2s; 
		
	}
	
	public static void N2O4toNO2(Bench b, double separatingProb) {
		
		Random rand = new Random(); 
		
		for (int i = 0; i<b.listN2O4.size(); i++) {
			
			Compound N2O4 = b.listN2O4.get(i); 
			double seed = rand.nextDouble(); 
			if (seed<separatingProb) {
				
				Stack<Compound> newNO2s = N2O4toNO2(N2O4); 
				b.listNO2.add(newNO2s.pop()); 
				b.listNO2.add(newNO2s.pop()); 
				b.listN2O4.remove(i); 
			}
		}
		
		return; 		
	}
	
}