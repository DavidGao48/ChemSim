import java.util.*; 
import java.math.*; 

import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.*;

//**************************************************************
//**************************************************************
//**************************************************************
//**************************************************************
public abstract class Material {
	
	//fields and defaults***********************************************
	
	private double mass; 
	private double volume; 
	private double temperature;
	private Molecule unit; //the repeated molecule in this material
	
	public Material(){
		
		this.mass = -1; //using -1 to emphasize unknownness 
		this.volume = -1;
		this.temperature = -1; 
		this.unit = null; 
		
	}
	
	public Material(double mass, double volume, double initialTemperature, Molecule unit) {
		
		this.mass = mass; 
		this.volume = volume; 
		this.temperature = initialTemperature; 
		this.unit = unit; 
		
	}
	
	//*******************************************************************	
	//sets gets**********************************************************
	
		//Defensive code needed!!! 
	public void set_mass(double mass) {this.mass = mass; }
	public void set_volume(double volume) {this.volume = volume; }
	public void set_temperature(double temperature) {this.temperature = temperature; }
	public void set_unit(Molecule unit) {this.unit = unit; }
	
	public double get_mass() {return this.mass; }
	public double get_volume() {return this.volume; }
	public double get_temperature() {return this.temperature; }
	public double get_density() {return this.mass/this.volume; }
	public Molecule get_unit() {return this.unit; }
	
	//*******************************************************************
	//behaviorals********************************************************
	
	
	
	//*******************************************************************
}//end Material

//**************************************************************
//**************************************************************
//**************************************************************
//**************************************************************
class Solid extends Material {
	
	//fields and defaults***********************************************
	
	private double meltingPoint; //Kelvin
	private String color; 
	private String crystalType; //ionic, pmolecular, npmolecular, metallic, network
	private boolean maleability; //consider switching to double if beyond AP level
	
	public Solid() {
		
		super(); 
		
		this.meltingPoint = -1; 
		this.color = null; 
		this.crystalType = null; 
		this.maleability = false; //maleability assumed false
		
	}
	
	public Solid(Material source) {
		
		super(source.get_mass(), source.get_volume(), source.get_temperature(), source.get_unit()); 
		
		this.meltingPoint = -1; 
		this.color = null; 
		this.crystalType = null; 
		this.maleability = false; //maleability assumed false
		
	}
	
	public Solid(double mass, double volume, double initialTemperature, Molecule unit, double meltingPoint, String color, String crystalType, boolean maleability) {
		
		super(mass, volume, initialTemperature, unit); 
		
		this.meltingPoint = meltingPoint; 
		this.color = color; 
		this.crystalType = crystalType; 
		this.maleability = maleability; 
		
	}
	
	
	//*******************************************************************	
	//sets gets**********************************************************
	
		//Defensive code needed
	public void set_meltingPoint(double meltingPoint) {this.meltingPoint = meltingPoint; }
	public void set_color(String color) {this.color = color; }
	public void set_crystalType(String crystalType) {this.crystalType = crystalType; }
	public void set_maleability(boolean maleability) {this.maleability = maleability; }
	
	public double get_meltingPoint() {return this.meltingPoint; }
	public String get_color() {return this.color; }
	public String get_crystalType() {return this.crystalType; }
	public boolean get_maleability() {return this.maleability; }
	
	//*******************************************************************
	//behaviorals********************************************************
	
	
	
	//*******************************************************************
}//end Solid

//**************************************************************
//**************************************************************
//**************************************************************

class SaltCrystal extends Solid {
	
	//fields and defaults***********************************************
	
	private double latticeEnergy; //the work needed to separate ions
	private Ion cation; //the positive ion
	private Ion anion; //the negative ion
	
	public SaltCrystal() {
		
		super(); 
		this.set_crystalType("ionic"); //all salts are ionic crystals
		
		this.latticeEnergy = -1; 
		this.cation = null; 
		this.anion = null; 
		
	}
	
	public SaltCrystal(double latticeEnergy, Ion cation, Ion anion) {
		
		super(); 
		this.set_crystalType("ionic"); //all salts are ionic crystals

		this.latticeEnergy = latticeEnergy; 
		this.cation = cation; 
		this.anion = anion; 
		
	}
	
	//*******************************************************************	
	//sets gets**********************************************************

		//Defensive code needed
	public void set_latticeEnergy(double latticeEnergy) {this.latticeEnergy = latticeEnergy; }
	public void set_cation(Ion cation) {this.cation = cation; }
	public void set_anion(Ion anion) {this.anion = anion; }
	
	public double get_latticeEnergy() {return this.latticeEnergy; }
	public Ion get_cation() {return this.cation; }
	public Ion get_anion() {return this.anion; }

	//*******************************************************************
	//behaviorals********************************************************



	//*******************************************************************
}

//**************************************************************
//**************************************************************
//**************************************************************
//**************************************************************

class Liquid extends Material {
	
	//fields and defaults***********************************************
	
	private Intermolecular intermolecular; //the defining interaction in this liquid
	private double freezingPoint; //in Kelvin
	private double evaporationPoint; //in Kelvin
	private String color; 
	
	public Liquid() {
		
		super(); 
		
		this.intermolecular = null; 
		this.freezingPoint = -1; 
		this.evaporationPoint = -1; 
		this.color = null; 
		
	}
	
	public Liquid(Material source) {
		
		super(source.get_mass(), source.get_volume(), source.get_temperature(), source.get_unit()); 
		
		this.intermolecular = null; 
		this.freezingPoint = -1; 
		this.evaporationPoint = -1; 
		this.color = null; 
		
	}
	
	public Liquid(double mass, double volume, double initialTemperature, Molecule unit, Intermolecular intermolecular, double freezingPoint, double evaporationPoint, String color) {
		
		super(mass, volume, initialTemperature, unit); 
		
		this.intermolecular = intermolecular; 
		this.freezingPoint = freezingPoint; 
		this.evaporationPoint = evaporationPoint; 
		this.color = color; 
		
	}
	
	//*******************************************************************	
	//sets gets**********************************************************

		//Defensive code needed
	public void set_intermolecular(Intermolecular intermolecular) {this.intermolecular = intermolecular; }
	public void set_freezingPoint(double freezingPoint) {this.freezingPoint = freezingPoint; }
	public void set_evaporationPoint(double evaporationPoint) {this.evaporationPoint = evaporationPoint; }
	public void set_color(String color) {this.color = color; }
	
	public Intermolecular get_intermolecular() {return this.intermolecular; }
	public double get_freezingPoint() {return this.freezingPoint; }
	public double get_evaporationPoint() {return this.evaporationPoint; }
	public String get_color() {return this.color; }

	//*******************************************************************
	//behaviorals********************************************************

	

	//*******************************************************************
}

//**************************************************************
//**************************************************************
//**************************************************************

class LSolvent extends Liquid {
	
	//fields and defaults***********************************************
	
	
	
	//*******************************************************************	
	//sets gets**********************************************************



	//*******************************************************************
	//behaviorals********************************************************



	//*******************************************************************
}

abstract class Gas extends Material {
	//abstract because a gas is just an arraylist of molecules. 
	//therefore the methods in this class only describe how 
	//gas molecules might interact
	
	//fields and defaults***********************************************
	
	
	
	//*******************************************************************	
	//sets gets**********************************************************



	//*******************************************************************
	//behaviorals********************************************************

	public static void Brownian(ArrayList<Compound> gas) {
		
		Random rand = new Random(); 
		
		for (int i = 0; i<gas.size(); i++) {
			
			Compound thismolecule = gas.get(i); 
			
			PVector dpos = new PVector(rand.nextDouble()*2-1, rand.nextDouble()*2-1); 
			
			if (thismolecule.get_x()<=2) {dpos.set_x(1);}
			else if (thismolecule.get_x()>=Bench.BEAKER_WIDTH-2) {dpos.set_x(-1);}
			if(thismolecule.get_y()<=2) {dpos.set_y(1);}
			else if (thismolecule.get_y()>=Bench.BEAKER_DEPTH-2) {dpos.set_y(-1);}
			
			double dtheta = -1+2*rand.nextDouble(); 
			
			thismolecule.move(dpos, dtheta);
			
		}
		
	}
	

	//*******************************************************************

}

class Molecule extends Material {
	
	//fields and defaults***********************************************
	
	private double radius; 
	private PVector pos; 
	private PVector dipole;
	
	public Molecule() {
		
		super(); 
		
		this.radius = -1; 
		this.pos = new PVector(); 
		this.dipole = new PVector(); 
		
	}
	
	public Molecule(double mass, double radius, PVector pos, PVector dipole) {
		
		super(mass, -1, -1, null); 
		
		this.radius = radius; 
		this.pos = pos; 
		this.dipole = dipole; 
		
	}
	
	//*******************************************************************	
	//sets gets**********************************************************
	
	public void set_radius(double radius) {this.radius = radius; }
	public void set_pos(PVector pos) {this.pos = pos; }
	public void set_x(double x) {this.pos.set_x(x); }
	public void set_y(double y) {this.pos.set_y(y); }
	public void set_dipole(PVector dipole) {this.dipole = dipole; }
	public void set_dipole(double mag, double theta){
		this.dipole.set_x(mag*Math.cos(theta)); 
		this.dipole.set_y(mag*Math.sin(theta));
		} 
	
	public double get_radius() {return this.radius; }
	public PVector get_pos() {return this.pos; }
	public double get_x() {return this.pos.get_x(); }
	public double get_y() {return this.pos.get_y(); }
	public PVector get_dipole() {return this.dipole; }
	public double get_dipoleMag() {return this.dipole.get_mag(); }
	public double get_dipoleTheta() {return this.dipole.get_theta(); }
	
	//*******************************************************************
	//behaviorals********************************************************

	

	//*******************************************************************
}

class Atom extends Molecule {
	
	//fields and defaults***********************************************
	
	private int number; 
	private int group; //we will use group=0 to denote transition elements. 
	private double partCharge; 
	private double electronegativity; 
	private int countLonePair; 
	private Compound in; 
	
	private Color color = Color.BLACK; 
	
	public Atom() {
		
		super(); 
		this.set_dipole(new PVector(0,0));//atoms have 0 dipoles 
		
		this.number = -1; 
		this.group=-1; 
		this.partCharge = 0; 
		this.electronegativity = 0; 
		this.countLonePair = -1; 
		this.in = null; 
		
	}
	
	public Atom(double mass, double radius, PVector pos, int number, int group, double partCharge, double electronegativity, int countLonePair, Compound in) {
		
		super();
		this.set_dipole(new PVector(0,0));//atoms have 0 dipoles
		this.set_mass(mass); 
		this.set_radius(radius);
		this.set_pos(pos); 
		
		this.number = number; 
		this.group = group; 
		this.partCharge = partCharge; 
		this.electronegativity = electronegativity; 
		this.countLonePair = countLonePair; 
		this.in = in; 
		
	}
	
	//*******************************************************************	
	//sets gets**********************************************************

	public void set_number(int number) {this.number = number; }
	public void set_group(int group) {this.group = group; }
	public void set_partCharge(double partCharge) {this.partCharge = partCharge; }
	public void set_electronegativity(double electronegativity) {this.electronegativity = electronegativity; }
	public void set_countLonePair(int countLonePair) {this.countLonePair = countLonePair; }
	public void set_inCompound(Compound in) {this.in = in; }
	
	public int get_number() {return this.number; }
	public int get_group() {return this.group; }
	public double get_partCharge() {return this.partCharge; }
	public double get_electronegativity() {return this.electronegativity; }
	public int get_countLonePair() {return this.countLonePair; }
	public Compound get_inCompound() {return this.in; }
	
	//*******************************************************************
	//behaviorals********************************************************

	public void image(Graphics g) {
		
		g.setColor(this.color);
		g.fillOval(
				(int)(Bench.BEAKER_UNIT*(Bench.BEAKER_X+this.get_pos().get_x()-this.get_radius())),
				(int)((Bench.BEAKER_UNIT*(Bench.BEAKER_Y+this.get_pos().get_y()-this.get_radius()))),
				(int) (Bench.BEAKER_UNIT*2*this.get_radius()),
				(int) (Bench.BEAKER_UNIT*2*this.get_radius())
				); 
		
	}

	public void move(PVector dpos) {
		
		this.get_pos().plus(dpos);;
		
	}

	//*******************************************************************
}

class Ion extends Molecule {
	
	//fields and defaults***********************************************
	

	private Atom sourceAtom; 
	private boolean isPolyatomic; 
	private int charge; 
	
	public Ion() {
		
		super(); 
		
		this.sourceAtom=null; 
		this.isPolyatomic = false; 
		this.charge=0; 
		
	}
	
	public Ion(Atom sourceAtom, int charge) {
		
		super(sourceAtom.get_mass(), sourceAtom.get_radius(), sourceAtom.get_pos(), new PVector(0,0)); 
		
		this.sourceAtom=sourceAtom; 
		this.charge = charge; 
	}
	
	//*******************************************************************	
	//sets gets**********************************************************

	public void set_sourceAtom(Atom sourceAtom) {this.sourceAtom = sourceAtom; }
	public void set_charge(int charge) {this.charge = charge; }
	public void set_polyatomic(boolean isPolyatomic) {this.isPolyatomic = isPolyatomic; }

	public Atom get_sourceAtom() {return this.sourceAtom; }
	public int get_charge() {return this.charge; }
	public boolean is_polyatomic() {return this.isPolyatomic; }
	
	//*******************************************************************
	//behaviorals********************************************************



	//*******************************************************************
}

class Compound extends Molecule {
	
	//fields and defaults***********************************************
	
	private int size; 
	private Atom[] components; 
	private ArrayList<Intramolecular> bonds; 
	
	public Compound (int size) {
		
		super(); 
		
		this.size = size; 
		this.components = new Atom[size]; 
		this.bonds = new ArrayList<Intramolecular>(); 
		
	}
	
	public Compound (Atom[] components, ArrayList<Intramolecular> bonds) {
		
		super(); 
		
		this.size = components.length; 
		this.set_components(components);
		this.set_bonds(bonds);
		
	}
	
	//*******************************************************************	
	//sets gets**********************************************************

	public void set_components(Atom[] components) {
		
		this.components = components; 
		this.size = components.length; 
		
		double gMM = 0; 
		for (int i = 0; i<components.length; i++) {gMM+=components[i].get_mass(); }
		this.set_mass(gMM);
	}
	
	public void set_bonds(ArrayList<Intramolecular> bonds) {
		
		this.bonds = bonds; 
		
		PVector dipole = new PVector(); 
		//for (int i = 0; i<bonds.size(); i++) {
			//dipole.plus(bonds.get(i).get_dipole()); 
			//}
		this.set_dipole(dipole);
		
	}
	
	public int get_size() {return this.size; }
	public Atom[] get_components() {return this.components; }
	public ArrayList<Intramolecular> get_bonds() {return this.bonds; }

	//*******************************************************************
	//behaviorals********************************************************

	public void image(Graphics g) {
		
		g.setColor(Color.BLUE);
		for (int i = 0; i<this.components.length; i++) {
			this.components[i].image(g);
		}
		for (int i = 0; i<this.bonds.size(); i++) {
			this.bonds.get(i).image(g);
		}
	}
	
	public void move(PVector dpos, double rotation) {
		
		this.set_pos(PVector.add(this.get_pos(),dpos));
		
		for (int i = 0; i<this.components.length; i++) {
			
			Atom thisAtom = this.components[i]; 
			thisAtom.move(dpos);
			
			
			PVector thisAtomArm = PVector.subtract(thisAtom.get_pos(),this.get_pos());
			
			double armX = thisAtomArm.get_x(); 
			double armY = thisAtomArm.get_y(); 
		
			thisAtom.move(new PVector(armX*Math.cos(rotation)-armX-armY*Math.sin(rotation),armX*Math.sin(rotation)-armY+armY*Math.cos(rotation)));
			
		}
	
		
	}

	//*******************************************************************
}



