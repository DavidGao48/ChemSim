import java.awt.Color;
import java.awt.Dimension;

import javax.swing.*;

public class Lab {

	private JFrame room; 
	private JPanel bench; 
	
	//Construction of the lab
	
	public Lab(String title) {
		
		this.room = new JFrame(); 
		this.room.setVisible(true);
		this.room.setSize(2000,1000);
		this.room.setBackground(Color.WHITE);
		this.room.setLocationRelativeTo(null);
		this.room.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		this.room.setTitle(title);
		
		this.room.getContentPane().add(new Bench());
	}
	
}
