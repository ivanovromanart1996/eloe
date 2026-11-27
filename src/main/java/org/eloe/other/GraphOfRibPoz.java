package org.eloe.other;

import java.awt.Frame;
import java.awt.Graphics;

class GraphOfRibPoz extends Frame {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	
	private int genesNumber;
	private int[] a;
	private int[] b;
	private int[] c;
	private int[] d;
	private int[] e;
	
	public GraphOfRibPoz(final String s, final int n, final int[] aa, final int[] bb, final int[] cc, final int[] dd, final int[] ee) {
		super(s);
		setBounds(0, 0, 800, 230);
				
		this.genesNumber = n;
		this.a = aa;
		this.b = bb;
		this.c = cc;
		this.d = dd;
		this.e = ee;
		
		setVisible(true);
	}
	
	public void paint(final Graphics g) {
		double k = 770.0 / this.genesNumber;
		int len = this.a.length;
		
		g.drawLine(10 + (int) (1 * k), 47, 10 + (int) (this.genesNumber * k), 47);
		for (int i = 0; i < len; i++) {
			g.drawOval(10 + (int) (this.a[i] * k), 45, 4, 4);
		}
		
		g.drawLine(1, 64, getWidth(), 64);
		g.drawLine(1, 66, getWidth(), 66);
		len = this.b.length;
		
		g.drawLine(10 + (int) (1 * k), 87, 10 + (int) (this.genesNumber * k), 87);
		for (int i = 0; i < len; i++) {
			g.drawOval(10 + (int) (this.b[i] * k), 85, 4, 4);
		}
		
		g.drawLine(1, 104, getWidth(), 104);
		g.drawLine(1, 106, getWidth(), 106);
		len = this.c.length;
		
		g.drawLine(10 + (int) (1 * k), 127, 10 + (int) (this.genesNumber * k), 127);
		for (int i = 0; i < len; i++) {
			g.drawOval(10 + (int) (this.c[i] * k), 125, 4, 4);
		}
		
		g.drawLine(1, 144, getWidth(), 144);
		g.drawLine(1, 146, getWidth(), 146);
		len = this.d.length;
		
		g.drawLine(10 + (int) (1 * k), 167, 10 + (int) (this.genesNumber * k), 167);
		for (int i = 0; i < len; i++) {
			g.drawOval(10 + (int) (this.d[i] * k), 165, 4, 4);
		}
		
		g.drawLine(1, 184, getWidth(), 184);
		g.drawLine(1, 186, getWidth(), 186);
		len = this.e.length;
		
		g.drawLine(10 + (int) (1 * k), 207, 10 + (int) (this.genesNumber * k), 207);
		for (int i = 0; i < len; i++) {
			g.drawOval(10 + (int) (this.e[i] * k), 205, 4, 4);
		}
	}
}