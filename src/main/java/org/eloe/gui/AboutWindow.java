package org.eloe.gui;

import java.awt.BorderLayout;
import java.awt.FlowLayout;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JPanel;
import javax.swing.border.EmptyBorder;
import javax.swing.JLabel;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;

public class AboutWindow extends JDialog {

	private static final long serialVersionUID = 1L;
	private final JPanel contentPanel = new JPanel();
	
	/**
	 * Create the dialog.
	 */
	public AboutWindow() {
		setResizable(false);
		initialize();
	}
	
	private void initialize() {
		setTitle("About");
		setBounds(300, 300, 340, 160);
		getContentPane().setLayout(new BorderLayout());
		contentPanel.setBorder(new EmptyBorder(5, 5, 5, 5));
		getContentPane().add(contentPanel, BorderLayout.CENTER);
		contentPanel.setLayout(null);
		{
			JLabel lblThisProgrammIs = new JLabel("The programm was created for calculation of");
			lblThisProgrammIs.setBounds(10, 11, 304, 14);
			contentPanel.add(lblThisProgrammIs);
		}
		{
			JLabel lblNewLabel = new JLabel("elongation efficiency index (EEI).");
			lblNewLabel.setBounds(10, 36, 304, 14);
			contentPanel.add(lblNewLabel);
		}
		{
			JLabel lblWarlord = new JLabel("\u00A9 Warlord, V. S. Sokolov, 1.02.14");
			lblWarlord.setBounds(10, 64, 304, 14);
			contentPanel.add(lblWarlord);
		}
		{
			JPanel buttonPane = new JPanel();
			getContentPane().add(buttonPane, BorderLayout.SOUTH);
			buttonPane.setLayout(new FlowLayout(FlowLayout.CENTER, 5, 5));
			{
				JButton okButton = new JButton("OK");
				okButton.addActionListener(new ActionListener() {
					public void actionPerformed(ActionEvent e) {
						AboutWindow.this.dispose();
					}
				});
				okButton.setActionCommand("OK");
				buttonPane.add(okButton);
				getRootPane().setDefaultButton(okButton);
			}
		}
	}
	

}
