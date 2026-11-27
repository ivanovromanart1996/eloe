package org.eloe.gui;

import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JMenuBar;
import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.JButton;
import javax.swing.JOptionPane;
import javax.swing.JLabel;
import javax.swing.JTextField;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.JProgressBar;

import org.apache.log4j.Logger;

import java.awt.EventQueue;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;

import org.eloe.gui.task.Task;
import org.eloe.parameters.ParametersStructure;

//
public class MainWindow {
	
	private static final Logger LOG = Logger.getLogger(MainWindow.class);
	
	private JFrame frmEeicalculator;
	private JButton btnStart;
	private JButton btnStop;
	private JButton btnParameters;
	private AboutWindow dialogAbout;
	private ParametersWindow dialogParameters;
	private JTextField use;
	private ParametersStructure ps;
	private File sourceFolder;
	private File resultFolder;
	private JTextArea console;
	private JTextField result;
	private JProgressBar progressBar;
	private Thread treadForTask = null;
	
	public MainWindow() {
		initialize();
	}

	public static void main(String[] args) {
		LOG.info("Application start.");
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					MainWindow window = new MainWindow();
					window.frmEeicalculator.setVisible(true);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		});
	}

	private void initialize() {
		ps = new ParametersStructure();
		
		frmEeicalculator = new JFrame();
		frmEeicalculator.addWindowListener(new WindowAdapter() {
			@Override
			public void windowClosing(WindowEvent arg0) {
				treadForTask = null;
				frmEeicalculator.dispose();
				LOG.info("Application stop.");
			}
		});
		frmEeicalculator.setResizable(false);
		frmEeicalculator.setTitle("EEI-Calculator");
		frmEeicalculator.setBounds(100, 100, 450, 450);
		
		JMenuBar menuBar = new JMenuBar();
		frmEeicalculator.setJMenuBar(menuBar);
		
		JMenu mnFile = new JMenu("File");
		menuBar.add(mnFile);
		
		JMenuItem mntmStop = new JMenuItem("Stop");
		mntmStop.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				LOG.error("\"Stop\" button was pushed!");
				if (treadForTask != null && treadForTask.isAlive()) {
					progressBar.setIndeterminate(false);
					treadForTask.stop();
					treadForTask= null;
					JOptionPane.showMessageDialog(btnStart.getParent(), "Warning! The thread was stopped."
							+ "\nYou have to restart the application for its correct work.", "Warning", JOptionPane.WARNING_MESSAGE);
					frmEeicalculator.dispose();
					LOG.info("Application stop.");
				}
				btnStart.setEnabled(true);
				btnParameters.setEnabled(true);
			}
		});
		mnFile.add(mntmStop);
		
		JMenuItem mntmExit = new JMenuItem("Exit");
		mntmExit.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				treadForTask = null;
				frmEeicalculator.dispose();
				LOG.info("Application stop.");
			}
		});
		mnFile.add(mntmExit);
		
		JMenu mnAbout = new JMenu("About");
		menuBar.add(mnAbout);
		
		JMenuItem mntmAboutEeicalculator = new JMenuItem("About EEI-Calculator");
		mntmAboutEeicalculator.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				dialogAbout = new AboutWindow();
				dialogAbout.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
				dialogAbout.setModal(true);
				dialogAbout.setVisible(true);
			}
		});
		mnAbout.add(mntmAboutEeicalculator);
		
		frmEeicalculator.getContentPane().setLayout(null);
		
		btnStart = new JButton("Start");//=============================================================================<Start
		btnStart.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				if ((sourceFolder != null) && (resultFolder != null)) {
					LOG.info("\"Start\" button was pushed.");
					btnStart.setEnabled(false);
					btnParameters.setEnabled(false);
					//btnStop.setEnabled(true);
					console.setText(null);
					progressBar.setIndeterminate(true);
					Task newtask = new Task(ps, sourceFolder, resultFolder, btnStart, btnParameters, progressBar, console);
					treadForTask = new Thread(newtask, "NewTask");
					treadForTask.setPriority(Thread.NORM_PRIORITY);
					treadForTask.start();
				}
				else {
					JOptionPane.showMessageDialog(btnStart.getParent(), "At first you must choose \"source\" and \"result\" directories!", 
							"Warning", JOptionPane.INFORMATION_MESSAGE);
				}
			}
		});
		btnStart.setBounds(152, 90, 140, 46);
		frmEeicalculator.getContentPane().add(btnStart);
		
		btnStop = new JButton("Stop");
		btnStop.setEnabled(false);
		btnStop.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				treadForTask= null;
				btnStop.setEnabled(false);
				btnStart.setEnabled(true);
				btnParameters.setEnabled(true);
			}
		});
		btnStop.setBounds(299, 90, 140, 46);
		frmEeicalculator.getContentPane().add(btnStop);
		
		btnParameters = new JButton("Parameters");
		btnParameters.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				try {
					dialogParameters = new ParametersWindow(ps);
					dialogParameters.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
					dialogParameters.setModal(true);
					dialogParameters.setVisible(true);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		});
		btnParameters.setBounds(5, 90, 140, 46);
		frmEeicalculator.getContentPane().add(btnParameters);
		
		final JLabel lblDirectorySources = new JLabel("Choose directory with organism/s:");
		lblDirectorySources.setBounds(5, 5, 394, 14);
		frmEeicalculator.getContentPane().add(lblDirectorySources);
		
		use = new JTextField();
		use.setEditable(false);
		use.setBounds(5, 20, 405, 20);
		use.setColumns(10);
		frmEeicalculator.getContentPane().add(use);
		
		JButton btnSources = new JButton("...");
		btnSources.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				JFileChooser fc = new JFileChooser(System.getProperty("user.dir"));
				fc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
				int returnVal = fc.showOpenDialog(frmEeicalculator);
				if (returnVal == JFileChooser.APPROVE_OPTION) {
					sourceFolder = fc.getSelectedFile();
					use.setText(sourceFolder.toString());
				}
			}
		});
		btnSources.setBounds(416, 19, 23, 23);
		frmEeicalculator.getContentPane().add(btnSources);
		
		JLabel lblDirectoryResult = new JLabel("Choose directory for results:");
		lblDirectoryResult.setBounds(5, 45, 394, 14);
		frmEeicalculator.getContentPane().add(lblDirectoryResult);
		
		result = new JTextField();
		result.setEditable(false);
		result.setBounds(5, 60, 405, 20);
		result.setColumns(10);
		frmEeicalculator.getContentPane().add(result);
		
		JButton btnResult = new JButton("...");
		btnResult.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				JFileChooser fc = new JFileChooser(System.getProperty("user.dir"));
				fc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
				int returnVal = fc.showOpenDialog(frmEeicalculator);
				if (returnVal == JFileChooser.APPROVE_OPTION) {
					resultFolder = fc.getSelectedFile();
					result.setText(resultFolder.toString());
				}
			}
		});
		btnResult.setBounds(416, 59, 23, 23);
		frmEeicalculator.getContentPane().add(btnResult);
		
		JScrollPane scrollPane = new JScrollPane();
		scrollPane.setBounds(5, 145, 434, 230);
		frmEeicalculator.getContentPane().add(scrollPane);
		
		console = new JTextArea();
		console.setEditable(false);
		scrollPane.setViewportView(console);
		
		progressBar = new JProgressBar();
		progressBar.setBounds(5, 380, 434, 14);
		progressBar.setIndeterminate(false);
		frmEeicalculator.getContentPane().add(progressBar);
	
	}
}
