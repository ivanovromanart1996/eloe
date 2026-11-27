package org.eloe.gui.task;

import java.io.File;

import javax.swing.JButton;
import javax.swing.JProgressBar;
import javax.swing.JTextArea;

import org.apache.log4j.Logger;

import org.eloe.engine.EEIEngine;
import org.eloe.parameters.ParametersStructure;

public class Task implements Runnable {
	
	private static final Logger LOG = Logger.getLogger(Task.class);
	
	private ParametersStructure ps;
	private File sourceFolder;
	private File resultFolder;
	private JButton btnStart;
	private JButton btnParameters;
	private JProgressBar progressBar;
	private JTextArea console;
	
	public Task(ParametersStructure ps1, File sourceFolder, File resultFolder, JButton btnStart1, JButton btnParameters1, JProgressBar progressBar1, JTextArea console1) {
		this.ps = ps1;
		this.sourceFolder = sourceFolder;
		this.resultFolder = resultFolder;
		this.btnStart = btnStart1;
		this.btnParameters = btnParameters1;
		this.progressBar = progressBar1;
		this.console = console1;
	}
	
	@Override
	public void run() {
		try {
			EEIEngine.engine(ps, sourceFolder, resultFolder, console);
			btnStart.setEnabled(true);
			btnParameters.setEnabled(true);
			progressBar.setIndeterminate(false);
		} catch (NumberFormatException e) {
			e.printStackTrace();
			LOG.error("", e);
		}
	}
	
}
