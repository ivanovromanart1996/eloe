package org.eloe.console;

import java.io.File;

import javax.swing.JTextArea;

import org.apache.log4j.Logger;

import org.eloe.engine.EEIEngine;
import org.eloe.parameters.ParametersStructure;

public class ConsoleTask implements Runnable {
	private static final Logger LOG = Logger.getLogger(ConsoleTask.class);
	
	private ParametersStructure ps;
	private File sourceFolder;
	private File resultFolder;

	public ConsoleTask(ParametersStructure ps, File sourceFolder, File resultFolder) {
		super();
		this.ps = ps;
		this.sourceFolder = sourceFolder;
		this.resultFolder = resultFolder;
	}

	@Override
	public void run() {
		try {
			JTextArea consoleStub = new JTextArea();//заглушка, которая по-большому счёту нам и не нужна
			//UPD: проверила код. Не должно быть проблем от использования этой заглушки
			EEIEngine.engine(ps, sourceFolder, resultFolder, consoleStub);
		} catch (NumberFormatException e) {
			e.printStackTrace();
			LOG.error("", e);
		}
	}

}
