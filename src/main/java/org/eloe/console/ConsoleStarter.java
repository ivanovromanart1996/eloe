package org.eloe.console;

import java.awt.EventQueue;
import java.io.File;
import org.apache.log4j.Logger;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import org.eloe.parameters.ParametersStructure;

/**
 * @author yasana
 * Пусковик из консоли
 */
public class ConsoleStarter {

	private static final Logger LOG = Logger.getLogger(ConsoleStarter.class);
	
	private ParametersStructure ps;
	private File sourceFolder;
	private File resultFolder;
	private Thread treadForTask = null;
	
	private void initialize() {
		ps = new ParametersStructure();
	}
	
	public ConsoleStarter() {
		initialize();
	}

	public void start() {
		if ((sourceFolder != null) && (resultFolder != null)) {
			LOG.info("\"Start\" method was invoked.");
			ConsoleTask newtask = new ConsoleTask(ps, sourceFolder, resultFolder);
			treadForTask = new Thread(newtask, "NewTask");
			treadForTask.setPriority(Thread.NORM_PRIORITY);
			treadForTask.start();
		}
		else {
			LOG.info("At first you must choose \"source\" and \"result\" directories!");
		}
	}
	
	//args[0] - use
	//args[1] - res
	public static void main(String[] args) {
		LOG.info("Application start.");
		String use = "use", res = "res";
		Path currentDirectory = Paths.get("").toAbsolutePath();//Get Current Path
		if (args.length>0) {
			if (args.length == 1) {
				LOG.info("Wrong number of arguments! Two arguments are needed: source folder and result folder");
			}
			use = args[0];
			res = args[1];
		}
		//Строим абсолютные пути, находя нужные папки
		final Path usePath = currentDirectory.resolve(use);
		final Path resPath = currentDirectory.resolve(res);
		if (Files.notExists(usePath)) {
			LOG.info("Error! There is no source folder!");
		}
		if (Files.notExists(resPath)) {
			LOG.info("Error! There is no result folder!");
		}
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					ConsoleStarter starter = new ConsoleStarter();
					starter.sourceFolder = usePath.toFile();
					starter.resultFolder = resPath.toFile();
					starter.start();
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		});
	}
}
