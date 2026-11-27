package org.eloe.gui;

import java.awt.BorderLayout;
import java.awt.FlowLayout;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JPanel;
import javax.swing.border.EmptyBorder;
import javax.swing.JOptionPane;
import javax.swing.JTabbedPane;
import javax.swing.JLabel;
import javax.swing.JTextField;
import javax.swing.JCheckBox;
import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;

import org.eloe.parameters.ParametersStructure;

import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import javax.swing.JRadioButton;

public class ParametersWindow extends JDialog {
	
	private final static Pattern CODONPATTERN = Pattern.compile("([atgcATGC]{3}\\s*;\\s*)*([atgcATGC]{3}\\s*;?\\s*){1}");

	private static final long serialVersionUID = 1L;
	private final JPanel contentPanel = new JPanel();
	private JTextField lFlang;
	private JTextField rFlang;
	private JTextField distBOper;
	private JTextField minLenGen;
	private JTextField startCodons;
	private JTextField midlCodons;
	private JTextField stopCodons;
	private JCheckBox checkStart;
	private JCheckBox deleteGwBC;
	private JCheckBox checkStop;
	private JCheckBox usePseudoGenes;
	private JCheckBox autocounting;
	private JSpinner maxNumbCod;
	private JSpinner unUseCodonsL1;
	private JSpinner unUseCodonsR1;
	private JSpinner smin1;
	private JSpinner smax1;
	private JSpinner lmin1;
	private JSpinner lmax1;
	private JSpinner smin2;
	private JSpinner smax2;
	private JSpinner lmin2;
	private JSpinner lmax2;
	private JTextField minloopEnergy;
	private JCheckBox makeFileFreq;
	private JCheckBox makeFileWholeGenome;
	private JCheckBox searchOrthologs;
	private JSpinner number_orthologs;
	private JCheckBox countLCIij;
	private JCheckBox countRepeatsInGenes;
	private JCheckBox showRibPos;
	private JRadioButton rdbtnNumberOfGenes;
	private JRadioButton rdbtnPercentageOfAll;
	private JTextField NumbOfGenes;
	private JSpinner PercentageOfGenes;
	private ParametersStructure psNew;
	private ParametersStructure psDefault;
	
	
	public ParametersWindow(ParametersStructure ps1) {
		setResizable(false);
		psDefault = new ParametersStructure();
		psNew = ps1;
		setTitle("Parameters");
		setBounds(200, 200, 450, 320);
		getContentPane().setLayout(new BorderLayout());
		contentPanel.setBorder(new EmptyBorder(5, 5, 5, 5));
		getContentPane().add(contentPanel, BorderLayout.CENTER);
		contentPanel.setLayout(null);
		{
			JTabbedPane tabbedPane = new JTabbedPane(JTabbedPane.TOP);
			tabbedPane.setBounds(0, 0, 444, 259);
			tabbedPane.setTabLayoutPolicy(JTabbedPane.SCROLL_TAB_LAYOUT);
			contentPanel.add(tabbedPane);
			//Extraction tab
			{
				JPanel panelExtraction = new JPanel();
				tabbedPane.addTab("Extraction", null, panelExtraction, null);
				panelExtraction.setLayout(null);
				{
					JLabel lblLenghtOfLeft = new JLabel("Lenght of left flang:");
					lblLenghtOfLeft.setBounds(5, 5, 328, 14);
					panelExtraction.add(lblLenghtOfLeft);
				}
				{
					JLabel lblLenghtOfRight = new JLabel("Lenght of right flang:");
					lblLenghtOfRight.setBounds(5, 25, 328, 14);
					panelExtraction.add(lblLenghtOfRight);
				}
				{
					JLabel lblMaxDistanceBetween = new JLabel("Max distance between cistrons:");
					lblMaxDistanceBetween.setBounds(5, 45, 328, 14);
					panelExtraction.add(lblMaxDistanceBetween);
				}
				{
					JLabel lblMinLenghtOf = new JLabel("Min lenght of gene:");
					lblMinLenghtOf.setBounds(5, 65, 328, 14);
					panelExtraction.add(lblMinLenghtOf);
				}
				{
					JLabel lblStartCodons = new JLabel("Start codons:");
					lblStartCodons.setBounds(5, 105, 328, 14);
					panelExtraction.add(lblStartCodons);
				}
				{
					JLabel lblBadCodons = new JLabel("Bad codons:");
					lblBadCodons.setBounds(5, 145, 328, 14);
					panelExtraction.add(lblBadCodons);
				}
				{
					JLabel lblStopCodons = new JLabel("Stop codons:");
					lblStopCodons.setBounds(5, 185, 328, 14);
					panelExtraction.add(lblStopCodons);
				}
				
				lFlang = new JTextField();
				lFlang.setText(psNew.lFlang+"");
				lFlang.setBounds(343, 2, 86, 20);
				panelExtraction.add(lFlang);
				lFlang.setColumns(10);
				
				rFlang = new JTextField();
				rFlang.setText(psNew.rFlang+"");
				rFlang.setBounds(343, 22, 86, 20);
				panelExtraction.add(rFlang);
				rFlang.setColumns(10);
				
				distBOper = new JTextField();
				distBOper.setText(psNew.distBOper+"");
				distBOper.setBounds(343, 42, 86, 20);
				panelExtraction.add(distBOper);
				distBOper.setColumns(10);
				
				minLenGen = new JTextField();
				minLenGen.setText(psNew.minLenGen+"");
				minLenGen.setBounds(343, 62, 86, 20);
				panelExtraction.add(minLenGen);
				minLenGen.setColumns(10);

				checkStart = new JCheckBox("Check the presence of start codons");
				checkStart.addActionListener(new ActionListener() {
					public void actionPerformed(ActionEvent e) {
						startCodons.setEditable(checkStart.isSelected());
					}
				});
				checkStart.setBounds(0, 80, 333, 23);
				checkStart.setSelected(psNew.checkStart);
				panelExtraction.add(checkStart);
				
				deleteGwBC = new JCheckBox("Delete genes with bad codons");
				deleteGwBC.addActionListener(new ActionListener() {
					public void actionPerformed(ActionEvent e) {
						midlCodons.setEditable(deleteGwBC.isSelected());
					}
				});
				deleteGwBC.setBounds(0, 120, 333, 23);
				deleteGwBC.setSelected(psNew.deleteGwBC);
				panelExtraction.add(deleteGwBC);
				
				checkStop = new JCheckBox("Check the presence of stop codons");
				checkStop.addActionListener(new ActionListener() {
					public void actionPerformed(ActionEvent e) {
						stopCodons.setEditable(checkStop.isSelected());
					}
				});
				checkStop.setBounds(0, 160, 333, 23);
				checkStop.setSelected(psNew.checkStop);
				panelExtraction.add(checkStop);
				
				startCodons = new JTextField();
				{
					String[] arr1 = psNew.startCodons;
					String str1 = "";
					for (String s : arr1) {
						str1 += s.trim()+"; ";
					}
					startCodons.setText(str1.trim());
				}
				startCodons.setBounds(343, 102, 86, 20);
				startCodons.setEditable(checkStart.isSelected());
				startCodons.setColumns(10);
				panelExtraction.add(startCodons);
				
				midlCodons = new JTextField();
				{
					String[] arr1 = psNew.midlCodons;
					String str1 = "";
					for (String s : arr1) {
						str1 += s.trim()+"; ";
					}
					midlCodons.setText(str1.trim());
				}
				midlCodons.setBounds(343, 142, 86, 20);
				midlCodons.setEditable(deleteGwBC.isSelected());
				midlCodons.setColumns(10);
				panelExtraction.add(midlCodons);
								
				stopCodons = new JTextField();
				{
					String[] arr1 = psNew.stopCodons;
					String str1 = "";
					for (String s : arr1) {
						str1 += s.trim()+"; ";
					}
					stopCodons.setText(str1.trim());
				}
				stopCodons.setBounds(343, 182, 86, 20);
				stopCodons.setEditable(checkStop.isSelected());
				stopCodons.setColumns(10);
				panelExtraction.add(stopCodons);
								
				usePseudoGenes = new JCheckBox("Use pseudo genes in calculation");
				usePseudoGenes.setBounds(0, 200, 333, 23);
				usePseudoGenes.setSelected(psNew.usePseudoGenes);
				panelExtraction.add(usePseudoGenes);
			}
			//The Calculation tab.
			{
				JPanel panelCalculation = new JPanel();
				tabbedPane.addTab("Calculation", null, panelCalculation, null);
				panelCalculation.setLayout(null);
				
				autocounting = new JCheckBox("Use auto counting for flangs");
				autocounting.addActionListener(new ActionListener() {
					public void actionPerformed(ActionEvent e) {
						maxNumbCod.setEnabled(autocounting.isSelected());
						unUseCodonsL1.setEnabled(!autocounting.isSelected());
						unUseCodonsR1.setEnabled(!autocounting.isSelected());
					}
				});
				autocounting.setBounds(0, 0, 328, 23);
				autocounting.setSelected(psNew.autocounting);
				panelCalculation.add(autocounting);
				
				{
					JLabel lblMaxNumberOf = new JLabel("Max number of unuse codons on flangs:");
					lblMaxNumberOf.setBounds(5, 25, 323, 14);
					panelCalculation.add(lblMaxNumberOf);
				}
				{
					JLabel lblNumberOfUnuse = new JLabel("Number of unuse codons on the left flang:");
					lblNumberOfUnuse.setBounds(5, 45, 323, 14);
					panelCalculation.add(lblNumberOfUnuse);
				}
				{
					JLabel lblNumberOfUnuse_1 = new JLabel("Number of unuse codons on the right flang:");
					lblNumberOfUnuse_1.setBounds(5, 65, 323, 14);
					panelCalculation.add(lblNumberOfUnuse_1);
				}
				{
					JLabel lblTrainingSample = new JLabel("Training sample of genes:");
					lblTrainingSample.setBounds(5, 85, 323, 14);
					panelCalculation.add(lblTrainingSample);
				}
				
				maxNumbCod = new JSpinner();
				maxNumbCod.setModel(new SpinnerNumberModel(psNew.maxNumbCod, 0, 20, 1));
				maxNumbCod.setBounds(343, 22, 86, 20);
				maxNumbCod.setEnabled(autocounting.isSelected());
				panelCalculation.add(maxNumbCod);
				
				unUseCodonsL1 = new JSpinner();
				unUseCodonsL1.setModel(new SpinnerNumberModel(psNew.unUseCodonsL1, 0, 10, 1));
				unUseCodonsL1.setBounds(343, 42, 86, 20);
				unUseCodonsL1.setEnabled(!autocounting.isSelected());
				panelCalculation.add(unUseCodonsL1);
				
				unUseCodonsR1 = new JSpinner();
				unUseCodonsR1.setModel(new SpinnerNumberModel(psNew.unUseCodonsR1, 0, 10, 1));
				unUseCodonsR1.setBounds(343, 62, 86, 20);
				unUseCodonsR1.setEnabled(!autocounting.isSelected());
				panelCalculation.add(unUseCodonsR1);
								
				rdbtnNumberOfGenes = new JRadioButton("Number of genes");
				rdbtnNumberOfGenes.addActionListener(new ActionListener() {
					public void actionPerformed(ActionEvent arg0) {
						rdbtnPercentageOfAll.setSelected(!rdbtnNumberOfGenes.isSelected());
						NumbOfGenes.setEditable(rdbtnNumberOfGenes.isSelected());
						PercentageOfGenes.setEnabled(rdbtnPercentageOfAll.isSelected());
					}
				});
				rdbtnNumberOfGenes.setBounds(0, 105, 323, 23);
				rdbtnNumberOfGenes.setSelected(psNew.number_percentage);
				panelCalculation.add(rdbtnNumberOfGenes);
				
				rdbtnPercentageOfAll = new JRadioButton("Percentage of all genes (%)");
				rdbtnPercentageOfAll.addActionListener(new ActionListener() {
					public void actionPerformed(ActionEvent e) {
						rdbtnNumberOfGenes.setSelected(!rdbtnPercentageOfAll.isSelected());
						NumbOfGenes.setEditable(rdbtnNumberOfGenes.isSelected());
						PercentageOfGenes.setEnabled(rdbtnPercentageOfAll.isSelected());
					}
				});
				rdbtnPercentageOfAll.setBounds(0, 130, 328, 23);
				rdbtnPercentageOfAll.setSelected(!psNew.number_percentage);
				panelCalculation.add(rdbtnPercentageOfAll);
				
				NumbOfGenes = new JTextField();
				NumbOfGenes.setText(psNew.teachingGenesNumber+"");
				NumbOfGenes.setBounds(343, 106, 86, 20);
				NumbOfGenes.setEditable(rdbtnNumberOfGenes.isSelected());
				NumbOfGenes.setColumns(10);
				panelCalculation.add(NumbOfGenes);
				
				PercentageOfGenes = new JSpinner();
				PercentageOfGenes.setModel(new SpinnerNumberModel(1, 1, 50, 1));
				PercentageOfGenes.setValue(psNew.teachingGenesPercentage);
				PercentageOfGenes.setBounds(343, 131, 86, 20);
				PercentageOfGenes.setEnabled(rdbtnPercentageOfAll.isSelected());
				panelCalculation.add(PercentageOfGenes);
				
			}
			//"LCI" tab
			{
				JPanel panelLCI = new JPanel();
				tabbedPane.addTab("LCI", null, panelLCI, null);
				panelLCI.setLayout(null);
				{
					JLabel lblForCountingOf = new JLabel("For counting LCI1");
					lblForCountingOf.setBounds(5, 5, 328, 14);
					panelLCI.add(lblForCountingOf);
				}
				{
					JLabel lblMinRepeatsLenght = new JLabel("Min repeat's lenght:");
					lblMinRepeatsLenght.setBounds(5, 25, 328, 14);
					panelLCI.add(lblMinRepeatsLenght);
				}
				{
					JLabel lblMaxRepeatsLenght = new JLabel("Max repeat's lenght:");
					lblMaxRepeatsLenght.setBounds(5, 45, 328, 14);
					panelLCI.add(lblMaxRepeatsLenght);
				}
				{
					JLabel lblMinDistanceBetween = new JLabel("Min distance between repeats:");
					lblMinDistanceBetween.setBounds(5, 65, 328, 14);
					panelLCI.add(lblMinDistanceBetween);
				}
				{
					JLabel lblMaxDistanceBetween_1 = new JLabel("Max distance between repeats:");
					lblMaxDistanceBetween_1.setBounds(5, 85, 328, 14);
					panelLCI.add(lblMaxDistanceBetween_1);
				}
				
				smin1 = new JSpinner();
				smin1.setModel(new SpinnerNumberModel(psNew.smin1, 2, 10, 1));
				smin1.setBounds(343, 22, 86, 20);
				panelLCI.add(smin1);
				
				smax1 = new JSpinner();
				smax1.setModel(new SpinnerNumberModel(psNew.smax1, 2, 20, 1));
				smax1.setBounds(343, 42, 86, 20);
				panelLCI.add(smax1);
				
				lmin1 = new JSpinner();
				lmin1.setModel(new SpinnerNumberModel(psNew.lmin1, 3, 50, 1));
				lmin1.setBounds(343, 62, 86, 20);
				panelLCI.add(lmin1);
				
				lmax1 = new JSpinner();
				lmax1.setModel(new SpinnerNumberModel(psNew.lmax1, 3, 50, 1));
				lmax1.setBounds(343, 82, 86, 20);
				panelLCI.add(lmax1);
				
				{
					JLabel lblForCountingOf_1 = new JLabel("For counting LCI2");
					lblForCountingOf_1.setBounds(5, 105, 328, 14);
					panelLCI.add(lblForCountingOf_1);
				}
				{
					JLabel lblMinRepeatsLenght_1 = new JLabel("Min repeat's lenght:");
					lblMinRepeatsLenght_1.setBounds(5, 125, 328, 14);
					panelLCI.add(lblMinRepeatsLenght_1);
				}
				{
					JLabel lblMaxRepeatsLenght_1 = new JLabel("Max repeat's lenght:");
					lblMaxRepeatsLenght_1.setBounds(5, 145, 328, 14);
					panelLCI.add(lblMaxRepeatsLenght_1);
				}
				{
					JLabel lblMinDistanceBetween_1 = new JLabel("Min distance between repeats:");
					lblMinDistanceBetween_1.setBounds(5, 165, 328, 14);
					panelLCI.add(lblMinDistanceBetween_1);
				}
				{
					JLabel lblMaxDistanceBetween_2 = new JLabel("Max distance between repeats:");
					lblMaxDistanceBetween_2.setBounds(5, 185, 328, 14);
					panelLCI.add(lblMaxDistanceBetween_2);
				}
				
				smin2 = new JSpinner();
				smin2.setModel(new SpinnerNumberModel(psNew.smin2, 2, 10, 1));
				smin2.setBounds(343, 122, 86, 20);
				panelLCI.add(smin2);
				
				smax2 = new JSpinner();
				smax2.setModel(new SpinnerNumberModel(psNew.smax2, 2, 20, 1));
				smax2.setBounds(343, 142, 86, 20);
				panelLCI.add(smax2);
				
				lmin2 = new JSpinner();
				lmin2.setModel(new SpinnerNumberModel(psNew.lmin2, 3, 50, 1));
				lmin2.setBounds(343, 162, 86, 20);
				panelLCI.add(lmin2);
				
				lmax2 = new JSpinner();
				lmax2.setModel(new SpinnerNumberModel(psNew.lmax2, 3, 50, 1));
				lmax2.setBounds(343, 182, 86, 20);
				panelLCI.add(lmax2);
				
				{
					JLabel lblMinEnergyFor = new JLabel("Min energy for loop exist (kcal/mol):");
					lblMinEnergyFor.setBounds(5, 205, 328, 14);
					panelLCI.add(lblMinEnergyFor);
				}
				
				minloopEnergy = new JTextField();
				minloopEnergy.setText(psNew.minloopEnergy+"");
				minloopEnergy.setBounds(343, 202, 86, 20);
				panelLCI.add(minloopEnergy);
				minloopEnergy.setColumns(10);
				
			}
			//"User" tab
			{
				JPanel panelUser = new JPanel();
				tabbedPane.addTab("User", null, panelUser, null);
				panelUser.setLayout(null);
				{
					JLabel lblHereWillBe = new JLabel("Here will be user's parameters");
					lblHereWillBe.setBounds(5, 5, 328, 14);
					panelUser.add(lblHereWillBe);
				}
			}
			//"Results" tab
			{
				JPanel panelResults = new JPanel();
				tabbedPane.addTab("Results", null, panelResults, null);
				panelResults.setLayout(null);

				makeFileFreq = new JCheckBox("Make file with codon frequencies");
				makeFileFreq.setBounds(0, 0, 328, 23);
				makeFileFreq.setSelected(psNew.makeFileFreq);
				panelResults.add(makeFileFreq);
				
				makeFileWholeGenome = new JCheckBox("Make file with whole genome");
				makeFileWholeGenome.setBounds(0, 20, 328, 23);
				makeFileWholeGenome.setSelected(psNew.makeFileWholeGenom);
				panelResults.add(makeFileWholeGenome);
				
				searchOrthologs = new JCheckBox("Search orthologs");
				searchOrthologs.addActionListener(new ActionListener() {
					public void actionPerformed(ActionEvent arg0) {
						number_orthologs.setEnabled(searchOrthologs.isSelected());
					}
				});
				searchOrthologs.setBounds(0, 40, 328, 23);
				searchOrthologs.setSelected(psNew.searchOrthologs);
				panelResults.add(searchOrthologs);
				
				{
					JLabel lblNumberOfGenes = new JLabel("Number of genes for orthologs search:");
					lblNumberOfGenes.setBounds(5, 65, 323, 14);
					panelResults.add(lblNumberOfGenes);
				}
				
				number_orthologs = new JSpinner();
				number_orthologs.setEnabled(searchOrthologs.isSelected());
				number_orthologs.setModel(new SpinnerNumberModel(psNew.number_orthologs, 1, null, 1));
				number_orthologs.setBounds(343, 62, 86, 20);
				panelResults.add(number_orthologs);
								
				countLCIij = new JCheckBox("Count LCI for individual nucleotides");
				countLCIij.setBounds(0, 80, 328, 23);
				countLCIij.setSelected(psNew.countLCIij);
				panelResults.add(countLCIij);
				
				countRepeatsInGenes = new JCheckBox("Count average number of repeats in all genes");
				countRepeatsInGenes.setBounds(0, 100, 328, 23);
				countRepeatsInGenes.setSelected(psNew.countRepeatsInGenes);
				panelResults.add(countRepeatsInGenes);
				
				showRibPos = new JCheckBox("Show ribosomal positions in sorted lists of genes");
				showRibPos.setBounds(0, 120, 328, 23);
				showRibPos.setSelected(psNew.showRibPos);
				panelResults.add(showRibPos);
				
			}
		}
		{
			JPanel buttonPane = new JPanel();
			buttonPane.setLayout(new FlowLayout(FlowLayout.RIGHT));
			getContentPane().add(buttonPane, BorderLayout.SOUTH);
			{
				final JButton okButton = new JButton("OK");
				okButton.addActionListener(new ActionListener() {
					public void actionPerformed(ActionEvent e) {
						try {
							//"Extraction"
							{
								psNew.lFlang = Integer.parseInt(lFlang.getText());
								psNew.lFlangTrue = Integer.parseInt(lFlang.getText());
								if (psNew.lFlang < 0) throw new NumberFormatException("lFlang < 0");
								psNew.rFlang = Integer.parseInt(rFlang.getText());
								psNew.rFlangTrue = Integer.parseInt(rFlang.getText());
								if (psNew.rFlang < 0) throw new NumberFormatException("rFlang < 0");
								psNew.distBOper = Integer.parseInt(distBOper.getText());
								if (psNew.distBOper < 0) throw new NumberFormatException("distBOper < 0");
								psNew.minLenGen = Integer.parseInt(minLenGen.getText());
								if (psNew.minLenGen < 0) throw new NumberFormatException("minLenGen < 0");
								psNew.checkStart = checkStart.isSelected();
								psNew.deleteGwBC = deleteGwBC.isSelected();
								psNew.checkStop = checkStop.isSelected();
								{
									String str1 = startCodons.getText();
									Matcher m = CODONPATTERN.matcher(str1);
							        if (!m.matches()) {
							        	if (checkStart.isSelected()) throw new NumberFormatException("startCodons are wrong");
							        	else {
							        		psNew.startCodons = psDefault.startCodons;
							        		StringBuilder str2 = new StringBuilder();
											for (String s : psDefault.startCodons) {
												str2.append(s+";");
											}
											startCodons.setText(str2.toString());
											startCodons.setEditable(checkStart.isSelected());
							        	}
							        }
							        else {
							        	String[] arr1 = str1.trim().split(";");
										psNew.startCodons = arr1;
							        }
								}
								{
									String str1 = midlCodons.getText();
									Matcher m = CODONPATTERN.matcher(str1);
							        if (!m.matches()) {
							        	if (deleteGwBC.isSelected()) throw new NumberFormatException("midlCodons are wrong");
							        	else {
							        		psNew.midlCodons = psDefault.midlCodons;
							        		StringBuilder str2 = new StringBuilder();
											for (String s : psDefault.midlCodons) {
												str2.append(s+";");
											}
											midlCodons.setText(str2.toString());
											midlCodons.setEditable(deleteGwBC.isSelected());
							        	}
							        }
							        else {
							        	String[] arr1 = str1.trim().split(";");
										psNew.midlCodons = arr1;
							        }
								}
								{
									String str1 = stopCodons.getText();
									Matcher m = CODONPATTERN.matcher(str1);
							        if (!m.matches()) {
							        	if (checkStop.isSelected()) throw new NumberFormatException("stopCodons are wrong");
							        	else {
							        		psNew.stopCodons = psDefault.stopCodons;
							        		StringBuilder str2 = new StringBuilder();
											for (String s : psDefault.stopCodons) {
												str2.append(s+";");
											}
											stopCodons.setText(str2.toString());
											stopCodons.setEditable(checkStop.isSelected());
							        	}
							        }
							        else {
							        	String[] arr1 = str1.trim().split(";");
										psNew.stopCodons = arr1;
							        }
								}
								psNew.usePseudoGenes = usePseudoGenes.isSelected();
								
							}
							//"Calculation"
							{
								psNew.autocounting = autocounting.isSelected();
								psNew.maxNumbCod = (Integer)maxNumbCod.getValue();
								psNew.unUseCodonsL1 = (Integer)unUseCodonsL1.getValue();
								psNew.unUseCodonsR1 = (Integer)unUseCodonsR1.getValue();
								psNew.number_percentage = rdbtnNumberOfGenes.isSelected();
								psNew.teachingGenesNumber = Integer.parseInt(NumbOfGenes.getText());
								if (psNew.teachingGenesNumber <= 0) throw new NumberFormatException("teachingGenesNumber <= 0");
								psNew.teachingGenesPercentage = (Integer)PercentageOfGenes.getValue();
							}
							//"LCI"
							{
								psNew.smin1 = (Integer)smin1.getValue();
								psNew.smax1 = (Integer)smax1.getValue();
								if (psNew.smin1 > psNew.smax1) throw new NumberFormatException("smin1 > smax1");
								psNew.lmin1 = (Integer)lmin1.getValue();
								psNew.lmax1 = (Integer)lmax1.getValue();
								if (psNew.lmin1 > psNew.lmax1) throw new NumberFormatException("lmin1 > lmax1");
								psNew.smin2 = (Integer)smin2.getValue();
								psNew.smax2 = (Integer)smax2.getValue();
								if (psNew.smin2 > psNew.smax2) throw new NumberFormatException("smin2 > smax2");
								psNew.lmin2 = (Integer)lmin2.getValue();
								psNew.lmax2 = (Integer)lmax2.getValue();
								if (psNew.lmin2 > psNew.lmax2) throw new NumberFormatException("lmin2 > lmax2");
								psNew.minloopEnergy = Double.parseDouble(minloopEnergy.getText());
							}
							//"User"
							{
								//Currently empty.
							}
							//"Results"
							{
								psNew.makeFileFreq = makeFileFreq.isSelected();
								psNew.makeFileWholeGenom = makeFileWholeGenome.isSelected();
								psNew.searchOrthologs = searchOrthologs.isSelected();
								psNew.number_orthologs = (Integer)number_orthologs.getValue();
								psNew.countLCIij = countLCIij.isSelected();
								psNew.countRepeatsInGenes = countRepeatsInGenes.isSelected();
								psNew.showRibPos = showRibPos.isSelected();
							}
							ParametersWindow.this.dispose();
						} catch (NumberFormatException nfe) {
							nfe.printStackTrace();
							JOptionPane.showMessageDialog(okButton, "Warning! Some parameters aren't correct.\n" +
									nfe.getMessage()+".", "Warning", JOptionPane.WARNING_MESSAGE);
						}						
					}
				});
				
				JButton btnDefault = new JButton("Default");
				btnDefault.addActionListener(new ActionListener() {
					public void actionPerformed(ActionEvent e) {
						//"Extraction"
						{
							lFlang.setText(psDefault.lFlang+"");
							rFlang.setText(psDefault.rFlang+"");
							distBOper.setText(psDefault.distBOper+"");
							minLenGen.setText(psDefault.minLenGen+"");
							checkStart.setSelected(psDefault.checkStart);
							deleteGwBC.setSelected(psDefault.deleteGwBC);
							checkStop.setSelected(psDefault.checkStop);
							{
								StringBuilder str1 = new StringBuilder();
								for (String s : psDefault.startCodons) {
									str1.append(s+";");
								}
								startCodons.setText(str1.toString());
								startCodons.setEditable(checkStart.isSelected());
							}
							{
								StringBuilder str1 = new StringBuilder();
								for (String s : psDefault.midlCodons) {
									str1.append(s+";");
								}
								midlCodons.setText(str1.toString());
								midlCodons.setEditable(deleteGwBC.isSelected());
							}
							{
								StringBuilder str1 = new StringBuilder();
								for (String s : psDefault.stopCodons) {
									str1.append(s+";");
								}
								stopCodons.setText(str1.toString());
								stopCodons.setEditable(checkStop.isSelected());
							}
							usePseudoGenes.setSelected(psDefault.usePseudoGenes);
						}
						//"Calculation"
						{
							autocounting.setSelected(psDefault.autocounting);
							maxNumbCod.setValue(psDefault.maxNumbCod);
							maxNumbCod.setEnabled(autocounting.isSelected());
							unUseCodonsL1.setValue(psDefault.unUseCodonsL1);
							unUseCodonsL1.setEnabled(!autocounting.isSelected());
							unUseCodonsR1.setValue(psDefault.unUseCodonsR1);
							unUseCodonsR1.setEnabled(!autocounting.isSelected());
							rdbtnNumberOfGenes.setSelected(psDefault.number_percentage);
							rdbtnPercentageOfAll.setSelected(!psDefault.number_percentage);
							NumbOfGenes.setText(psDefault.teachingGenesNumber+"");
							NumbOfGenes.setEditable(rdbtnNumberOfGenes.isSelected());
							PercentageOfGenes.setValue(psDefault.teachingGenesPercentage);
							PercentageOfGenes.setEnabled(rdbtnPercentageOfAll.isSelected());
						}
						//"LCI"
						{
							smin1.setValue(psDefault.smin1);
							smax1.setValue(psDefault.smax1);
							lmin1.setValue(psDefault.lmin1);
							lmax1.setValue(psDefault.lmax1);
							smin2.setValue(psDefault.smin2);
							smax2.setValue(psDefault.smax2);
							lmin2.setValue(psDefault.lmin2);
							lmax2.setValue(psDefault.lmax2);
							minloopEnergy.setText(psDefault.minloopEnergy+"");
						}
						//"User"
						{
							//Currently empty.
						}
						//"Results"
						{
							makeFileFreq.setSelected(psDefault.makeFileFreq);
							makeFileWholeGenome.setSelected(psDefault.makeFileWholeGenom);
							searchOrthologs.setSelected(psDefault.searchOrthologs);
							number_orthologs.setValue(psDefault.number_orthologs);
							number_orthologs.setEnabled(searchOrthologs.isSelected());
							countLCIij.setSelected(psDefault.countLCIij);
							countRepeatsInGenes.setSelected(psDefault.countRepeatsInGenes);
							showRibPos.setSelected(psDefault.showRibPos);
						}
					}
				});
				buttonPane.add(btnDefault);
				okButton.setActionCommand("OK");
				buttonPane.add(okButton);
				getRootPane().setDefaultButton(okButton);
			}
			{
				JButton cancelButton = new JButton("Cancel");
				cancelButton.addActionListener(new ActionListener() {
					public void actionPerformed(ActionEvent e) {
						ParametersWindow.this.dispose();
					}
				});
				cancelButton.setActionCommand("Cancel");
				buttonPane.add(cancelButton);
			}
		}
	}
}
