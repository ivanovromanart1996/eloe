package org.eloe.parameters;

import org.apache.log4j.Logger;

public class Parameters implements PExtraction, PProcessing, PLCI1, PLCI2, PTa, PGeneToString, PAutoCounting, PEEI, POrthologs, PFreq {
	
	private static final Logger LOG = Logger.getLogger(Parameters.class);
	
	private String context = "CDS ";
	private int lFlang = 600;
	private int lFlangTrue = 600;
	private int rFlang = 600;
	private int rFlangTrue = 600;
	private int distBOper = 40;
	private int minLenGen = 90;
	private boolean checkStart = false;
	private String[] startCodons = {"atg", "gtg", "ttg"};
	private boolean deleteGwBC = false;
	private String[] midlCodons = {"tag", "taa"};
	private boolean checkStop = false;
	private String[] stopCodons = {"tag", "taa"};
	private boolean usePseudoGenes = false;
	private boolean makeFileFreq = false;
	private boolean makeFileWholeGenom = false;
	private boolean searchOrthologs = false;
	private int number_orthologs = 25;
	private boolean autocounting = false;
	private int maxNumbCod = 20;
	private int unUseCodonsL1 = 1;
	private int unUseCodonsR1 = 1;
	private int smin1 = 3;
	private int smax1 = 6;
	private int lmin1 = 3;
	private int lmax1 = 50;
	private int smin2 = 3;
	private int smax2 = 6;
	private int lmin2 = 3;
	private int lmax2 = 50;
	private double minloopEnergy = 0.0;
	private boolean countLCIij = false;
	private boolean countRepeatsInGenes = false;
	private boolean number_percentage = true;
	private int teachingGenesNumber = 150;
	private int teachingGenesPercentage = 10;
	private double minKc = 0.1;
	private double ve = 1.;
	private int forSigma1 = 2;
	private int forSigma2 = 2;
	private boolean useUserParam1 = false;
	private double m11 = 1.;
	private double sigma11 = 0.1;
	private double m22 = 1.;
	private double sigma22 = 0.1;
	private boolean useUserParam2 = false;
	private double tmin1 = 0.1;
	private double tmax1 = 1.;
	private boolean useRibGenesAsHighExpressed1 = true;
	private boolean useOtherGenesAsHighExpressed1 = false;
	
	public Parameters(final ParametersStructure ps) {
		setContext(ps.context);
		setlFlang(ps.lFlang);
		setlFlangTrue(ps.lFlangTrue);
		setrFlang(ps.rFlang);
		setrFlangTrue(ps.rFlangTrue);
		setDistBOper(ps.distBOper);
		setMinLenGen(ps.minLenGen);
		setCheckStart(ps.checkStart);
		setStartCodons(ps.startCodons);
		setDeleteGwBC(ps.deleteGwBC);
		setMidlCodons(ps.midlCodons);
		setCheckStop(ps.checkStop);
		setStopCodons(ps.stopCodons);
		setUsePseudoGenes(ps.usePseudoGenes);
		setMakeFileFreq(ps.makeFileFreq);
		setMakeFileWholeGenom(ps.makeFileWholeGenom);
		setSearchOrthologs(ps.searchOrthologs);
		setNumber_orthologs(ps.number_orthologs);
		setAutocounting(ps.autocounting);
		setMaxNumbCod(ps.maxNumbCod);
		setUnUseCodonsL1(ps.unUseCodonsL1);
		setUnUseCodonsR1(ps.unUseCodonsR1);
		setSmin1(ps.smin1);
		setSmax1(ps.smax1);
		setLmin1(ps.lmin1);
		setLmax1(ps.lmax1);
		setSmin2(ps.smin2);
		setSmax2(ps.smax2);
		setLmin2(ps.lmin2);
		setLmax2(ps.lmax2);
		setMinloopEnergy(ps.minloopEnergy);
		setCountLCIij(ps.countLCIij);
		setCountRepeatsInNoRib(ps.countRepeatsInGenes);
		setNumberPercentage(ps.number_percentage);
		setTeachingGenesNumber(ps.teachingGenesNumber);
		setTeachingGenesPercentage(ps.teachingGenesPercentage);
		setMinKc(ps.minKc);
		setVe(ps.ve);
		setForSigma1(ps.forSigma1);
		setForSigma2(ps.forSigma2);
		setUseUserParam1(ps.useUserParam1);
		setM11(ps.m11);
		setSigma11(ps.sigma11);
		setM22(ps.m22);
		setSigma22(ps.sigma22);
		setUseUserParam2(ps.useUserParam2);
		setTmin1(ps.tmin1);
		setTmax1(ps.tmin1);
		setUseRibGenesAsHighExpressed1(ps.useRibGenesAsHighExpressed1);
		setUseOtherGenesAsHighExpressed1(ps.useOtherGenesAsHighExpressed1);
	}

	/* (non-Javadoc)
	 * @see org.eloe.parameters.PExtraction#getContext()
	 */
	@Override
	public String getContext() {
		return context;
	}

	public void setContext(String context) {
		this.context = context;
	}

	/* (non-Javadoc)
	 * @see org.eloe.parameters.PExtraction#getlFlang()
	 */
	@Override
	public int getlFlang() {
		return lFlang;
	}

	public void setlFlang(int lFlang) {
		this.lFlang = lFlang;
	}

	/* (non-Javadoc)
	 * @see org.eloe.parameters.PGeneToString#getlFlangTrue()
	 */
	@Override
	public int getlFlangTrue() {
		return lFlangTrue;
	}

	public void setlFlangTrue(int lFlangTrue) {
		this.lFlangTrue = lFlangTrue;
	}

	/* (non-Javadoc)
	 * @see org.eloe.parameters.PExtraction#getrFlang()
	 */
	@Override
	public int getrFlang() {
		return rFlang;
	}

	public void setrFlang(int rFlang) {
		this.rFlang = rFlang;
	}

	/* (non-Javadoc)
	 * @see org.eloe.parameters.PGeneToString#getrFlangTrue()
	 */
	@Override
	public int getrFlangTrue() {
		return rFlangTrue;
	}

	public void setrFlangTrue(int rFlangTrue) {
		this.rFlangTrue = rFlangTrue;
	}

	/* (non-Javadoc)
	 * @see org.eloe.parameters.PExtraction#getDistBOper()
	 */
	@Override
	public int getDistBOper() {
		return distBOper;
	}

	public void setDistBOper(int distBOper) {
		this.distBOper = distBOper;
	}

	/* (non-Javadoc)
	 * @see org.eloe.parameters.PExtraction#getMinLenGen()
	 */
	@Override
	public int getMinLenGen() {
		return minLenGen;
	}

	public void setMinLenGen(int minLenGen) {
		this.minLenGen = minLenGen;
	}

	/* (non-Javadoc)
	 * @see org.eloe.parameters.PProcessing#isCheckStart()
	 */
	@Override
	public boolean isCheckStart() {
		return checkStart;
	}

	public void setCheckStart(boolean checkStart) {
		this.checkStart = checkStart;
	}

	/* (non-Javadoc)
	 * @see org.eloe.parameters.PProcessing#getStartCodons()
	 */
	@Override
	public String[] getStartCodons() {
		return startCodons;
	}

	public void setStartCodons(String[] startCodons) {
		this.startCodons = startCodons;
	}

	/* (non-Javadoc)
	 * @see org.eloe.parameters.PProcessing#isDeleteGwBC()
	 */
	@Override
	public boolean isDeleteGwBC() {
		return deleteGwBC;
	}

	public void setDeleteGwBC(boolean deleteGwBC) {
		this.deleteGwBC = deleteGwBC;
	}

	/* (non-Javadoc)
	 * @see org.eloe.parameters.PProcessing#getMidlCodons()
	 */
	@Override
	public String[] getMidlCodons() {
		return midlCodons;
	}

	public void setMidlCodons(String[] midlCodons) {
		this.midlCodons = midlCodons;
	}

	/* (non-Javadoc)
	 * @see org.eloe.parameters.PProcessing#isCheckStop()
	 */
	@Override
	public boolean isCheckStop() {
		return checkStop;
	}

	public void setCheckStop(boolean checkStop) {
		this.checkStop = checkStop;
	}

	/* (non-Javadoc)
	 * @see org.eloe.parameters.PProcessing#getStopCodons()
	 */
	@Override
	public String[] getStopCodons() {
		return stopCodons;
	}

	public void setStopCodons(String[] stopCodons) {
		this.stopCodons = stopCodons;
	}

	/* (non-Javadoc)
	 * @see org.eloe.parameters.PExtraction#isUsePseudoGenes()
	 * @see org.eloe.parameters.PProcessing#isUsePseudoGenes()
	 */
	@Override
	public boolean isUsePseudoGenes() {
		return usePseudoGenes;
	}

	public void setUsePseudoGenes(boolean usePseudoGenes) {
		this.usePseudoGenes = usePseudoGenes;
	}

	/* (non-Javadoc)
	 * @see org.eloe.parameters.PFreq#isMakeFileFreq()
	 */
	@Override
	public boolean isMakeFileFreq() {
		return makeFileFreq;
	}

	public void setMakeFileFreq(boolean makeFileFreq) {
		this.makeFileFreq = makeFileFreq;
	}

	/* (non-Javadoc)
	 * @see org.eloe.parameters.PExtraction#isMakeFileWholeGenom()
	 */
	@Override
	public boolean isMakeFileWholeGenome() {
		return makeFileWholeGenom;
	}

	public void setMakeFileWholeGenom(boolean makeFileWholeGenom) {
		this.makeFileWholeGenom = makeFileWholeGenom;
	}

	/* (non-Javadoc)
	 * @see org.eloe.parameters.POrthologs#isSearchOrthologs()
	 */
	@Override
	public boolean isSearchOrthologs() {
		return searchOrthologs;
	}

	public void setSearchOrthologs(boolean searchOrthologs) {
		this.searchOrthologs = searchOrthologs;
	}

	/* (non-Javadoc)
	 * @see org.eloe.parameters.POrthologs#getNumber_orthologs()
	 */
	@Override
	public int getNumber_orthologs() {
		return number_orthologs;
	}

	public void setNumber_orthologs(int number_orthologs) {
		this.number_orthologs = number_orthologs;
	}

	/* (non-Javadoc)
	 * @see org.eloe.parameters.PAutoCounting#isAutocounting()
	 */
	@Override
	public boolean isAutocounting() {
		return autocounting;
	}

	public void setAutocounting(boolean autocounting) {
		this.autocounting = autocounting;
	}

	/* (non-Javadoc)
	 * @see org.eloe.parameters.PExtraction#getMaxNumbCod()
	 * @see org.eloe.parameters.PAutoCounting#getMaxNumbCod()
	 */
	@Override
	public int getMaxNumbCod() {
		return maxNumbCod;
	}

	public void setMaxNumbCod(int maxNumbCod) {
		this.maxNumbCod = maxNumbCod;
	}

	/* (non-Javadoc)
	 * @see org.eloe.parameters.PEEI#getUnUseCodonsL1()
	 * @see org.eloe.parameters.PFreq#getUnUseCodonsL1()
	 */
	/* (non-Javadoc)
	 * 
	 */
	@Override
	public int getUnUseCodonsL1() {
		return unUseCodonsL1;
	}

	/* (non-Javadoc)
	 * @see org.eloe.parameters.PAutoCounting#getUnUseCodonsL1()
	 */
	@Override
	public void setUnUseCodonsL1(int unUseCodonsL1) {
		this.unUseCodonsL1 = unUseCodonsL1;
	}

	/* (non-Javadoc)
	 * @see org.eloe.parameters.PEEI#getUnUseCodonsR1()
	 * @see org.eloe.parameters.PFreq#getUnUseCodonsR1()
	 */
	@Override
	public int getUnUseCodonsR1() {
		return unUseCodonsR1;
	}

	/* (non-Javadoc)
	 * @see org.eloe.parameters.PAutoCounting#getUnUseCodonsR1()
	 */
	@Override
	public void setUnUseCodonsR1(int unUseCodonsR1) {
		this.unUseCodonsR1 = unUseCodonsR1;
	}

	/* (non-Javadoc)
	 * @see org.eloe.parameters.PLCI1#getSmin1()
	 */
	@Override
	public int getSmin1() {
		return smin1;
	}

	public void setSmin1(int smin1) {
		this.smin1 = smin1;
	}

	/* (non-Javadoc)
	 * @see org.eloe.parameters.PLCI1#getSmax1()
	 */
	@Override
	public int getSmax1() {
		return smax1;
	}

	public void setSmax1(int smax1) {
		this.smax1 = smax1;
	}

	/* (non-Javadoc)
	 * @see org.eloe.parameters.PLCI1#getLmin1()
	 */
	@Override
	public int getLmin1() {
		return lmin1;
	}

	public void setLmin1(int lmin1) {
		this.lmin1 = lmin1;
	}

	/* (non-Javadoc)
	 * @see org.eloe.parameters.PLCI1#getLmax1()
	 */
	@Override
	public int getLmax1() {
		return lmax1;
	}

	public void setLmax1(int lmax1) {
		this.lmax1 = lmax1;
	}

	/* (non-Javadoc)
	 * @see org.eloe.parameters.PLCI2#getSmin2()
	 */
	@Override
	public int getSmin2() {
		return smin2;
	}

	public void setSmin2(int smin2) {
		this.smin2 = smin2;
	}

	/* (non-Javadoc)
	 * @see org.eloe.parameters.PLCI2#getSmax2()
	 */
	@Override
	public int getSmax2() {
		return smax2;
	}

	public void setSmax2(int smax2) {
		this.smax2 = smax2;
	}

	/* (non-Javadoc)
	 * @see org.eloe.parameters.PLCI2#getLmin2()
	 */
	@Override
	public int getLmin2() {
		return lmin2;
	}

	public void setLmin2(int lmin2) {
		this.lmin2 = lmin2;
	}

	/* (non-Javadoc)
	 * @see org.eloe.parameters.PLCI2#getLmax2()
	 */
	@Override
	public int getLmax2() {
		return lmax2;
	}

	public void setLmax2(int lmax2) {
		this.lmax2 = lmax2;
	}

	/* (non-Javadoc)
	 * @see org.eloe.parameters.PLCI2#getMinloopEnergy()
	 */
	@Override
	public double getMinloopEnergy() {
		return minloopEnergy;
	}

	public void setMinloopEnergy(double minloopEnergy) {
		this.minloopEnergy = minloopEnergy;
	}

	public boolean isCountLCIij() {
		return countLCIij;
	}

	public void setCountLCIij(boolean countLCIij) {
		this.countLCIij = countLCIij;
	}
	
	public boolean isCountRepeatsInNoRib() {
		return countRepeatsInGenes;
	}

	public void setCountRepeatsInNoRib(boolean countRepeatsInNoRib) {
		this.countRepeatsInGenes = countRepeatsInNoRib;
	}

	/* (non-Javadoc)
	 * @see org.eloe.parameters.PTa#getNumberPercentage()
	 */
	@Override
	public boolean getNumberPercentage() {
		return number_percentage;
	}

	public void setNumberPercentage(boolean number_percentage) {
		this.number_percentage = number_percentage;
	}

	@Override
	public int getTeachingGenesNumber() {
		return teachingGenesNumber;
	}

	public void setTeachingGenesNumber(int teachingGenes) {
		this.teachingGenesNumber = teachingGenes;
	}
	
	@Override
	public int getTeachingGenesPercentage() {
		return teachingGenesPercentage;
	}
	
	public void setTeachingGenesPercentage(int teachingGenesPercentage) {
		this.teachingGenesPercentage = teachingGenesPercentage;
	}

	/* (non-Javadoc)
	 * @see org.eloe.parameters.PTa#getMinKc()
	 * @see org.eloe.parameters.PAutoCounting#getMinKc()
	 * @see org.eloe.parameters.PEEI#getMinKc()
	 */
	@Override
	public double getMinKc() {
		return minKc;
	}

	public void setMinKc(double minKc) {
		this.minKc = minKc;
	}

	/* (non-Javadoc)
	 * @see org.eloe.parameters.PTa#getVe()
	 * @see org.eloe.parameters.PAutoCounting#getVe()
	 * @see org.eloe.parameters.PEEI#getVe()
	 */
	@Override
	public double getVe() {
		return ve;
	}

	public void setVe(double ve) {
		this.ve = ve;
	}

	/* (non-Javadoc)
	 * @see org.eloe.parameters.PTa#getForSigma1()
	 */
	@Override
	public int getForSigma1() {
		return forSigma1;
	}

	public void setForSigma1(int forSigma1) {
		this.forSigma1 = forSigma1;
	}

	/* (non-Javadoc)
	 * @see org.eloe.parameters.PTa#getForSigma2()
	 */
	@Override
	public int getForSigma2() {
		return forSigma2;
	}

	public void setForSigma2(int forSigma2) {
		this.forSigma2 = forSigma2;
	}

	public boolean isUseUserParam1() {
		return useUserParam1;
	}

	public void setUseUserParam1(boolean useUserParam1) {
		this.useUserParam1 = useUserParam1;
	}

	public double getM11() {
		return m11;
	}

	public void setM11(double m11) {
		this.m11 = m11;
	}

	public double getSigma11() {
		return sigma11;
	}

	public void setSigma11(double sigma11) {
		this.sigma11 = sigma11;
	}

	public double getM22() {
		return m22;
	}

	public void setM22(double m22) {
		this.m22 = m22;
	}

	public double getSigma22() {
		return sigma22;
	}

	public void setSigma22(double sigma22) {
		this.sigma22 = sigma22;
	}

	/* (non-Javadoc)
	 * @see org.eloe.parameters.PEEI#isUseUserParam2()
	 */
	@Override
	public boolean isUseUserParam2() {
		return useUserParam2;
	}

	public void setUseUserParam2(boolean useUserParam2) {
		this.useUserParam2 = useUserParam2;
	}

	/* (non-Javadoc)
	 * @see org.eloe.parameters.PEEI#getTmin1()
	 */
	@Override
	public double getTmin1() {
		return tmin1;
	}

	public void setTmin1(double tmin1) {
		this.tmin1 = tmin1;
	}

	/* (non-Javadoc)
	 * @see org.eloe.parameters.PEEI#getTmax1()
	 */
	@Override
	public double getTmax1() {
		return tmax1;
	}

	public void setTmax1(double tmax1) {
		this.tmax1 = tmax1;
	}

	/* (non-Javadoc)
	 * @see org.eloe.parameters.PExtraction#isUseRibGenesAsHighExpressed1()
	 */
	@Override
	public boolean isUseRibGenesAsHighExpressed1() {
		return useRibGenesAsHighExpressed1;
	}

	public void setUseRibGenesAsHighExpressed1(boolean useRibGenesAsHighExpressed1) {
		this.useRibGenesAsHighExpressed1 = useRibGenesAsHighExpressed1;
	}

	/* (non-Javadoc)
	 * @see org.eloe.parameters.PExtraction#isUseOtherGenesAsHighExpressed1()
	 */
	@Override
	public boolean isUseOtherGenesAsHighExpressed1() {
		return useOtherGenesAsHighExpressed1;
	}

	public void setUseOtherGenesAsHighExpressed1(
			boolean useOtherGenesAsHighExpressed1) {
		this.useOtherGenesAsHighExpressed1 = useOtherGenesAsHighExpressed1;
	}	
	
	//TODO
	public final int teachingGenes(final int number) {
		//Подсчет количества генов в выборке, если указано точное количество
		int number2;
		if (this.number_percentage) {
			number2 = this.teachingGenesNumber;
			//Подсчет количества генов в выборке, если указано точное количество, 
			//но удвоенное оно больше общего числа генов в геноме
			if (2 * number2 > number) {
				//Берем 10%
				number2 = (int) Math.round(number / 10);
				//Вывод в консоль сообщения об ошибке
				LOG.warn("WARNING! Number of genes is bigger than half of number of genes in genome! (For counting 10% of genome were taken)");
			}
		} else { //Подсчет количества генов в выборке, если указано количество в процентах от полного числа генов в геноме
			number2 = (int) Math.round(number * this.teachingGenesPercentage / 100.0);
		}
		return number2;
	}
}
