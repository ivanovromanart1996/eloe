package org.eloe.parameters;

public interface PAutoCounting {

	boolean isAutocounting();

	int getMaxNumbCod();

	void setUnUseCodonsL1(int unUseCodonsL1);

	void setUnUseCodonsR1(int unUseCodonsR1);

	double getMinKc();

	double getVe();

}