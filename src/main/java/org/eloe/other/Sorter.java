package org.eloe.other;

import java.util.Collections;
import java.util.Comparator;

import org.eloe.organism.Gene;
import org.eloe.organism.Organism;

public class Sorter {
	//Comparator for eei
	private static final Comparator<Gene> eeiUp = new Comparator<Gene>() {
		@Override
		public int compare(Gene arg0, Gene arg1) {
			int a = 0;
			double delta = arg0.getEei()-arg1.getEei();
			if (delta < 0) {
				a = -1;
			} else {
				if (delta > 0) {
					a = 1;
				}
			}
			return a;
		}
	};
	
	//Comparator for eei1
	private static final Comparator<Gene> eei1Up = new Comparator<Gene>() {
		@Override
		public int compare(Gene arg0, Gene arg1) {
			int a = 0;
			double delta = arg0.getEeiMass()[0]-arg1.getEeiMass()[0];
			if (delta < 0) {
				a = -1;
			} else {
				if (delta > 0) {
					a = 1;
				}
			}
			return a;
		}
	};
	
	//Comparator for eei2
	private static final Comparator<Gene> eei2Up = new Comparator<Gene>() {
		@Override
		public int compare(Gene arg0, Gene arg1) {
			int a = 0;
			double delta = arg0.getEeiMass()[1]-arg1.getEeiMass()[1];
			if (delta < 0) {
				a = -1;
			} else {
				if (delta > 0) {
					a = 1;
				}
			}
			return a;
		}
	};
	
	//Comparator for eei3
	private static final Comparator<Gene> eei3Up = new Comparator<Gene>() {
		@Override
		public int compare(Gene arg0, Gene arg1) {
			int a = 0;
			double delta = arg0.getEeiMass()[2]-arg1.getEeiMass()[2];
			if (delta < 0) {
				a = -1;
			} else {
				if (delta > 0) {
					a = 1;
				}
			}
			return a;
		}
	};
	
	//Comparator for eei4
	private static final Comparator<Gene> eei4Up = new Comparator<Gene>() {
		@Override
		public int compare(Gene arg0, Gene arg1) {
			int a = 0;
			double delta = arg0.getEeiMass()[3]-arg1.getEeiMass()[3];
			if (delta < 0) {
				a = -1;
			} else {
				if (delta > 0) {
					a = 1;
				}
			}
			return a;
		}
	};
	
	//Comparator for eei5
	private static final Comparator<Gene> eei5Up = new Comparator<Gene>() {
		@Override
		public int compare(Gene arg0, Gene arg1) {
			int a = 0;
			double delta = arg0.getEeiMass()[4]-arg1.getEeiMass()[4];
			if (delta < 0) {
				a = -1;
			} else {
				if (delta > 0) {
					a = 1;
				}
			}
			return a;
		}
	};
	
	//Comparator for gene number
	private static final Comparator<Gene> geneNumberUp = new Comparator<Gene>() {
		@Override
		public int compare(Gene arg0, Gene arg1) {
			int a = 0;
			double delta = arg0.getNumber()-arg1.getNumber();
			if (delta < 0) {
				a = -1;
			} else {
				if (delta > 0) {
					a = 1;
				}
			}
			return a;
		}
	};
	
	//Comparator for lci1
	private static final Comparator<Gene> lci1Up = new Comparator<Gene>() {
		@Override
		public int compare(Gene arg0, Gene arg1) {
			int a = 0;
			double delta = arg0.getLci1()-arg1.getLci1();
			if (delta < 0) {
				a = -1;
			} else {
				if (delta > 0) {
					a = 1;
				}
			}
			return a;
		}
	};
	
	//Comparator for lci2
	private static final Comparator<Gene> lci2Up = new Comparator<Gene>() {
		@Override
		public int compare(Gene arg0, Gene arg1) {
			int a = 0;
			double delta = arg0.getLci2()-arg1.getLci2();
			if (delta < 0) {
				a = -1;
			} else {
				if (delta > 0) {
					a = 1;
				}
			}
			return a;
		}
	};
	
	//Comparator for tai
	private static final Comparator<Gene> taiUp = new Comparator<Gene>() {
		@Override
		public int compare(Gene arg0, Gene arg1) {
			int a = 0;
			double delta = arg0.getTai()-arg1.getTai();
			if (delta < 0) {
				a = -1;
			} else {
				if (delta > 0) {
					a = 1;
				}
			}
			return a;
		}
	};
	
	//Comparator for tei
	private static final Comparator<Gene> teiUp = new Comparator<Gene>() {
		@Override
		public int compare(Gene arg0, Gene arg1) {
			int a = 0;
			double delta = arg0.getTei()-arg1.getTei();
			if (delta < 0) {
				a = -1;
			} else {
				if (delta > 0) {
					a = 1;
				}
			}
			return a;
		}
	};
	
	//Comparator for p_i1
	private static final Comparator<Gene> pi1Up = new Comparator<Gene>() {
		@Override
		public int compare(Gene arg0, Gene arg1) {
			int a = 0;
			double delta = arg0.getPi1()-arg1.getPi1();
			if (delta < 0) {
				a = -1;
			} else {
				if (delta > 0) {
					a = 1;
				}
			}
			return a;
		}
	};
	
	//Comparator for p_i2
	private static final Comparator<Gene> pi2Up = new Comparator<Gene>() {
		@Override
		public int compare(Gene arg0, Gene arg1) {
			int a = 0;
			double delta = arg0.getPi2()-arg1.getPi2();
			if (delta < 0) {
				a = -1;
			} else {
				if (delta > 0) {
					a = 1;
				}
			}
			return a;
		}
	};
	
	//Comparator for number of repetitions
	private static final Comparator<Gene> repeatsUp = new Comparator<Gene>() {
		@Override
		public int compare(Gene arg0, Gene arg1) {
			int a = 0;
			double delta = arg0.getRepeats()-arg1.getRepeats();
			if (delta < 0) {
				a = -1;
			} else {
				if (delta > 0) {
					a = 1;
				}
			}
			return a;
		}
	};
	
	//TODO Sort by ascending EEI.
	public static void sortup(Organism organism) {
		switch (organism.getIndex()) {
			case 0: Collections.sort(organism.getGenes(), eei1Up); break;
			case 1: Collections.sort(organism.getGenes(), eei2Up); break;
			case 2: Collections.sort(organism.getGenes(), eei3Up); break;
			case 3: Collections.sort(organism.getGenes(), eei4Up); break;
			case 4: Collections.sort(organism.getGenes(), eei5Up); break;
		}
	}
	
	//Sort by descending EEI.
	public static void sortdown(Organism organism) {
		switch (organism.getIndex()) {
			case 0: Collections.sort(organism.getGenes(), Collections.reverseOrder(eei1Up)); break;
			case 1: Collections.sort(organism.getGenes(), Collections.reverseOrder(eei2Up)); break;
			case 2: Collections.sort(organism.getGenes(), Collections.reverseOrder(eei3Up)); break;
			case 3: Collections.sort(organism.getGenes(), Collections.reverseOrder(eei4Up)); break;
			case 4: Collections.sort(organism.getGenes(), Collections.reverseOrder(eei5Up)); break;
		}
	}
	
	//Sort by ascending EEI.
	public static void sortupEEIN(Organism organism, int n) {
		switch (n) {
			case 1: Collections.sort(organism.getGenes(), eei1Up); break;
			case 2: Collections.sort(organism.getGenes(), eei2Up); break;
			case 3: Collections.sort(organism.getGenes(), eei3Up); break;
			case 4: Collections.sort(organism.getGenes(), eei4Up); break;
			case 5: Collections.sort(organism.getGenes(), eei5Up); break;
			default: System.out.println("Wrong parameter n = "+n+"! It must be {1, 2, 3, 4, 5}.");
		}
	}
	
	//Sort by descending EEI.
	public static void sortdownEEIN(Organism organism, int n) {
		switch (n) {
			case 1: Collections.sort(organism.getGenes(), Collections.reverseOrder(eei1Up)); break;
			case 2: Collections.sort(organism.getGenes(), Collections.reverseOrder(eei2Up)); break;
			case 3: Collections.sort(organism.getGenes(), Collections.reverseOrder(eei3Up)); break;
			case 4: Collections.sort(organism.getGenes(), Collections.reverseOrder(eei4Up)); break;
			case 5: Collections.sort(organism.getGenes(), Collections.reverseOrder(eei5Up)); break;
			default: System.out.println("Wrong parameter n = "+n+"! It must be {1, 2, 3, 4, 5}.");
		}
	}
	
	//Sort by descending EEI.
	public static void sorteeidown(Organism organism) {
		Collections.sort(organism.getGenes(), Collections.reverseOrder(eeiUp));
	}
	
	//Sort the list by ascending gene number.
	public static void sortNup(Organism organism) {
		Collections.sort(organism.getGenes(), geneNumberUp);
	}
	
	//Sort by descending gene number.
	public static void sortNdown(Organism organism) {
		Collections.sort(organism.getGenes(), Collections.reverseOrder(geneNumberUp));
	}
	
	//Sort by ascending LCI1.
	public static void sortLCI1up(Organism organism) {
		Collections.sort(organism.getGenes(), lci1Up);
	}
	
	//Sort by descending LCI1.
	public static void sortLCI1down(Organism organism) {
		Collections.sort(organism.getGenes(), Collections.reverseOrder(lci1Up));
	}
	
	//Sort by ascending LCI2.
	public static void sortLCI2up(Organism organism) {
		Collections.sort(organism.getGenes(), lci2Up);
	}
	
	//Sort by descending LCI2.
	public static void sortLCI2down(Organism organism) {
		Collections.sort(organism.getGenes(), Collections.reverseOrder(lci2Up));
	}
	
	//Sort by ascending Tai.
	public static void sortTaup(Organism organism) {
		Collections.sort(organism.getGenes(), taiUp);
	}
	
	//Sort by descending Tai.
	public static void sortTadown(Organism organism) {
		Collections.sort(organism.getGenes(), Collections.reverseOrder(taiUp));
	}
	
	//Sort by ascending Tei.
	public static void sortTeup(Organism organism) {
		Collections.sort(organism.getGenes(), teiUp);
	}
	
	//Sort by descending Tei.
	public static void sortTedown(Organism organism) {
		Collections.sort(organism.getGenes(), Collections.reverseOrder(teiUp));
	}
	
	//Sort by ascending p_i1.
	public static void sortP_i1up(Organism organism) {
		Collections.sort(organism.getGenes(), pi1Up);
	}
	
	//Sort by ascending p_i2.
	public static void sortP_i2up(Organism organism) {
		Collections.sort(organism.getGenes(), pi2Up);
	}
	
	//Sorting the list for increasing the number of repeats in the gene
	public static void sortRepeatsup(Organism organism) {
		Collections.sort(organism.getGenes(), repeatsUp);
	}
}
