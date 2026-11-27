package org.eloe.other;

import java.util.Collections;
import java.util.Comparator;

import org.eloe.organism.Gene;
import org.eloe.organism.Organism;

public class Sorter {
	//Компаратор для eei
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
	
	//Компаратор для eei1
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
	
	//Компаратор для eei2
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
	
	//Компаратор для eei3
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
	
	//Компаратор для eei4
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
	
	//Компаратор для eei5
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
	
	//Компаратор для номера гена
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
	
	//Компаратор для lci1
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
	
	//Компаратор для lci2
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
	
	//Компаратор для tai
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
	
	//Компаратор для tei
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
	
	//Компаратор для p_i1
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
	
	//Компаратор для p_i2
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
	
	//Компаратор для количества повторов
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
	
	//TODO Сортировка списка по увеличению EEI гена
	public static void sortup(Organism organism) {
		switch (organism.getIndex()) {
			case 0: Collections.sort(organism.getGenes(), eei1Up); break;
			case 1: Collections.sort(organism.getGenes(), eei2Up); break;
			case 2: Collections.sort(organism.getGenes(), eei3Up); break;
			case 3: Collections.sort(organism.getGenes(), eei4Up); break;
			case 4: Collections.sort(organism.getGenes(), eei5Up); break;
		}
	}
	
	//Сортировка списка по уменьшению EEI гена
	public static void sortdown(Organism organism) {
		switch (organism.getIndex()) {
			case 0: Collections.sort(organism.getGenes(), Collections.reverseOrder(eei1Up)); break;
			case 1: Collections.sort(organism.getGenes(), Collections.reverseOrder(eei2Up)); break;
			case 2: Collections.sort(organism.getGenes(), Collections.reverseOrder(eei3Up)); break;
			case 3: Collections.sort(organism.getGenes(), Collections.reverseOrder(eei4Up)); break;
			case 4: Collections.sort(organism.getGenes(), Collections.reverseOrder(eei5Up)); break;
		}
	}
	
	//Сортировка списка по увеличению EEI гена
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
	
	//Сортировка списка по уменьшению EEI гена
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
	
	//Сортировка списка по уменьшению EEI гена
	public static void sorteeidown(Organism organism) {
		Collections.sort(organism.getGenes(), Collections.reverseOrder(eeiUp));
	}
	
	//Сортировка списка по увеличению номера гена
	public static void sortNup(Organism organism) {
		Collections.sort(organism.getGenes(), geneNumberUp);
	}
	
	//Сортировка списка по уменьшению номера гена
	public static void sortNdown(Organism organism) {
		Collections.sort(organism.getGenes(), Collections.reverseOrder(geneNumberUp));
	}
	
	//Сортировка списка по увеличению LCI1 гена
	public static void sortLCI1up(Organism organism) {
		Collections.sort(organism.getGenes(), lci1Up);
	}
	
	//Сортировка списка по уменьшению LCI1 гена
	public static void sortLCI1down(Organism organism) {
		Collections.sort(organism.getGenes(), Collections.reverseOrder(lci1Up));
	}
	
	//Сортировка списка по увеличению LCI2 гена
	public static void sortLCI2up(Organism organism) {
		Collections.sort(organism.getGenes(), lci2Up);
	}
	
	//Сортировка списка по уменьшению LCI2 гена
	public static void sortLCI2down(Organism organism) {
		Collections.sort(organism.getGenes(), Collections.reverseOrder(lci2Up));
	}
	
	//Сортировка списка по увеличению Tai гена
	public static void sortTaup(Organism organism) {
		Collections.sort(organism.getGenes(), taiUp);
	}
	
	//Сортировка списка по уменьшению Tai гена
	public static void sortTadown(Organism organism) {
		Collections.sort(organism.getGenes(), Collections.reverseOrder(taiUp));
	}
	
	//Сортировка списка по увеличению Tei гена
	public static void sortTeup(Organism organism) {
		Collections.sort(organism.getGenes(), teiUp);
	}
	
	//Сортировка списка по уменьшению Tei гена
	public static void sortTedown(Organism organism) {
		Collections.sort(organism.getGenes(), Collections.reverseOrder(teiUp));
	}
	
	//Сортировка списка по увеличению p_i1 гена
	public static void sortP_i1up(Organism organism) {
		Collections.sort(organism.getGenes(), pi1Up);
	}
	
	//Сортировка списка по увеличению p_i2 гена
	public static void sortP_i2up(Organism organism) {
		Collections.sort(organism.getGenes(), pi2Up);
	}
	
	//Сортировка списка по увеличению колическва повторов в гене
	public static void sortRepeatsup(Organism organism) {
		Collections.sort(organism.getGenes(), repeatsUp);
	}
}
