# Dehydration Pathways for Glucose and Cellobiose During Fast Pyrolysis

Mckay W. Easton,\* John J. Nash, and Hilkka I. Kenttãmaa

Department of Chemistry, Purdue University, West Lafayette, Indiana 47907, United States

Supporting Information

ABSTRACT: A full understanding of all possible elementary reactions applicable to cellulose fast pyrolysis is key to developing a comprehensive kinetic model for fast pyrolysis of cellulose. Since water is an observed product of fast pyrolysis of cellulose, the energetics of the dehydration reactions of cellulose were explored computationally by using density functional theory. Glucose and cellobiose were selected as the cellulose model compounds. The four water loss mechanisms studied are Maccoll elimination, Pinacol ring contraction, cyclic Grob fragmentation, and alcohol condensation, some of which have not been considered previously in the literature. Levoglucosan formation via alcohol condensation has the

lowest calculated free-energy barrier $(50.4\mathrm{kcalmol}^{-1})$ for glucose dehydration. All other water loss reactions have calculated free-energy barriers greater than $60~\mathrm{kcal~mol}^{-1}$ . Cellobiose dehydration shows similar trends to those of glucose, suggesting that these reactions are applicable to glucooligosaccharides with higher degrees of polymerization. Secondary reactions of dehydrated glucose and dehydrated cellobiose via retro-Diels-Alder and aldol rearrangement mechanisms are also explored computationally.

![](images/89ad36ef6f1bef7acfbe5b22605be1fbd973ae47a65a220394aac7e00a5017ae.jpg)

# INTRODUCTION

Recent energy concerns have fueled research into alternative renewable energy sources, including biomass. The conversion of lignocellulosic biomass by fast pyrolysis followed by catalytic hydrodeoxygenation is an attractive option for the generation of fuels from biomass while utilizing the existing liquid fuel infrastructure. In addition to fuel production, fast pyrolysis of biomass can subsidize the current biomass refining industry for the production of fine chemicals by downstream processing of the pyrolysis vapors. Pyrolysis varies from other biorefining techniques that preserve the cellulosic fibers for use in industries such as the pulp and paper industry, since it depolymerizes the biomass into smaller molecules. The $\mathrm{H}_2\mathrm{Bioil}$ process introduced recently in the literature utilizes pyrolysis with very fast heating to achieve higher yields of liquid products and less char formation compared to standard pyrolysis heating rates. This approach has nearly twice the carbon recovery potential of cellulosic ethanol production by fermentation, with a theoretical carbon efficiency of $\sim 70\%$ . Carbon efficiencies greater than $50\%$ have already been achieved.

Cellulose is the most abundant biopolymer in biomass, making up between $40\%$ and $60\%$ of the total biomass by weight, depending on the feedstock. Therefore, understanding of the fast pyrolysis process of cellulose is of great interest. While early mechanistic models of cellulose depolymerization during fast pyrolysis used lumped kinetic parameters and unknown "active cellulose" intermediates, more recent kinetic

models utilize the current understanding of the elementary reactions that may occur during pyrolysis. $^{10}$ The elementary reactions and specific mechanisms of cellulose fast pyrolysis are under active debate. $^{2,10-16}$ In recent experiments, glucose (GCL) and cellobiose were considered constituents of the "active cellulose" intermediate of cellulose pyrolysis, which is likely to undergo dehydration. $^{17}$ Water as a product from pyrolysis accounts for $41\mathrm{wt}\%$ for glucose, $^{18}$ $24\mathrm{wt}\%$ for cellobiose, $^{18}$ and up to $\sim 40\mathrm{wt}\%$ for cellulose, $^{4,18,19}$ indicating that dehydration must occur. Furthermore, water loss reactions have been linked to char formation at temperatures above $200^{\circ}\mathrm{C}$ , of which intramolecular dehydration reactions were considered intrinsic to formation of carbonyl, double-bonded, and cyclic ether structures. $^{20}$

Several water loss mechanisms for glucosaccharides have been reported in the literature, including 1,2-dehydration,[10,17,19,21,22] Pinacol rearrangement,[19,23] "Tandem Alkaline Pinacol Rearrangement/Retro Aldol Fragmentation" (TAPRRAF),[24] cyclic Grob fragmentation (CGF),[19,24] and generalized alcohol condensation.[25] These mechanisms have been used to explain the cleavages of C-C and C-O bonds that lead to carbohydrate fragments such as formic acid, formaldehyde, glycolaldehyde, and erythrose, as well as many other oxygenated hydrocarbons. Formation of some fast pyrolysis

Received: March 8, 2018

Revised: September 14, 2018

Published: September 14, 2018

![](images/18a5c70aa76357f0950ac9e2ae8f52d3979565d2efb5ea53f4f28e3e404985bf.jpg)

![](images/50708c5aac22125280c9c32cb4710ad60f64d3196991f14a0f013d5c8c250191.jpg)

![](images/734e589675eb8db8a12248362f7804f1968f9581787facbd228ec59d497e40be.jpg)  
Figure 1. Numbering scheme for glucose and cellobiose. For cellobiose, atoms on the nonreducing end are indicated by primed numbers.

products has also been explained by a retro-aldol condensation reaction mechanism $^{26}$ and a retro-ene mechanism. $^{24}$ The retro-Diels-Alder mechanism has been suggested to be of importance during cellulose pyrolysis. $^{10,23,27}$ Some computational work has been done on select reactions of cellobiose, $^{11,12,14,28,29}$ as well as glucose, $^{19,21,25,30}$ but a comprehensive survey on the dehydration and retro-Diels-Alder decomposition of cellobiose is notably absent in the literature. We report here a computational study on Maccoll elimination, $^{31}$ Pinacol rearrangement, $^{32}$ cyclic Grob fragmentation, $^{33}$ condensation cyclization, retro-aldol condensation, and retro-Diels-Alder reactions of glucose and cellobiose. Although intermolecular interactions may offer a catalyzing or inhibiting effect for many of these reactions, $^{25,34}$ this study focuses solely on intramolecular reactions as a first approximation to the complex reaction network of fast pyrolysis.

# COMPUTATIONAL METHODS

Conformational analysis was accomplished with the Macro-Model program as part of the Schrodinger software package $^{35}$ by using $\mathrm{MCMM}^{36}$ torsional sampling with the MMFF94s force field. $^{37}$ More accurate energies for the low-energy conformers were obtained by geometry optimization and energy calculations using density functional theory (DFT) with the Gaussian 09 program suite. $^{38}$ Note that numbering of the atoms in glucose and cellobiose (Figure 1) is according to convention.

It is important to select an appropriate conformation of cellobiose for calculations, since pyrolysis is an inherently multiphase process, and cellobiose can adopt different conformations depending on its phase. Extensive conformational studies on the rotamers of cellobiose have shown several low-energy structures. In this study, conformational analysis along with DFT optimization of cellobiose revealed two classes of relevant conformations: (1) the syn-form, with the hydroxymethylene groups on opposite sides of the rings relative to each other, which is the conformation most closely resembling the crystalline cellulose form I $\alpha$ present in native plant matter, and (2) the anti-form, where the hydroxymethylene groups are on the same sides of the rings relative to each other. On the one hand, crystallographic data on solid cellobiose indicate that, in the solid state, the molecule is in the syn-form. On the other hand, previous work has shown that the energetically most favorable conformation of cellobiose in the gas phase is the anti form. Conformational changes may occur when the particles enter the melt phase of pyrolysis. However, liquid-phase steered molecular dynamics (SMD) simulations with a continuum solvation model

indicate that in the liquid phase, the crystal-like syn conformation is lower in energy than the gas-phase conformation.40 Furthermore, conformational changes in cellulose become increasingly more difficult as the number of neighboring cellulose chains increases due to intermolecular forces between chains.34 For these reasons, the syn conformation of cellobiose was selected as the conformation that most closely reflects the structure of cellulose during pyrolysis. Transition-state geometries were constructed to mimic the syn conformation of cellobiose as the starting structure for the DFT calculations.

$\mathrm{B3LYP}^{47-50}$ has been used often in the literature to accurately predict molecular structures,[51] but in many cases, this functional underestimates the energies of transition-state structures.[52,53] M06-2X is among the best functionals for calculating barrier heights, particularly for unimolecular reactions and hydrogen atom transfer reactions.[54] M06-2X has been shown to accurately replicate results obtained using MP4(SDQ), and it has been used extensively in the literature for studying sugar chemistry.[11] The 6-311++G(d,p) basis set was used for all DFT calculations, as suggested in the literature.[55] Thus, unless otherwise noted, all DFT calculations were performed at the M06-2X/6-311++G(d,p) level of theory. All geometry optimizations were followed by a full Hessian calculation to guarantee that the minima had all positive eigenvalues and the transition states had exactly one negative eigenvalue. The coordinates of minima and transition-state structures can be found in the Supporting Information. Intrinsic reaction coordinate (IRC) calculations were used to confirm that reactant geometries associated with each transition state were in the syn conformation. Product geometries associated with each transition state were obtained from the IRC for further optimization. Thermal corrections for internal energy, enthalpy, and Gibbs free energy were calculated using statistical mechanics, assuming an ideal gas and using the harmonic oscillator approximation at a pressure of 1 atm. Where applicable, partial atomic charges were computed by fitting point charges to the electrostatic potential using the CHELPG algorithm.[56]

# RESULTS AND DISCUSSION

Maccoll Elimination. Maccoll elimination involves cleavage of a $\mathrm{C - O}$ bond with a concerted hydrogen atom transfer from a $\beta$ -carbon to the oxygen atom to form a $\mathrm{C - C}$ double bond and eliminate either $\mathrm{H}_2\mathrm{O}$ or a compound containing a hydroxyl group.[31,57] When referring to a Maccoll elimination, two numbers will be used to describe the carbon atoms involved in the $\mathrm{C - O}$ and $\mathrm{C - H}$ bond cleavages, respectively (Figure 2). Glucose can undergo Maccoll elimination via eight

![](images/22c04d0642cfa6114d519b50827dfcca5b85f055987f75d051cebceca5748793.jpg)  
Figure 2. Example of Maccoll elimination for GLC, depicting the "GLC32Mac" reaction (i.e., "3" indicates the carbon atom at which the C-O bond is broken, and "2" indicates the carbon atom at which the C-H bond is broken).

different reactions (Table 1), while cellobiose has 16 possible reactions (Table 2). Each of the Maccoll eliminations of glucose (or cellobiose) leads to an alkene, one (two) of which is (are) terminal alkenes and 7 (14) of which are cyclic alkenes that can undergo a retro-Diels-Alder reaction.

Tables 1 and 2 contain the reaction product and the calculated free-energy barriers and reaction free energies for glucose and cellobiose at room temperature and at a typical fast pyrolysis temperature $(600^{\circ}\mathrm{C})$ . In most cases, the increase in temperature causes a slight decrease in the barrier height $(\sim 2\mathrm{kcalmol}^{-1})$ and a large decrease in the reaction energy $(\sim 20\mathrm{kcalmol}^{-1})$ . The large change in reaction energy is due to the additional rotational and translational entropy resulting from decomposing $1\mathrm{mol}$ of reactant into $2\mathrm{mol}$ of product. The transition-state structures possess the same number of translational and rotational degrees of freedom as the reactants, and only the vibrational entropy plays a minor but detectable role in the temperature variance in the Gibbs free energy of the transition-state structures. In two cases (CB12Mac and CB34Mac), the relative transition-state energy increases at the elevated temperature. In both of these cases, the leaving water molecule acts as a hydrogen-bond donor, which constrains some of the normal modes of vibration in the molecule. Although the effect is small, the constrained geometry causes an elevation in the vibrational temperatures, leading to less populated vibrational states, and an overall reduction in the vibrational entropy relative to other transition

states. This accounts for the slight increase in barrier height for these two reactions at higher temperatures.

The lowest free-energy barrier among the various Maccoll elimination reactions of glucose is $65.9\mathrm{kcalmol}^{-1}$ calculated for GLC12Mac (Table 1), which involves a loss of a hydroxyl group from C1 and a loss of a hydrogen atom from C2. The products from this reaction are water and dehydrated glucose having a double bond between carbons C1 and C2 (Table 1), which may be an intermediate to levoglucosan or levoglucosenone formation.[28] The anomeric carbon benefits from conjugation with the electron-rich oxygen atom in the ring, leading to a partial charge of $0.480\mathrm{C}$ in the transition-state structure of GLC12Mac, which is similar in structure to an oxocarbenium. The delocalization of charge likely plays a role in the stabilization of this transition state. By comparison, the reaction with the largest barrier is GLC21Mac, where the hydroxyl group is lost from C2, while the hydrogen atom is lost from C1. The carbon atom that undergoes $\mathrm{C - O}$ cleavage (C2) in reaction GLC21Mac has a substantially less delocalized partial charge of $0.724\mathrm{C}$ , while the carbon atom that undergoes $\mathrm{C - H}$ cleavage (C1) has a partial charge of $-0.217\mathrm{C}$ , which causes destabilization in the $\mathrm{C - O}$ bond between C1 and the ring oxygen. This electronic effect is likely the reason for the $14.5\mathrm{kcalmol}^{-1}$ difference in barrier heights between reaction GLC12Mac and reaction GLC21Mac.

Seven of the eight Maccoll elimination reactions for glucose result in a cyclic alkene. The GLC65Mac reaction is unique, because it forms a terminal alkene at the C5 position. It has the second-lowest barrier at $70.1\mathrm{kcalmol}^{-1}$ . The remaining Maccoll elimination reactions involve the loss of a hydroxyl group at C2, C3, or C4. The barrier heights for these five reactions are within $2.7\mathrm{kcalmol}^{-1}$ of each other, suggesting that these reactions are likely to be competitive. The barrier heights reported here coincide with other literature reports on 1,2-dehydration of glucose. For example, this report claims a Gibbs free energy of $65.9\mathrm{kcalmol}^{-1}$ for the transition state of reaction GLC12Mac, while Mayes et al. report a Gibbs free energy of $59.6\mathrm{kcalmol}^{-1}$ at a temperature of $500^{\circ}\mathrm{C}$ using the

Table 1. Calculated Transition-State $\left( {\Delta {G}^{ \ddagger  }}\right)$ and Reaction $\left( {\Delta {G}_{\mathrm{{rxn}}}}\right)$ Gibbs Free Energies (kcal mol ${}^{-1}$ ) for Maccoll Eliminations of Equatorial $\beta$ -D-Glucopyranose (Glucose) ${}^{a}$   

<table><tr><td rowspan="2">Reaction Label</td><td rowspan="2">Product Structure</td><td colspan="2">25°C</td><td colspan="2">600°C</td></tr><tr><td>ΔG‡</td><td>ΔGrxn</td><td>ΔG‡</td><td>ΔGrxn</td></tr><tr><td>GLC12Mac</td><td>HO HO OH</td><td>69.4</td><td>5.7</td><td>65.9</td><td>-17.9</td></tr><tr><td>GLC21Mac</td><td>HO HO OH</td><td>81.9</td><td>-0.7</td><td>80.4</td><td>-25.3</td></tr><tr><td>GLC23Mac</td><td>HO HO OH</td><td>74.6</td><td>1.1</td><td>72.3</td><td>-22.6</td></tr><tr><td>GLC32Mac</td><td>HO HO OH</td><td>77.4</td><td>4.5</td><td>73.6</td><td>-20.1</td></tr><tr><td>GLC34Mac</td><td>HO HO OH</td><td>75.9</td><td>4.4</td><td>74.3</td><td>-20.7</td></tr><tr><td>GLC43Mac</td><td>HO HO OH</td><td>75.5</td><td>-2.5</td><td>73.8</td><td>-25.7</td></tr><tr><td>GLC45Mac</td><td>HO HO OH</td><td>78.1</td><td>1.3</td><td>75.0</td><td>-23.3</td></tr><tr><td>GLC65Mac</td><td>HO HO OH</td><td>71.6</td><td>1.2</td><td>70.1</td><td>-22.2</td></tr></table>

${}^{a}$ The reaction labels indicate the reactant (GLC = glucose),the respective carbon atoms from which $\mathrm{{OH}}$ and $\mathrm{H}$ are lost during dehydration (numbering according to Figure 1),and the reaction type (Mac = Maccoll). Note that, in each case, water (not shown) is also a product.

Table 2. Calculated Transition-State $\left( {\Delta {G}^{ \ddagger  }}\right)$ and Reaction $\left( {\Delta {G}_{\mathrm{{rxn}}}}\right)$ Gibbs Free Energies (kcal mol ${}^{-1}$ ) for Maccoll Eliminations of Cellobiose ${}^{a}$   

<table><tr><td rowspan="2">Reaction Label</td><td rowspan="2">Product Structure(s)</td><td colspan="2">25°C</td><td colspan="2">600°C</td></tr><tr><td>ΔG‡</td><td>ΔGrxn</td><td>ΔG‡</td><td>ΔGrxn</td></tr><tr><td>CB12Mac</td><td></td><td>63.3</td><td>9.2</td><td>63.6</td><td>-15.3</td></tr><tr><td>CB21Mac</td><td></td><td>82.9</td><td>-1.3</td><td>80.6</td><td>-24.7</td></tr><tr><td>CB23Mac</td><td></td><td>72.7</td><td>-1.1</td><td>71.3</td><td>-25.8</td></tr><tr><td>CB32Mac</td><td></td><td>72.4</td><td>1.8</td><td>70.8</td><td>-23.6</td></tr><tr><td>CB34Mac</td><td></td><td>64.4</td><td>3.8</td><td>65.0</td><td>-23.5</td></tr><tr><td>CB43Mac</td><td></td><td>68.2</td><td>2.2</td><td>65.9</td><td>-24.1</td></tr><tr><td>CB45Mac</td><td></td><td>71.9</td><td>6.0</td><td>68.7</td><td>-21.8</td></tr><tr><td>CB65Mac</td><td></td><td>69.3</td><td>1.2</td><td>68.6</td><td>-23.1</td></tr><tr><td>CB1&#x27;2&#x27;Mac</td><td></td><td>66.9</td><td>10.4</td><td>62.5</td><td>-16.4</td></tr><tr><td>CB2&#x27;1&#x27;Mac</td><td></td><td>82.5</td><td>4.0</td><td>80.6</td><td>-20.5</td></tr><tr><td>CB2&#x27;3&#x27;Mac</td><td></td><td>75.9</td><td>2.1</td><td>73.4</td><td>-22.6</td></tr><tr><td>CB3&#x27;2&#x27;Mac</td><td></td><td>75.2</td><td>3.3</td><td>71.1</td><td>-22.2</td></tr><tr><td>CB3&#x27;4&#x27;Mac</td><td></td><td>74.6</td><td>4.7</td><td>72.7</td><td>-19.9</td></tr><tr><td>CB4&#x27;3&#x27;Mac</td><td></td><td>73.3</td><td>-0.3</td><td>70.7</td><td>-24.4</td></tr><tr><td>CB4&#x27;5&#x27;Mac</td><td></td><td>71.7</td><td>-0.8</td><td>66.8</td><td>-26.0</td></tr><tr><td>CB6&#x27;5&#x27;Mac</td><td></td><td>68.8</td><td>-0.1</td><td>67.2</td><td>-24.5</td></tr></table>

${}^{a}$ The reaction labels indicate the reactant $\left( {\mathrm{{CB}} = \text{cellobiose}}\right)$ ,the respective carbon atoms from which $\mathrm{{OH}}$ and $\mathrm{H}$ are lost during dehydration (numbering according to Figure 1),and the reaction type (Mac = Maccoll). Note that, in each case, water (not shown) is also a product.

M06-2X/6-311+G(2df,p) level of theory with an implicit solvent model using tetrahydrofuran (THF) as the solvent. $^{19}$ The major discrepancies are attributed to the additional hydrogen bonding provided by the gauche-gauche hydroxymethylene conformer that is used in literature as opposed to the gauche-trans conformer that is used in this report, as well as the use of an implicit solvent model in literature. Assary et al. report an activation energy of 60.5 kcal mol $^{-1}$ at a temperature of 25 °C using an implicit solvent model with the G4 level of theory. $^{30}$ The differences in energies between the current report and the report by Assary et al. must be a result of the difference in the levels of theory—including the use of solvation—since the conformations used ( ${}^{4}\mathrm{C}_{1}$ ring puckered gauche-trans rotamer) are identical across both reports. Similarly, this report claims a Gibbs free energy of 80.4 kcal mol $^{-1}$ for the transition state of reaction GLC21Mac, while

Mayes et al. report a Gibbs free energy of $79.6\mathrm{kcalmol}^{-1}$ , which also might be lower due to the difference in conformation at the hydroxymethylene group as well as the difference in implicit solvation (notably Assary did not consider this higher-energy reaction). All other available barrier heights show similar trends to those in Table 1, with the exception that the calculations performed by Assary with the solvent model are lower by $\sim 10\mathrm{kcalmol}^{-1}$ . Calculations on the 1,2-dehydration of cellobiose are currently absent in literature.

Comparison of the Maccoll eliminations for glucose and for cellobiose reveals some similarities in the trends for the barrier heights of analogous dehydration reactions (Table 2). For example, similar to glucose, the two lowest calculated barriers for cellobiose involve cleavage of a C-O bond at the anomeric carbon, slightly favoring cleavage of the C-O bond in the

![](images/4c9410e823d0774236c83a1f6338d2f1da49c4dcdd664749e7b8005ab23910ed.jpg)

![](images/94e311511750385fbbed088c06f77469432811bde1410cc7c53fce504b413c7f.jpg)  
Figure 3. Comparison of free-energy barriers $(\Delta G^{\ddagger})$ for Maccoll eliminations of GLC and each glucose ring in cellobiose (reducing end: "CB-R"; nonreducing end "CB-NR").   
Figure 4. Example of Pinacol ring contraction for GLC, depicting the "GLC12Pin" reaction (i.e., "1" and "2" indicate the carbon atoms at which the C-O bond is broken and a carbonyl group is formed, respectively).

nonreducing end (CB1'2'Mac; 62.5 kcal mol $^{-1}$ ; Table 2) over the reducing end of cellobiose (CB12Mac; 63.3 kcal mol $^{-1}$ ; Table 2). It is important to note that in the case of the nonreducing end, the reaction corresponds to a Maccoll elimination of the intact reducing end glucose ring by cleavage of the glycosidic bond.

There are three cases where Maccoll elimination of cellobiose leads to formation of a glucose molecule instead of a water molecule: CB43Mac, CB45Mac, and CB1'2'Mac. Although these reactions do not lead to elimination of water, they are still relevant for carbohydrate fast pyrolysis, since the glycosidic linkage is broken. On the one hand, the free-energy

Scheme 1. Mechanism for Water Loss and Ring Contraction of Glucose Protonated at O2

![](images/3c10407901086cd44566f6d42ba1aebb96f75a8d521212ce0b909c32b2ee49ed.jpg)

barriers for these reactions range between 62.5 and $68.7\mathrm{kcal}$ $\mathrm{mol}^{-1}$ , and they are similar in magnitude to other reports.[28] On the other hand, the bond dissociation energies (BDE) for homolysis and heterolysis of the glycosidic bond have been calculated to be 79.1 and $157.5\mathrm{kcal mol}^{-1}$ , respectively.[58] Hence, cleavage of the glycosidic bond is substantially more favorable when concerted with a hydrogen atom transfer as in Maccoll elimination. Cleavage of the aglyconic bond of cellobiose during the Maccoll reactions CB43Mac and CB45Mac leads to intact glucose from the nonreducing end of cellobiose. In contrast, the CB1'2'Mac reaction involves glycosidic bond cleavage, producing intact glucose from the reducing end. In our previous work,[14] glucose was demonstrated to be formed exclusively from the nonreducing end upon fast pyrolysis of cellobiose. On the basis of the reaction barriers calculated here, it would be expected that some glucose from the reducing end would be produced in fast pyrolysis experiments (via CB1'2'Mac; Table 2) if these were the dominant reactions for glucose formation. Thus, these reactions alone cannot fully explain glucose formation from pyrolysis of cellobiose.

The Maccoll elimination reactions occurring on the reducing and nonreducing ends of cellobiose show similar trends as glucose, but they tend to have slightly lower free-energy barriers than Maccoll elimination reactions of glucose (Figure 3). The only exceptions are GLC21Mac with a barrier that is $0.2\mathrm{kcalmol}^{-1}$ lower than the barriers of both CB21Mac and $\mathrm{CB2^{\prime}1^{\prime}Mac}$ , as well as GLC23Mac with a barrier that is 1.2 kcal mol $^{-1}$ lower than that of CB2'3'Mac. The decreases in transition-state energy can be attributed to the additional

Table 3. Calculated Transition State $\left( {\Delta {G}^{ \ddagger  }}\right)$ and Reaction $\left( {\Delta {G}_{\mathrm{{rxn}}}}\right)$ Gibbs Free Energies (kcal mol ${}^{-1}$ ) for Pinacol Ring Contractions of Equatorial $\beta$ -D-Glucopyranose (Glucose) ${}^{a}$   

<table><tr><td rowspan="2">Reaction Label</td><td rowspan="2">Product Structure</td><td colspan="2">25°C</td><td colspan="2">600°C</td></tr><tr><td>ΔG‡</td><td>ΔGrxn</td><td>ΔG‡</td><td>ΔGrxn</td></tr><tr><td>GLC12Pin</td><td>HO HO OH</td><td>72.9</td><td>3.9</td><td>71.1</td><td>-21.6</td></tr><tr><td>GLC21Pin</td><td>HO HO OH</td><td>68.7</td><td>2.5</td><td>67.1</td><td>-24.0</td></tr><tr><td>GLC23Pin</td><td>HO HO OH</td><td>74.2</td><td>-3.4</td><td>72.3</td><td>-27.6</td></tr><tr><td>GLC32Pin</td><td>HO HO OH</td><td>71.6</td><td>-5.1</td><td>70.5</td><td>-30.0</td></tr><tr><td>GLC34Pin</td><td>HO HO OH</td><td>72.0</td><td>-2.2</td><td>69.8</td><td>-27.5</td></tr><tr><td>GLC43Pin</td><td>HO HO OH</td><td>76.1</td><td>-2.0</td><td>74.5</td><td>-27.7</td></tr></table>

${}^{a}$ The reaction labels indicate the reactant $\left( {\mathrm{{GLC}} = \text{glucose}}\right)$ ,the respective carbon atoms at which a $\mathrm{C} - \mathrm{O}$ bond is broken and a $\mathrm{C} = \mathrm{O}$ bond is formed (numbering according to Figure 1),and the reaction type (Pin = Pinacol rearrangement). Note that, in each case, water (not shown) is also a product.

Table 4. Calculated Transition-State $\left( {\Delta {G}^{ \ddagger  }}\right)$ and Reaction $\left( {\Delta {G}_{\mathrm{{rxn}}}}\right)$ Gibbs Free Energies (kcal mol ${}^{-1}$ ) for Pinacol Ring Contractions of Cellobiose ${}^{a}$   

<table><tr><td rowspan="2">Reaction Label</td><td rowspan="2">Product Structure(s)</td><td colspan="2">25°C</td><td colspan="2">600°C</td></tr><tr><td>ΔG‡</td><td>ΔGrxn</td><td>ΔG‡</td><td>ΔGrxn</td></tr><tr><td>CB12Pin</td><td>HO HO HO HO HO OH</td><td>62.0</td><td>4.0</td><td>61.5</td><td>-20.6</td></tr><tr><td>CB21Pin</td><td>HO HO HO HO HO OH</td><td>70.1</td><td>7.1</td><td>68.2</td><td>-19.6</td></tr><tr><td>CB23Pin</td><td>HO HO HO HO HO OH</td><td>77.5</td><td>-0.2</td><td>75.2</td><td>-26.0</td></tr><tr><td>CB32Pin</td><td>HO HO HO HO HO OH</td><td>71.5</td><td>-2.9</td><td>72.2</td><td>-26.9</td></tr><tr><td>CB43Pin</td><td>HO HO HO HO HO OH</td><td>77.5</td><td>2.7</td><td>76.8</td><td>-28.1</td></tr><tr><td>CB1&#x27;2&#x27;Pin</td><td>HO HO HO HO HO OH</td><td>76.3</td><td>8.6</td><td>72.6</td><td>-22.0</td></tr><tr><td>CB2&#x27;3&#x27;Pin</td><td>HO HO HO HO HO OH</td><td>76.9</td><td>-1.5</td><td>74.1</td><td>-25.1</td></tr><tr><td>CB3&#x27;2&#x27;Pin</td><td>HO HO HO HO HO OH</td><td>78.1</td><td>-2.2</td><td>73.7</td><td>-25.7</td></tr><tr><td>CB3&#x27;4&#x27;Pin</td><td>HO HO HO HO HO OH</td><td>74.9</td><td>2.6</td><td>71.8</td><td>-25.8</td></tr><tr><td>CB4&#x27;3&#x27;Pin</td><td>HO HO HO HO HO OH</td><td>78.1</td><td>-2.7</td><td>75.7</td><td>-25.8</td></tr></table>

${}^{a}$ The reaction labels indicate the reactant $\left( {\mathrm{{CB}} = \text{cellobiose}}\right)$ ,the respective carbon atoms at which a $\mathrm{C} - \mathrm{O}$ bond is broken and a $\mathrm{C} = \mathrm{O}$ bond is formed (numbering according to Figure 1),and the reaction type (Pin = Pinacol ring contraction). Note that, in each case, water (not shown) is also a product.

![](images/d9f30ebd9c5138b9d46b96a0f2e8bb98fadea3fd7d59ae98a03cbac7854e6799.jpg)  
Figure 5. Generalized cyclic Grob fragmentation mechanism.

hydrogen bonding in cellobiose versus glucose. For example, the transition-state structure of GLC12Mac possesses three hydrogen bonds, while the transition-state structure of CB12Mac possesses seven hydrogen bonds. Interestingly, the transition-state structures of the "21" reactions shown in Figure 3 possess fewer hydrogen bonds than the other transition-state structures and show the least deviation between glucose and cellobiose. As discussed previously, the leaving water group in the transition-state structure of CB34Mac contributes additional hydrogen bonding, which may account for the more dramatic difference in barrier heights between analogous "34" reactions of glucose and cellobiose compared to other sets of reactions. These results indicate that the barrier heights for glucose and cellobiose may be adequate preliminary approximations for Maccoll eliminations of larger glucooligosaccharides but that the effect of hydrogen bonding should not be overlooked. Furthermore, we note that these reaction barriers may be lower under certain conditions; for example, in the presence of water,[25] under acidic conditions,[21,22] or in the presence of alkali metal cations.[22]

Pinacol Ring Contraction. While a Maccoll elimination reaction involves cleavage of a C-H bond, a second

dehydration reaction type considered here, termed a "Pinacol ring contraction" due to its similarity to the Pinacol rearrangement reaction,[24,32] instead involves cleavage of an O–H bond. In this mechanism, the hydrogen atom source is the OH group adjacent to the OH group that is eliminated. This reaction forms a carbonyl group while simultaneously contracting the six-membered pyranosyl ring to a five-membered furanosyl ring (Figure 4). The Pinacol ring contraction reactions of glucose produce six distinct products, including three sets of stereoisomers: (1) GLC12Pin and GLC21Pin, (2) GLC23Pin and GLC32Pin, and (3) GLC34Pin and GLC43Pin. The Gibbs free-energy barriers for Pinacol ring contraction of glucose are comparable at pyrolysis temperature and at room temperature, but the overall reactions are more exergonic by $24 - 27\mathrm{kcalmol}^{-1}$ at the higher temperature (Table 3). The reaction with the lowest free-energy barrier is GLC21Pin (Table 3). The product of this reaction, as well as the product from reaction of GLC12Pin, are intermediates to the formation of 5-hydroxymethylfurfural (HMF) via water losses at C3 and C4. Mayes et al. have also calculated the Gibbs free energy for the transition state of the 21Pin reaction from the $\alpha$ -glucopyranose anomer of glucose $(62.4\mathrm{kcalmol}^{-1})$ as a possible route to HMF,[19] but reaction energies for all other Pinacol reactions are not currently in the literature. The GLC21Pin reaction involves cleavage of the C-O bond at the C2 atom similar to the GLC21Mac reaction, which was the highest of all the Maccoll reactions for glucose. The reasoning for this is the localized charge left on C2, which

Table 5. Calculated Transition-State $\left( {\Delta {G}^{ \ddagger  }}\right)$ and Reaction $\left( {\Delta {G}_{\mathrm{{rxn}}}}\right)$ Gibbs Free Energies (kcal mol ${}^{-1}$ ) for Cyclic Grob Fragmentation of $\beta$ -D-Glucopyranose (Glucose) ${}^{a}$   

<table><tr><td rowspan="2">Reaction Label</td><td rowspan="2">Product Structure(s)</td><td colspan="2">25°C</td><td colspan="2">600°C</td></tr><tr><td>ΔG‡</td><td>ΔGrxn</td><td>ΔG‡</td><td>ΔGrxn</td></tr><tr><td>GLCax13Grob</td><td>HO
OH
OH</td><td>63.7</td><td>20.8</td><td>62.7</td><td>-9.9</td></tr><tr><td>GLCax24Grob</td><td>HO
OH
OH</td><td>95.0</td><td>19.0</td><td>92.1</td><td>-12.5</td></tr><tr><td>GLCax31Grob</td><td>HO
OH
OH</td><td>78.5</td><td>3.8</td><td>80.8</td><td>-23.6</td></tr><tr><td>GLCax42Grob</td><td>HO
OH
OH</td><td>81.0</td><td>17.3</td><td>83.6</td><td>-6.8</td></tr><tr><td>GLCax46Grob</td><td>O2CH2OH</td><td>82.3</td><td>13.9</td><td>82.7</td><td>-16.3</td></tr><tr><td>GLCeq46Grob</td><td>O2CH2OH</td><td>83.5</td><td>13.7</td><td>81.4</td><td>-17.7</td></tr><tr><td>GLCax53Grob</td><td>HO
OH
OH</td><td>68.7</td><td>27.8</td><td>68.6</td><td>-18.7</td></tr><tr><td>GLCax64Grob</td><td>HO
OH
OH</td><td>69.1</td><td>7.7</td><td>70.1</td><td>-22.7</td></tr></table>

${}^{a}$ The reaction labels indicate the reactant (GLCeq = equatorial glucose; GLCax = axial glucose), the respective carbon atoms at which a $\mathrm{C} - \mathrm{O}$ bond is broken and a $\mathrm{C} = \mathrm{O}$ bond is formed (numbering according to Figure 1),and the reaction type (Grob = cyclic Grob fragmentation). Note that all free-energy barriers are calculated relative to equatorial glucose, and in the case of GLCax53Grob transformation from an ether to a hydroxyl group occurs instead of $\mathrm{C} - \mathrm{O}$ bond cleavage. Note that, in each case, water (not shown) is also a product.

in the case of GLC21Pin is delocalized to the ring oxygen by the furan cyclization of the Pinacol reaction. This is clear by observing that the partial charge on C2 $(-0.149\mathrm{C})$ is lower than for glucose (0.016 C). All possible Pinacol ring contractions for glucose are shown in Table 3, including the Gibbs free-energy barriers.

The free-energy barriers for Pinacol rearrangement are larger than the most favorable Maccoll elimination reactions. However, the Pinacol rearrangement reaction has been reported to be catalyzed by acids,[32] including formic acid, which is known to be a product of biomass fast pyrolysis.[18] Furthermore, when an OH group in glucose is protonated, $\mathrm{H}_2\mathrm{O}$ will spontaneously cleave from the ring, accompanied by a ring contraction to stabilize the carbocation formed (Scheme 1). This spontaneous water elimination and ring contraction has been observed by others during molecular orbital geometry optimization of protonated glucose.[59]

In the case of cellobiose, the lowest free-energy barriers for Pinacol ring contraction and water elimination correspond to cleavage of the $\mathrm{C - O}$ bond at the anomeric carbon of each glucose ring (Table 4). For reaction CB12Pin, water is lost from the anomeric carbon of the reducing end with a free-energy barrier of $61.5\mathrm{kcalmol}^{-1}$

For the analogous reaction of C-O cleavage on the nonreducing end (i.e., CB1'2'Pin), the products are glucose (instead of water) and the same product as for the GLC12Pin reaction (Table 4), and the barrier is $72.6\mathrm{kcalmol}^{-1}$ . One other Pinacol ring contraction reaction for cellobiose produces glucose. While CB1'2'Pin produces glucose from the reducing end, CB43Pin produces glucose from the nonreducing end

with a barrier of $76.8\mathrm{kcalmol}^{-1}$ (Table 4). These barriers for glucose formation are larger than those for the Maccoll eliminations that produce glucose, suggesting that Pinacol ring contraction reactions are less likely pathways for formation of glucose. As with glucose, the free-energy barriers for all Pinacol ring contraction reactions of cellobiose are relatively independent of temperature, while the free energy of reaction decreases by $23 - 31\mathrm{kcalmol}^{-1}$ between standard temperature and the pyrolysis temperature (Table 4).

Cyclic Grob Fragmentation. The cyclic Grob fragmentation reaction is another dehydration mechanism, which applies to 1,3-diols such as glucose. $^{10,24}$ During the reaction, a C-C bond is broken and an alkene and an aldehyde are formed (Figure 5). This reaction cannot proceed if the two hydroxyl groups are not in close proximity. The $^1\mathrm{C}_4$ "axial" conformation of glucose, which is $6.7~\mathrm{kcal~mol^{-1}}$ higher in energy than the lowest-energy $^4\mathrm{C}_1$ "equatorial" conformation, possesses hydroxyl groups that are in axial positions, which causes the participating hydroxyl groups to be close enough to react. With one exception, the hydroxyl groups of the $^4\mathrm{C}_1$ equatorial conformation are not close enough to each other for reaction to occur. The exception is the cyclic Grob fragmentation of glucose in the $^4\mathrm{C}_1$ equatorial position, where the C-O bond is broken at C4 and a carbonyl forms at C6. This particular reaction is referred to as GLCeq46Grob to differentiate it from GLCax46Grob, which involves the axial conformation of glucose. Note that, in the calculations, the free energy required to "flip" the glucose ring from the more stable equatorial conformation to the axial conformation was ignored, because the barrier for glucose ring-flipping is 12.1 kcal

Table 6. Calculated Transition-State $\left( {\Delta {G}^{ \ddagger  }}\right)$ and Reaction $\left( {\Delta {G}_{\mathrm{{rxn}}}}\right)$ Gibbs Free Energies (kcal mol ${}^{-1}$ ) for Cyclic Grob Fragmentation of Cellobiose ${}^{a}$   

<table><tr><td rowspan="2">Reaction Label</td><td rowspan="2">Product Structure(s)</td><td colspan="2">25°C</td><td colspan="2">600°C</td></tr><tr><td>ΔG‡</td><td>ΔGrxn</td><td>ΔG‡</td><td>ΔGrxn</td></tr><tr><td>CB13Grob</td><td></td><td>64.4</td><td>22.6</td><td>63.4</td><td>-2.4</td></tr><tr><td>CB31Grob</td><td></td><td>77.4</td><td>0.2</td><td>79.7</td><td>-26.3</td></tr><tr><td>CB42Grob</td><td></td><td>81.1</td><td>22.6</td><td>81.8</td><td>-5.3</td></tr><tr><td>CBax46Grob</td><td></td><td>83.5</td><td>12.6</td><td>83.7</td><td>-37.1</td></tr><tr><td>CBeq46Grob</td><td></td><td>85.7</td><td>13.4</td><td>84.4</td><td>-37.4</td></tr><tr><td>CB53Grob</td><td></td><td>67.5</td><td>30.1</td><td>66.1</td><td>-3.3</td></tr><tr><td>CB1&#x27;3&#x27;Grob</td><td></td><td>67.2</td><td>29.2</td><td>65.7</td><td>-4.2</td></tr><tr><td>CB2&#x27;4&#x27;Grob</td><td></td><td>93.9</td><td>19.7</td><td>93.2</td><td>-8.8</td></tr><tr><td>CB4&#x27;2&#x27;Grob</td><td></td><td>82.2</td><td>21.0</td><td>83.0</td><td>-6.4</td></tr><tr><td>CB4&#x27;6&#x27;Grob</td><td></td><td>87.0</td><td>10.6</td><td>86.6</td><td>-36.3</td></tr><tr><td>CB5&#x27;3&#x27;Grob</td><td></td><td>61.6</td><td>28.2</td><td>61.8</td><td>-4.5</td></tr><tr><td>CBax6&#x27;4&#x27;Grob</td><td></td><td>75.8</td><td>19.8</td><td>76.9</td><td>-7.3</td></tr><tr><td>CBeq6&#x27;4&#x27;Grob</td><td></td><td>70.1</td><td>10.7</td><td>69.0</td><td>-18.5</td></tr></table>

${}^{a}$ The reaction labels indicate the reactant (CBeq = equatorial cellobiose; CBax = axial cellobiose),the respective carbon atoms at which a $\mathrm{C} - \mathrm{O}$ bond is broken and a $\mathrm{C} = \mathrm{O}$ bond is formed (numbering according to Figure 1),and the reaction type (Grob = cylic Grob fragmentation). Note that in the case of CB53Grob and CB5’ ${3}^{\prime }$ -Grob,transformation from an ether to a hydroxyl group occurs instead of $\mathrm{C} - \mathrm{O}$ bond cleavage.

$\mathrm{mol}^{-1}$ ,<sup>25</sup> ensuring that it is not the rate-determining step of this reaction pathway. Thus, the free-energy barriers for the cyclic Grob fragmentation reactions shown in Table 5 are all calculated relative to the lowest-energy equatorial conformation of glucose for comparison with other mechanisms.

Glucose can undergo cyclic Grob fragmentation via eight separate reactions (Table 5). Five of these reactions result in six-carbon acyclic products, while GLCeq46Grob, GLCax46-Grob, and GLC53Grob involve losses of one- or two-carbon fragments. For GLCeq46Grob and GLCax46Grob, a hydroxyl group is lost from C4 and a hydrogen atom from C6, which generates formaldehyde and a cyclic alkene that can further

undergo retro-Diels-Alder reaction. The free-energy barrier for GLCeq46Grob is $1.4\mathrm{kcalmol}^{-1}$ lower in energy than for GLCax46Grob, but both generate the same products. The overall lowest barrier among the studied cyclic Grob fragmentations is for GLC13Grob $(62.7\mathrm{kcalmol}^{-1})$ , and the highest barrier is that for GLC24Grob $(92.1\mathrm{kcalmol}^{-1})$ .The high variability in the Grob reaction barriers for glucose may be due to the partial charges on the carbon atoms in glucose. For example, the respective partial charges for carbon atoms C1 and C2 are 0.603 and 0.016 C. On the one hand, in the transition-state structure of reaction GLC13Grob, the anomeric carbon has a partial charge of 0.392 C, suggesting

![](images/cf4bc090f36e4162a5e87d0cf823f9aa6f26f7c10ad586093337142db8f79a4b.jpg)  
Figure 6. Depiction of the conformational changes that occur for cellobiose prior to dehydration via cyclic Grob fragmentation.

# Scheme 2. Generalized Aldol Rearrangement Mechanism

![](images/5f8f8b0ae7ed5d2f65d8d47996346f1a7ada270606eb9dae107329c8da3761d9.jpg)

that the transition-state structure benefits from additional electron donation from the ring oxygen to the anomeric carbon. On the other hand, in the transition state for reaction GLC24Grob, the charges on C1 and C2 are 0.637 and 0.078 C, which deviate from neutral glucose by 0.034 and 0.062 C, respectively. The similarity in charges between glucose and the transition state of GLC24Grob suggests that there is not a substantial stabilizing delocalization of charge as is observed with the transition-state structure of GLC13Grob. It is noteworthy that these barriers are similar in magnitude to those for Maccoll elimination.

In reaction GLC53Grob, a water molecule is not formed, since a hydrogen atom is transferred from C3 to the ring

oxygen to form a hydroxyl group. During this reaction, the O-C1 and C2-C3 bonds are broken to yield ethenediol, the tautomer of glycolaldehyde. However, this reaction, with a barrier of 68.6 kcal mol $^{-1}$ , is not competitive with the previously proposed reaction pathway for glycolaldehyde formation, $^{14}$ which requires only 49.2 kcal mol $^{-1}$ of free energy to proceed.

For cyclic Grob fragmentations of cellobiose, 13 separate reactions are reported here (Table 6). Nine of these result in the loss of a water molecule, while the other four result in the loss of glucose (CB42Grob and CB1'3'Grob) or glycolaldehyde (CB53Grob and CB5'3'Grob). Among those that lose water, the lowest barriers for reactions at the reducing and nonreducing ends, respectively, are CB13Grob (63.4 kcal mol $^{-1}$ ) and CBeq6'4'Grob (69.0 kcal mol $^{-1}$ ). Six of these reactions lead to a dehydrated cellobiose molecule, and the remaining three yield a low-molecular weight product (formaldehyde, ethenediol, or erythrose) and a glucoside. Note that, because the barrier for ring-flipping from equatorial to axial ring conformations is so low in glucose (12.1 kcal mol $^{-1}$ ), this conformational change was ignored for cellobiose—all free energies are calculated relative to the free energy of the cellobiose conformation described in the Computational Methods section. The conformations in Figure 6 show the ring-flipping intermediates that lead to the cyclic Grob fragmentation (i.e., the reducing end in the axial position for reactions at the reducing end and the nonreducing end in the axial position for reactions at the nonreducing end). As with the other water loss reactions, cyclic Grob fragmentation has free-energy barriers that are relatively independent of temperature. For those reaction barriers that increase at elevated temperatures, the transition-state geometries have shorter hydrogen bonds that lead to more constrained geometries and higher vibrational temperatures much like CB12Mac and CB34Mac, whose barrier heights also increase

![](images/ccc87b946e1144c48affbe4124c315750b5e067fe7b5af892c7a3e04de9f3934.jpg)  
Figure 7. Interconversion between Pinacol ring contraction products and cyclic Grob fragmentation products via aldol rearrangement. The labels next to the products/intermediates refer to the respective reactions in Table 3 and Table 5.

![](images/7521ba57b2cc0caf6de6aac4e7c043b17f1ef39fddb19dee0f921580b54ebef1.jpg)  
Figure 8. Calculated potential (free) energy surface $(\mathrm{kcal mol}^{-1})$ for the formation and interconversion of Pinacol ring contraction and cyclic Grob fragmentation products of glucose.

Table 7. Calculated Transition-State $(\Delta G^{\ddagger})$ Gibbs Free Energies (kcal mol $^{-1}$ ) for Transformations Between Pinacol Products and Grob Products Via Aldol Rearrangement   

<table><tr><td>reaction</td><td>ΔG‡ (forward)</td><td>ΔG‡ (reverse)</td></tr><tr><td>12Pin → 13Grob</td><td>37.3</td><td>25.6</td></tr><tr><td>21Pin → 13Grob</td><td>42.1</td><td>28.0</td></tr><tr><td>23Pin → 31Grob</td><td>50.3</td><td>46.3</td></tr><tr><td>23Pin → 24Grob</td><td>35.6</td><td>20.5</td></tr><tr><td>32Pin → 31Grob</td><td>34.5</td><td>28.0</td></tr><tr><td>32Pin → 24Grob</td><td>39.4</td><td>21.9</td></tr><tr><td>34Pin → 42Grob</td><td>41.9</td><td>21.2</td></tr><tr><td>43Pin → 42Grob</td><td>38.9</td><td>17.9</td></tr></table>

![](images/2635ed8f9f6b1645b90dfbc929a330ff34482cdec5b08f54f7c307c2a84156aa.jpg)  
Scheme 3. Generalized Alcohol Condensation Reaction Mechanism

at the elevated temperature due to a smaller vibrational entropic contribution to the Gibbs free energy.

Aldol Rearrangement. Some Pinacol and Grob reactions yield a product containing a carbonyl group that is in the $\gamma$ -position relative to a hydroxyl group (Table 3 and Table 5). These products can undergo aldol rearrangement as depicted in Scheme 2. For the products from reactions GLC12Pin and GLC21Pin, the hydroxyl group at C3 is in the $\gamma$ -position relative to the carbonyl group. After aldol rearrangement, the carbonyl group is at C3, and there is a double bond between C1 and C2. This compound corresponds to the product from the GLC13Grob reaction (Figure 7). The products of GLC23Pin and GLC32Pin (Figure 7) have two hydroxyl groups that are in the $\gamma$ -position relative to the carbonyl group. When these products undergo aldol rearrangement with the hydroxyl group at C1, the product is the same as that from the GLC24Grob reaction (Figure 7). These same products (i.e., GLC23Pin and GLC32Pin) can also undergo aldol rearrangement.

ment with the hydroxyl group at C4, resulting in the product from the GLC31Grob reaction. Finally, the carbonyl and the C3 hydroxyl groups in the products of the GLC34Pin and GLC43Pin reactions can react to yield the product of the GLC42Grob reaction. This reaction network is shown in Figure 7.

The potential (free) energy surface and free-energy barriers for the aldol rearrangement network are shown in Figure 8 and Table 7. Calculations on the transition states connecting the Pinacol products with the Grob products reveal that interconversion between these molecules is likely to occur at elevated temperatures. At $600^{\circ}\mathrm{C}$ , the highest free-energy barrier from a Pinacol product to a Grob product is $50.3\mathrm{kcalmol}^{-1}$ for the reaction from the GLC23Pin product to the GLC31Grob product, while the lowest barrier is $34.5\mathrm{kcalmol}^{-1}$ for the reaction from the GLC32Pin product to the GLC31Grob product.

Among the Pinacol and Grob reactions, GLC32Pin produces the molecule that is thermodynamically most stable $(-30.0\mathrm{kcalmol}^{-1}$ relative to glucose), with a free-energy barrier of $70.5\mathrm{kcalmol}^{-1}$ . In contrast, the transition state for GLC13Grob is $9.8\mathrm{kcalmol}^{-1}$ lower in energy, even though the Grob product is thermodynamically less stable $(-9.9\mathrm{kcalmol}^{-1}$ relative to glucose). However, this Grob product can readily convert to the thermodynamically more stable GLC12Pin product $(-21.6\mathrm{kcalmol}^{-1}$ relative to glucose) or GLC21Pin product $(-24.0\mathrm{kcalmol}^{-1}$ relative to glucose) with barriers of 25.6 and $28.0\mathrm{kcalmol}^{-1}$ , respectively.

Alcohol Condensation. The final class of dehydration reactions considered is generalized alcohol condensation, wherein an ether linkage replaces two hydroxyl groups (Scheme 3). The nature of these reactions is such that two separate transition states lead to the same product depending on which oxygen atom becomes the ether oxygen atom and which one leaves in the water molecule. An example of alcohol condensation is the well-known concerted levoglucosan formation via attack of the hydroxymethylene group of the glucose ring at the anomeric carbon to form a five-membered ring, which has been reported for glucose,[25] methyl $\beta$ -D-glucoside,[60] and methyl cellobiose.[11] In comparison with the

Table 8. Calculated Transition-State $\left( {\Delta {G}^{ \ddagger  }}\right)$ and Reaction $\left( {\Delta {G}_{\mathrm{{rxn}}}}\right)$ Gibbs Free Energies (kcal mol ${}^{-1}$ ) for Alcohol Condensation of Glucose ${}^{a}$   

<table><tr><td rowspan="2">Reaction Label</td><td rowspan="2">Product Structure(s)</td><td colspan="2">25°C</td><td colspan="2">600°C</td></tr><tr><td>ΔG‡</td><td>ΔGrxn</td><td>ΔG‡</td><td>ΔGrxn</td></tr><tr><td>GLC13Cond</td><td rowspan="2">HO
OH
OH</td><td>63.6</td><td>20.9</td><td>62.3</td><td>-0.2</td></tr><tr><td>GLC31Cond</td><td>83.9</td><td>20.9</td><td>82.0</td><td>-0.2</td></tr><tr><td>GLC16Cond</td><td rowspan="2">OH
OH
OH</td><td>48.1</td><td>-1.0</td><td>50.4</td><td>-17.9</td></tr><tr><td>GLC61Cond</td><td>76.4</td><td>-1.0</td><td>80.7</td><td>-17.9</td></tr><tr><td>GLC24Cond</td><td rowspan="2">HO
OH
OH</td><td>89.7</td><td>17.7</td><td>92.0</td><td>-2.4</td></tr><tr><td>GLC42Cond</td><td>80.3</td><td>17.7</td><td>78.9</td><td>-2.4</td></tr><tr><td>GLC36Cond</td><td rowspan="2">OH
OH
OH</td><td>89.0</td><td>2.3</td><td>93.8</td><td>-15.8</td></tr><tr><td>GLC63Cond</td><td>72.0</td><td>2.3</td><td>75.0</td><td>-15.8</td></tr><tr><td>GLC46Cond</td><td rowspan="2">OH
OH
OH</td><td>84.0</td><td>32.1</td><td>84.7</td><td>12.6</td></tr><tr><td>GLC64Cond</td><td>89.2</td><td>32.1</td><td>88.2</td><td>12.6</td></tr></table>

${}^{a}$ The reaction labels indicate the reactant $\left( {\mathrm{{GLC}} = \text{glucose}}\right)$ ,the respective carbon atoms at which a $\mathrm{C} - \mathrm{O}$ bond is broken and a $\mathrm{{CO}} - \mathrm{H}$ bond is broken (numbering according to Figure 1),and the reaction type (Cond = alcohol condensation). Note that each product can be formed in two ways depending on which OH group leaves as water and which OH group loses a hydrogen atom.

![](images/0a92256052c086199c1f01c0f819d0fb4bc0999840c38ebd7123d90e4eecd548.jpg)  
Figure 9. Calculated Gibbs free-energy barriers (kcal mol $^{-1}$ ) at 600 $^\circ \mathrm{C}$ for alcohol condensation reactions of cellobiose. The labels next to the reaction arrows indicate the reactant (CB = cellobiose; CBN = cellobiosan), the respective carbon atoms at which a C-O bond is broken and a C-O (ether) bond is formed (numbering according to Figure 1), and the reaction type (Cond = alcohol condensation).

# Scheme 4. Generalized Retro-Diels Alder Reaction

![](images/46fc3b3acd7bd47a5cb529cace1cf33fb50be8f643dd3b823af30ad71a56d96c.jpg)  
Each R group can be either carbon- or oxygen-centered.

calculated barrier of 50.4 for GLC16Cond in Table 8, Mayes et al. report an activation energy of $48.2\mathrm{kcalmol}^{-1}$ , while Seshadri et al. report activation energies of 46.1 and $53.2\mathrm{kcalmol}^{-1}$ for levoglucosan production from glucopyranose in the chair and boat conformations, respectively, calculated using the CBS-QB3 level of theory. The calculated barriers for reactions of this type are among the lowest reported here as well as in the literature on pyrolysis of cellulose model com

pounds, $^{11,25,60}$ which may help to explain the high yield of levoglucosan upon glucooligosaccharide pyrolysis. $^{4}$ However, this does not necessarily preclude other reactions from occurring during the pyrolysis of cellulose, since heat transfer effects can lead to thermal gradients, wherein higher barrier reactions may take place. $^{61}$ The incomplete conversion to levoglucosan and the production of char as a dehydration product suggest that at least some other dehydration pathways may be accessible during pyrolysis. $^{20,62,63}$

Free-energy barriers and reaction free energies for alcohol condensation of glucose are shown in Table 8. The free-energy barrier of $50.4\mathrm{kcalmol}^{-1}$ for the levoglucosan-producing GLC16Cond reaction makes it the most kinetically favorable water loss reaction discussed here. Levoglucosan formation via the loss of the oxygen atom on C1 (reaction GLC16Cond; Table 8) has a free-energy barrier that is $30.3\mathrm{kcalmol}^{-1}$ lower than the barrier for levoglucosan formation via the loss of the oxygen atom on C6 (reaction GLC61Cond; Table 8), which is in agreement with the literature.[25] Besides levoglucosan, another bicyclic five-membered ring product, formed by condensation of the C3 and C6 hydroxyl groups, as well as three other products that contain four-membered rings, can be formed by condensation of 1,3-diols (Table 8). Ring strain helps explain the large free-energy barriers for many of these reactions, with the exception of GLC13Cond, whose barrier is $62.3\mathrm{kcalmol}^{-1}$ (Table 8).

Aside from GLC16Cond, the high free-energy barriers and high reaction free energies for alcohol condensation reactions suggest that these reactions are less likely routes for loss of water than other dehydration mechanisms. Therefore, the only alcohol condensation reactions for cellobiose that were calculated are CB16Cond and CB1'6'Cond (Figure 9), which involve cleavage of a C-O bond at the anomeric carbon on the reducing end and nonreducing end, respectively, analogous to GLC16Cond. Following CB16Cond, an alcohol condensation of the resulting cellobiosan is also possible

Table 9. Calculated Transition-State $(\Delta G^{\ddagger})$ Gibbs Free Energies (kcal mol $^{-1}$ ) at $600^{\circ}\mathrm{C}$ for Retro-Diels Alder Reactions of Cyclic Alkene Fragments from Glucose Dehydration $^{a}$   

<table><tr><td>Reaction Label</td><td>Reactant Structure</td><td>Product Structures</td><td>ΔG‡</td></tr><tr><td>GLC12rDA</td><td>HO HO OH
HO HO
OH</td><td>HO CH=CH
HO CH=CH</td><td>47.1</td></tr><tr><td>GLC21rDA</td><td>HO HO
HO HO
OH</td><td>HO CH=CH
HO CH=CH</td><td>39.3</td></tr><tr><td>GLC23rDA</td><td>HO HO
HO HO
OH</td><td>HO CH=CH
HO CH=CH</td><td>43.1</td></tr><tr><td>GLC32rDA</td><td>HO HO
HO HO
OH</td><td>HO CH=CH
HO CH=CH</td><td>46.5</td></tr><tr><td>GLC34rDA</td><td>HO HO
HO HO
OH</td><td>HO CH=CH
HO CH=CH</td><td>47.2</td></tr><tr><td>GLC43rDA</td><td>HO HO
HO HO
OH</td><td>HO CH=CH
HO CH=CH</td><td>55.2</td></tr><tr><td>GLC45rDA</td><td>HO HO
HO HO
OH</td><td>HO CH=CH
HO CH=CH</td><td>43.1</td></tr><tr><td>GLC46rDA</td><td>HO HO
HO HO
OH</td><td>HO CH=CH
HO CH=CH</td><td>47.4</td></tr></table>

aThe reaction labels indicate the reactant $\mathrm{(GLC =}$ glucose), the respective carbon atoms at which an OH group and an H atom were lost during the previous dehydration reaction (numbering according to Figure 1), and the reaction type $\mathrm{(rDA =}$ retro-Diels Alder). Note that, for GLC46rDA, water loss occurred via GLC46Grob; for all other cases, water loss occurred via the corresponding Maccoll elimination reaction (e.g., GLC12rDA was preceded by GLC12Mac).

(CBN1'6'Cond; Figure 9). The barriers for the conformational change that must precede CB16Cond and $\mathrm{CBN1^{\prime}6^{\prime}Cond}$ are 10.4 and $9.2\mathrm{kcalmol^{-1}}$ , respectively, implying that these conformational changes likely occur at a much faster rate than the corresponding condensation reactions (see Supporting Information for potential energy surface). The barrier for $\mathrm{CB1^{\prime}6^{\prime}Cond}$ reported here (56.2 kcal $\mathrm{mol^{-1}}$ ) is similar to the barrier for methyl-cellobiose (54.4 kcal $\mathrm{mol^{-1}}$ ). Notably, the products from $\mathrm{CB1^{\prime}6^{\prime}Cond}$ do not agree with a previous experimental study on selectively $^{13}\mathrm{C}$ labeled cellobioses, because the glucose formed during fast pyrolysis of these molecules was demonstrated to originate exclusively from the nonreducing end. However, the yields of levoglucosans originating from each glucose ring in cellobiose are the same if cellobiose undergoes sequential alcohol condensations at C1 and C6 of each glucose ring, which is consistent with previous experimental results. This finding suggests that there may be something prohibitive to $1^{\prime}6^{\prime}$ condensation of cellobiose relative to other possible reactions.

Retro-Diels-Alder Reaction. Cyclic alkenes that are formed from Maccoll elimination and cyclic Grob reactions

can undergo a wide variety of subsequent reactions. The retro-Diels-Alder (rDA) reaction (Scheme 4) has been suggested to be important during cellulose pyrolysis.[27,64] The transition-state structures and reaction free energies were calculated for the rDA reaction of the Maccoll elimination products of glucose and cellobiose (Table 9 and Table 10).

Since the rDA reactions involve a highly conjugated system, the delocalization of electrons across the $\pi$ -bonds is an important contributor to the transition-state energies. The largest overlap of $\pi$ -orbitals occurs for atoms that are covalently bound, and there is less orbital overlap between the atoms whose covalent bonds cleave during the reaction, which interrupts the conjugation in the system. Interestingly, the rDA reactions that produce a one-carbon fragment possess some of the largest barriers: GLC34rDA has a barrier of 47.2 kcal mol $^{-1}$ , and GLC43rDA has a barrier of 55.2 kcal mol $^{-1}$ . The rDA reactions that produce a two-carbon fragment possess barriers that are lower than those that produce a one-carbon fragment: GLC23rDA has a barrier of 43.1 kcal mol $^{-1}$ , GLC32rDA has a barrier of 46.5 kcal mol $^{-1}$ , and GLC45rDA has a barrier of 43.1 kcal mol $^{-1}$ . Finally, the rDA reaction with the lowest barrier is GLC21rDA (39.3 kcal mol $^{-1}$ ), which produces two fragments each with three carbons. This suggests a possible trend between the stabilizing effect of the conjugation and the size of the resulting fragments. An exception to this trend is reaction GLC12rDA, which also produces two fragments with three carbons each, but has a reaction barrier of 47.1 kcal mol $^{-1}$ . In this case, the stabilizing effect of conjugation is dampened by the interaction of the two electron-donating hydroxyl groups that are $\alpha$ to each other. By comparison, reaction GLC21rDA has an electron-donating hydroxyl group that is $\beta$ to an electron-withdrawing acid group, which further enhances the conjugation. The effect of $\alpha$ -hydroxyl pairs on the reaction barrier can also be observed for the two reactions producing five-carbon fragments: the reaction with the higher transition-state energy (GLC43rDA) produces a fragment with two hydroxyl groups $\alpha$ to each other, while the other reaction (GLC34rDA) does not. Both reactions that produce four-carbon fragments (GLC23rDA and GLC32rDA) have hydroxyl groups that are $\alpha$ to each other and differ in barrier heights by only 3.4 kcal mol $^{-1}$ .

Cyclic alkenes arising from dehydration are likely to fragment further through a retro-Diels-Alder reaction, since the calculated barriers for rDA are relatively low—between 37.1 and $55.2\mathrm{kcalmol}^{-1}$ (Table 9 and Table 10)—compared with the dehydration reactions that precede them. The discrepancies between the barriers for corresponding reactions of glucose and cellobiose are due to the extent of hydrogen bonding present in the transition-state structure as was shown to be the case for Maccoll elimination. The additional glucose ring in cellobiose has little impact on the overall barriers for rDA compared to the corresponding glucose reactions beyond the change in the hydrogen bonding. An increase in the hydrogen bonding is associated with a lower barrier. However, a large fiber of cellulose may be less likely to undergo rDA if it is constrained in the crystal phase by its hydrogen-bonding network. As the fiber liquefies and depolymerizes during fast pyrolysis, these reactions may be important to the overall product distribution.

# CONCLUSIONS

Both glucose and cellobiose can lose water via several mechanisms, including Maccoll elimination, Pinacol rearrange-

Table 10. Calculated Transition-State $\left( {\Delta {G}^{ \ddagger  }}\right)$ Gibbs Free Energies (kcal mol ${}^{-1}$ ) at ${600}{}^{ \circ  }\mathrm{C}$ for Retro-Diels Alder Reactions of Cyclic Alkene Fragments from Cellobiose Dehydration ${}^{a}$   

<table><tr><td>Reaction Label</td><td>Reactant Structure</td><td>Product Structures</td><td>ΔG‡</td></tr><tr><td>CB12rDA</td><td>HO HO OH HO OH OH</td><td>HO HO OH HO OH</td><td>41.7</td></tr><tr><td>CB21rDA</td><td>HO HO OH HO OH OH</td><td>HO HO OH HO OH</td><td>40.2</td></tr><tr><td>CB23rDA</td><td>HO HO OH HO OH OH</td><td>HO HO OH HO OH</td><td>40.3</td></tr><tr><td>CB32rDA</td><td>HO HO OH HO OH OH</td><td>HO HO OH HO OH</td><td>53.7</td></tr><tr><td>CB34rDA</td><td>HO HO OH HO OH OH</td><td>HO HO OH HO OH</td><td>51.2</td></tr><tr><td>CB2&#x27;1&#x27;rDA</td><td>HO HO OH HO OH OH</td><td>HO HO OH HO OH</td><td>37.1</td></tr><tr><td>CB2&#x27;3&#x27;rDA</td><td>HO HO OH HO OH OH</td><td>HO HO OH HO OH</td><td>38.9b</td></tr><tr><td>CB3&#x27;2&#x27;rDA</td><td>HO HO OH HO OH OH</td><td>HO HO OH HO OH</td><td>47.8</td></tr><tr><td>CB3&#x27;4&#x27;rDA</td><td>HO HO OH HO OH OH</td><td>HO HO OH HO OH</td><td>48.4</td></tr><tr><td>CB4&#x27;3&#x27;rDA</td><td>HO HO OH HO OH OH</td><td>HO HO OH HO OH</td><td>43.9</td></tr><tr><td>CB4&#x27;5&#x27;rDA</td><td>HO HO OH HO OH OH</td><td>HO HO OH HO OH</td><td>45.0b</td></tr></table>

${}^{a}$ The reaction labels indicate the reactant $\left( {\mathrm{{CB}} = \text{cellobiose}}\right)$ ,the respective carbon atoms at which an $\mathrm{{OH}}$ group and an $\mathrm{H}$ atom were lost during the previous dehydration reaction (numbering according to Figure 1),and the reaction type (rDA = retro-Diels Alder). Note that,in all cases,water loss occurred via the corresponding Maccoll elimination reaction (e.g., CB12rDA was preceded by CB12Mac). Two of the rDA transition states (CB2′3′rDA and CB4′5′rDA) were found not to proceed via a one-step concerted mechanism but rather via a two-step mechanism with a loosely bound zwitterionic intermediate. ${}^{b}$ Barrier heights for rate-limiting step in two-step mechanism.

ment, cyclic Grob fragmentation, and alcohol condensation. The dehydration reactions of glucose and cellobiose under fast pyrolysis conditions (calculated at $600^{\circ}\mathrm{C}$ ) have barriers ranging from 50 to 94 kcal mol $^{-1}$ . For glucose, the lowest barriers for each water-loss mechanism are 65.9 kcal mol $^{-1}$ for Maccoll elimination, 67.1 kcal mol $^{-1}$ for Pinacol ring contraction, 62.7 kcal mol $^{-1}$ for cyclic Grob fragmentation, and 50.4 kcal mol $^{-1}$ for alcohol condensation. Levoglucosan

formation via alcohol condensation is therefore the most kinetically favored dehydration reaction for glucose. A new reaction network of aldol rearrangements unlocks the interconversion of Pinacol rearrangement and cyclic Grob fragmentation products. The energies for the transition states of these aldol rearrangements are as low as $4.5\mathrm{kcalmol}^{-1}$ relative to glucose. For cellobiose, the lowest barriers for each water-loss mechanism are $62.5\mathrm{kcalmol}^{-1}$ for Maccoll

elimination, $61.5\mathrm{kcalmol}^{-1}$ for Pinacol ring contraction, and $61.8\mathrm{kcalmol}^{-1}$ for cyclic Grob fragmentation. Retro-Diels-Alder reactions on the dehydration products of glucose and cellobiose have much lower barriers than the dehydration reactions themselves: $39.3\mathrm{kcalmol}^{-1}$ for dehydrated glucose rDA and $37.1\mathrm{kcalmol}^{-1}$ for dehydrated cellobiose rDA. Thus, if the water-loss products are formed during fast pyrolysis, it is likely that they continue to react via this mechanism. It is also important to note that the reaction barriers are relatively independent of temperature (i.e., 25 vs $600^{\circ}\mathrm{C}$ ) and the reactions are more exergonic at higher temperatures, as expected. The similarity of the reactions for glucose and cellobiose suggests that these reactions should provide insights into the overall fast pyrolysis process for oligosaccharides and cellulose.

# ASSOCIATED CONTENT

# Supporting Information

The Supporting Information is available free of charge on the ACS Publications website at DOI: 10.1021/acs.jpca.8b02312.

Full Gaussian reference, full potential free-energy surface for levoglucosan formation from cellobiose, and XYZ coordinates of transition states and minima (PDF)

# AUTHOR INFORMATION

ORCID

Mckay W. Easton: 0000-0003-2754-4799

Hilkka I. Kenttämaa: 0000-0001-8988-6984

Notes

The authors declare no competing financial interest.

# ACKNOWLEDGMENTS

This work was supported as part of the Center for Direct Catalytic Conversion of Biomass to Biofuels (C3Bio), an Energy Frontier Research Center (EFRC) funded by the U.S. Department of Energy, Office of Science, Office of Basic Energy Sciences under Award No. DE-SC0000997.

# REFERENCES

(1) Mohan, D.; Pittman, C. U.; Steele, P. H. Pyrolysis of Wood/Biomass for Bio-Oil: A Critical Review. Energy Fuels 2006, 20, 848-889.   
(2) Mettler, M. S.; Paulsen, A. D.; Vlachos, D. G.; Dauenhauer, P. J. Pyrolytic Conversion of Cellulose to Fuels: Levoglucosan Deoxygenation via Elimination and Cyclization within Molten Biomass. Energy Environ. Sci. 2012, 5, 7864-7868.   
(3) Agrawal, R.; Singh, N. R. Synergistic Routes to Liquid Fuel for a Petroleum-Deprived Future. AIChE J. 2009, 55, 1898–1905.   
(4) Venkatakrishnan, V. K.; Degenstein, J. C.; Smeltz, A. D.; Delgass, W. N.; Agrawal, R.; Ribeiro, F. H. High-Pressure Fast-Pyrolysis, Fast-Hydroxylation and Catalytic Hydrodeoxygenation of Cellulose: Production of Liquid Fuel from Biomass. Green Chem. 2014, 16, 792-802.   
(5) Wang, W.; Shi, Y.; Cui, Y.; Li, X. Catalytic Fast Pyrolysis of Cellulose for Increasing Contents of Furans and Aromatics in Biofuel Production. J. Anal. Appl. Pyrolysis 2018, 131, 93-100.   
(6) Mettler, M. S.; Vlachos, D. G.; Dauenhauer, P. J. Top Ten Fundamental Challenges of Biomass Pyrolysis for Biofuels. Energy Environ. Sci. 2012, 5, 7797-7809.   
(7) Brown, T. R.; Wright, M. M.; Brown, R. C. Estimating Profitability of Two Biochar Production Scenarios: Slow Pyrolysis vs Fast Pyrolysis. Biofuels, Bioprod. Biorefin. 2011, S, 54-68.   
(8) Parsell, T.; Yohe, S.; Degenstein, J.; Jarrell, T.; Klein, I.; Gencer, E.; Hewetson, B.; Hurt, M.; Kim, J. I.; Choudhari, H.; et al. A

Synergistic Biorefinery Based on Catalytic Conversion of Lignin Prior to Cellulose Starting from Lignocellulosic Biomass. Green Chem. 2015, 17, 1492-1499.   
(9) Bradbury, A. G. W.; Sakai, Y.; Shafizadeh, F. A Kinetic Model for Pyrolysis of Cellulose. J. Appl. Polym. Sci. 1979, 23, 3271-3280.   
(10) Vinu, R.; Broadbelt, L. J. A Mechanistic Model of Fast Pyrolysis of Glucose-Based Carbohydrates to Predict Bio-Oil Composition. Energy Environ. Sci. 2012, 5, 9808-9826.   
(11) Mayes, H. B.; Broadbelt, L. J. Unraveling the Reactions That Unravel Cellulose. J. Phys. Chem. A 2012, 116, 7098-7106.   
(12) Liang, X.; Montoya, A.; Haynes, B. S. Local Site Selectivity and Conformational Structures in the Glycosidic Bond Scission of Cellobiose. J. Phys. Chem. B 2011, 115, 10682-10691.   
(13) Agarwal, V.; Dauenhauer, P. J.; Huber, G. W.; Auerbach, S. M. Ab Initio Dynamics of Cellulose Pyrolysis: Nascent Decomposition Pathways at 327 and $600^{\circ}\mathrm{C}$ . J. Am. Chem. Soc. 2012, 134, 14958-14972.   
(14) Degenstein, J. C.; Murria, P.; Easton, M.; Sheng, H.; Hurt, M.; Dow, A. R.; Gao, J.; Nash, J. J.; Agrawal, R.; Delgass, W. N.; et al. Fast Pyrolysis of 13C-Labeled Cellobioses: Gaining Insights into the Mechanisms of Fast Pyrolysis of Carbohydrates. J. Org. Chem. 2015, 80, 1909-1914.   
(15) Degenstein, J. C.; Hurt, M.; Murria, P.; Easton, M.; Choudhari, H.; Yang, L.; Riedeman, J.; Carlsen, M. S.; Nash, J. J.; Agrawal, R.; et al. Mass Spectrometric Studies of Fast Pyrolysis of Cellulose. Eur. J. Mass Spectrom. 2015, 21, 321-326.   
(16) Hurt, M. R.; Degenstein, J. C.; Gawecki, P.; Morton II, D. J.; Vinueza, N. R.; Yang, L.; Agrawal, R.; Delgass, W. N.; Ribeiro, F. H.; Kenttämaa, H. I. On-Line Mass Spectrometric Methods for the Determination of the Primary Products of Fast Pyrolysis of Carbohydrates and for Their Gas-Phase Manipulation. Anal. Chem. 2013, 85, 10927-10934.   
(17) Wang, S.; Guo, X.; Liang, T.; Zhou, Y.; Luo, Z. Mechanism Research on Cellulose Pyrolysis by Py-GC/MS and Subsequent Density Functional Theory Studies. Bioresour. Technol. 2012, 104, 722-728.   
(18) Patwardhan, P. R.; Satrio, J. A.; Brown, R. C.; Shanks, B. H. Product Distribution from Fast Pyrolysis of Glucose-Based Carbohydrates. J. Anal. Appl. Pyrolysis 2009, 86, 323-330.   
(19) Mayes, H. B.; Nolte, M. W.; Beckham, G. T.; Shanks, B. H.; Broadbelt, L. J. The Alpha-Bet(a) of Glucose Pyrolysis: Computational and Experimental Investigations of 5-Hydroxymethylfurfural and Levoglucosan Formation Reveal Implications for Cellulose Pyrolysis. ACS Sustainable Chem. Eng. 2014, 2, 1461-1473.   
(20) Xin, S.; Yang, H.; Chen, Y.; Yang, M.; Chen, L.; Wang, X.; Chen, H. Chemical Structure Evolution of Char during the Pyrolysis of Cellulose. J. Anal. Appl. Pyrolysis 2015, 116, 263-271.   
(21) Lin, X.; Qu, Y.; Lv, Y.; Xi, Y.; Phillips, D. L.; Liu, C. The First Dehydration and the Competing Reaction Pathways of Glucose Homogeneously and Heterogeneously Catalyzed by Acids. Phys. Chem. Chem. Phys. 2013, 15, 2967-2982.   
(22) Nimlos, M. R.; Blanksby, S. J.; Ellison, G. B.; Evans, R. J. Enhancement of 1,2-Dehydration of Alcohols by Alkali Cations and Protons: A Model for Dehydration of Carbohydrates. J. Anal. Appl. Pyrolysis 2003, 66, 3-27.   
(23) Zhang, M.; Geng, Z.; Yu, Y. Density Functional Theory (DFT) Study on the Pyrolysis of Cellulose: The Pyran Ring Breaking Mechanism. Comput. Theor. Chem. 2015, 1067, 13-23.   
(24) Paine, J. B.; Pithawalla, Y. B.; Naworal, J. D. Carbohydrate Pyrolysis Mechanisms from Isotopic Labeling: Part 2. The Pyrolysis of d-Glucose: General Disconnective Analysis and the Formation of C1 and C2 Carbonyl Compounds by Electrocyclic Fragmentation Mechanisms. J. Anal. Appl. Pyrolysis 2008, 82, 10-41.   
(25) Seshadri, V.; Westmoreland, P. R. Concerted Reactions and Mechanism of Glucose Pyrolysis and Implications for Cellulose Kinetics. J. Phys. Chem. A 2012, 116, 11997-12013.   
(26) Matsuoka, S.; Kawamoto, H.; Saka, S. Retro-Aldol-Type Fragmentation of Reducing Sugars Preferentially Occurring in

Polyether at High Temperature: Role of the Ether Oxygen as a Base Catalyst. J. Anal. Appl. Pyrolysis 2012, 93, 24-32.   
(27) Zhou, X.; Nolte, M. W.; Mayes, H. B.; Shanks, B. H.; Broadbelt, L. J. Experimental and Mechanistic Modeling of Fast Pyrolysis of Neat Glucose-Based Carbohydrates. 1. Experiments and Development of a Detailed Mechanistic Model. Ind. Eng. Chem. Res. 2014, 53, 13274-13289.   
(28) Assary, R. S.; Curtiss, L. A. Thermochemistry and Reaction Barriers for the Formation of Levoglucosenone from Cellobiose. ChemCatChem 2012, 4, 200-205.   
(29) Zhang, Y.; Liu, C.; Chen, X. Unveiling the Initial Pyrolytic Mechanisms of Cellulose by DFT Study. J. Anal. Appl. Pyrolysis 2015, 113, 621-629.   
(30) Assary, R. S.; Curtiss, L. A. Comparison of Sugar Molecule Decomposition through Glucose and Fructose: A High-Level Quantum Chemical Study. Energy Fuels 2012, 26, 1344-1352.   
(31) Maccoll, A. Heterolysis and Pyrolysis of Alkyl Halides in Gas Phase. Chem. Rev. 1969, 69, 33-60.   
(32) Collins, C. J. The Pinacol Rearrangement. Q. Rev., Chem. Soc. 1960, 14, 357-377.   
(33) Grob, C. A.; Baumann, W. Die 1,4-Eliminierung Unter Fragmentierung. Helv. Chim. Acta 1955, 38, 594-610.   
(34) Hosoya, T.; Sakaki, S. Levoglucosan Formation from Crystalline Cellulose: Importance of a Hydrogen Bonding Network in the Reaction. ChemSusChem 2013, 6, 2356-2368.   
(35) MacroModel, Version 10.3; Schrödinger, LLC: New York, NY, 2014.   
(36) Chang, G.; Guida, W. C.; Still, W. C. An Internal-Coordinate Monte Carlo Method for Searching Conformational Space. J. Am. Chem. Soc. 1989, 111, 4379-4386.   
(37) Halgren, T. A. MMFF VI. MMFF94s Option for Energy Minimization Studies. J. Comput. Chem. 1999, 20, 720-729.   
(38) Full reference in Supporting Information.   
(39) Schnupf, U.; Momany, F. A. Rapidly Calculated DFT Relaxed Iso-Potential $\phi/\psi$ Maps: $\beta$ -Cellobiose. Cellulose 2011, 18, 859-887.   
(40) French, A. D.; Johnson, G. P.; Cramer, C. J.; Csonka, G. I. Conformational Analysis of Cellobiose by Electronic Structure Theories. Carbohydr. Res. 2012, 350, 68-76.   
(41) Barrows, S. E.; Dulles, F. J.; Cramer, C. J.; French, A. D.; Truhlar, D. G. Relative Stability of Alternative Chair Forms and Hydroxymethyl Conformations of $\beta$ -d-Glucopyranose. Carbohydr. Res. 1995, 276, 219-251.   
(42) Dowd, M. K.; French, A. D.; Reilly, P. J. Conformational Analysis of the Anomeric Forms of Sophorose, Laminarabiose, and Cellobiose Using MM3. Carbohydr. Res. 1992, 233, 15-34.   
(43) Nishiyama, Y.; Langan, P.; Chanzy, H. Crystal Structure and Hydrogen-Bonding System in Cellulose 1 Beta from Synchrotron X-Ray and Neutron Fiber Diffraction. J. Am. Chem. Soc. 2002, 124, 9074-9082.   
(44) Jacobson, R. A.; Wunderlich, J. A.; Lipscomb, W. N. The Crystal and Molecular Structure of Cellobiose. Acta Crystallogr. 1961, 14, 598-607.   
(45) Dauenhauer, P. J.; Colby, J. L.; Balonek, C. M.; Suszynski, W. J.; Schmidt, L. D. Reactive Boiling of Cellulose for Integrated Catalysis through an Intermediate Liquid. Green Chem. 2009, 11, 1555-1561.   
(46) Teixeira, A. R.; Mooney, K. G.; Kruger, J. S.; Williams, C. L.; Suszynski, W. J.; Schmidt, L. D.; Schmidt, D. P.; Dauenhauer, P. J. Aerosol Generation by Reactive Boiling Ejection of Molten Cellulose. Energy Environ. Sci. 2011, 4, 4306-4321.   
(47) Becke, A. D. Density-Functional Thermochemistry. III. The Role of Exact Exchange. J. Chem. Phys. 1993, 98, 5648.   
(48) Lee, C.; Yang, W.; Parr, R. G. Development of the Colle-Salvetti Correlation-Energy Formula into a Functional of the Electron Density. Phys. Rev. B: Condens. Matter Mater. Phys. 1988, 37, 785-789.   
(49) Vosko, S. H.; Wilk, L.; Nusair, M. Accurate Spin-Dependent Electron Liquid Correlation Energies for Local Spin Density Calculations: A Critical Analysis. Can. J. Phys. 1980, 58, 1200-1211.

(50) Stephens, P. J.; Devlin, F. J.; Chabalowski, C. F.; Frisch, M. J. Ab Initio Calculation of Vibrational Absorption and Circular Dichroism Spectra Using Density Functional Force Fields. J. Phys. Chem. 1994, 98, 11623-11627.   
(51) Sladkovičová, M.; Mach, P.; Smrčok, L.; Rundlöf, H. DFT and Neutron Diffraction Study of 1,6-Anhydro-β-D-Glucopyranose (Levoglucosan). Cent. Eur. J. Chem. 2007, S, 55–70.   
(52) Lynch, B. J.; Truhlar, D. G. How Well Can Hybrid Density Functional Methods Predict Transition State Geometries and Barrier Heights? J. Phys. Chem. A 2001, 10S, 2936-2941.   
(53) Zhao, Y.; González-García, N.; Truhlar, D. G. Benchmark Database of Barrier Heights for Heavy Atom Transfer, Nucleophilic Substitution, Association, and Unimolecular Reactions and Its Use to Test Theoretical Methods. J. Phys. Chem. A 2005, 109, 2012-2018.   
(54) Zhao, Y.; Truhlar, D. G. The M06 Suite of Density Functionals for Main Group Thermochemistry, Thermochemical Kinetics, Noncovalent Interactions, Excited States, and Transition Elements: Two New Functionals and Systematic Testing of Four M06-Class Functionals and 12 Other Functionals. Theor. Chem. Acc. 2008, 120, 215-241.   
(55) Strati, G. L.; Willett, J. L.; Momany, F. A. DFT/Ab Initio Study of Hydrogen Bonding and Conformational Preference in Model Cellobiose Analogs Using B3LYP/6-311++G**. Carbohydr. Res. 2002, 337, 1851-1859.   
(56) Breneman, C. M.; Wiberg, K. B. Determining Atom-Centered Monopoles from Molecular Electrostatic Potentials. The Need for High Sampling Density in Formamide Conformational Analysis. J. Comput. Chem. 1990, 11, 361-373.   
(57) Jarvis, M. W.; Daily, J. W.; Carstensen, H.-H.; Dean, A. M.; Sharma, S.; Dayton, D. C.; Robichaud, D. J.; Nimlos, M. R. Direct Detection of Products from the Pyrolysis of 2-Phenethyl Phenyl Ether. J. Phys. Chem. A 2011, 115, 428-438.   
(58) Zhang, X.; Li, J.; Yang, W.; Blasiak, W. Formation Mechanism of Levoglucosan and Formaldehyde during Cellulose Pyrolysis. Energy Fuels 2011, 25, 3739-3746.   
(59) Stubbs, J. M.; Marx, D. Aspects of Glycosidic Bond Formation in Aqueous Solution: Chemical Bonding and the Role of Water. Chem. - Eur. J. 2005, 11, 2651-2659.   
(60) Hosoya, T.; Nakao, Y.; Sato, H.; Kawamoto, H.; Sakaki, S. Thermal Degradation of Methyl $\beta$ -d-Glucoside. A Theoretical Study of Plausible Reaction Mechanisms. J. Org. Chem. 2009, 74, 6891-6894.   
(61) Okekunle, P. O.; Pattanotai, T.; Watanabe, H.; Okazaki, K. Numerical and Experimental Investigation of Intra-Particle Heat Transfer and Tar Decomposition during Pyrolysis of Wood Biomass. J. Therm. Sci. Technol. 2011, 6, 360-375.   
(62) Broido, A.; Nelson, M. A. Char Yield on Pyrolysis of Cellulose. Combust. Flame 1975, 24, 263-268.   
(63) Ball, R.; McIntosh, A. C.; Brindley, J. The Role of Char-Forming Processes in the Thermal Decomposition of Cellulose. Phys. Chem. Chem. Phys. 1999, 1, 5035-5043.   
(64) Richards, G. N. Glycolaldehyde from Pyrolysis of Cellulose. J. Anal. Appl. Pyrolysis 1987, 10, 251-255.