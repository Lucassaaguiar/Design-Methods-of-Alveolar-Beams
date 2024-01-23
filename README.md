[EN](#Design-Methods-of-Alveolar-Beams)   |   [PT](#Dimensionamento-de-Vigas-Alveolares)

# Design Methods of Alveolar Beams
Computational implementation for design methods of cellular steel beams, aiming to identify failure modes and determine critical failure loads. The formulation is based on the works of Veríssimo et al. (2012), Fares et al. (2016), Ward (1990), Grilo et al. (2018), and Annex N (1998).

The codes were developed using Python version 3.8.8, and the program structure adheres to the principles of Object-Oriented Programming (OOP). These codes are organized into two distinct modules: Design Procedures of Alveolar Beams (DPAB) and Numerical Analysis of Alveolar Beams (NAAB). The DPAB is responsible for declaring the class, defining its attributes, and specifying its methods. Meanwhile, the NAAB is tasked with creating the object (alveolar beams) and conducting the numerical analysis. This involves incrementing loads and employing the relevant attributes and methods associated with that specific object.

For more information, see [Aguiar (2023)](https://lume.ufrgs.br/handle/10183/259192#).

# Examples

Examples | $d$ [mm] | $t_w$ [mm] | $t_f$ [mm] | $b_f$ [mm] | $b_w$ [mm] | $d_g$ [mm] | $L$ [m] | $D_0$ [mm] | $f_y$ [MPa] | $E$ [GPa] | $\nu$
--- | --- | --- | --- |--- |--- |--- |--- |--- |--- |--- |--- 
01 | 30.30 | 0.51 | 0.57 | 10.10 | 8.18 | 45.45 | 4.55 | 27.27 | 345 | 200 | 0.30 
02 | 40.30 | 0.70 | 1.12 | 14.00 | 8.06 | 60.45 | 12.09 | 40.30 | 345 | 200 | 0.30 
03 | 52.90 | 0.97 | 1.36 | 16.60 | 31.74 | 79.35 | 3.97 | 63.48 | 345 | 200 | 0.30 

# References 

AGUIAR, L. A. Análise paramétrica sobre procedimentos para dimensionamento de vigas alveolares de aço. 180 p. Thesis (master) — Federal University of Rio Grande do Sul, 2023.

ASSOCIAÇÃO BRASILEIRA DE NORMAS TÉCNICAS (ABNT). NBR 8800: Projeto de estruturas de aço e de estruturas mistas de aço e concreto de edifícios. Rio de Janeiro, 2008. 247 p.

CEN - EUROPEAN COMMITTEE FOR STANDARDIZATION. ENV 1993-1-1. Eurocode 3: Design of Steel Structures - Part 1-1 - General rules and rules for building. Amendment A2: N: Openings in webs (draft). Europe, 1998.

FARES, S.; COULSON, J. DINEHART, D. AISC design guide 31: Castellated and cellular beam design, American Institute of Steel Construction, USA, 2016.

GRILO, L. F.; FAKURY, R. H.; VERÍSSIMO, G. S. Design procedure for the web-post buckling of steel cellular beams. Journal of Constructional Steel Research, v. 148, p. 525–541, 2018.

VERÍSSIMO, G. S.; VIEIRA, W. B.; SILVEIRA, E. G.; RIBEIRO, J. C. L.; PAES, J. L. R.; BEZERRA, E. M.; Castro e Silva, A. L. R. Estados limites aplicáveis às vigas alveolares de aço. Revista de Estrutura de Aço, v. 2, n. 2, p. 126-144, 2013.

WARD, J. K. Design of Composite and Non-Composite Cellular Beams: Publication nº100. Ascot, UK, 1990.

---

# Dimensionamento de Vigas Alveolares

A implementação computacional foi realizada para dimensionar e verificar vigas alveolares de aço, com o objetivo de identificar modos e cargas críticas de falha. A formulação é baseada nos trabalhos de Veríssimo et al. (2012), Fares et al. (2016), Ward (1990), Grilo et al. (2018) e Anexo N (1998).

Os códigos foram desenvolvidos utilizando a versão 3.8.8 do Python, e a estrutura do programa adere aos princípios da Programação Orientada a Objetos (OOP). Esses códigos estão organizados em dois módulos distintos: Design Procedures of Alveolar Beams (DPAB) e Numerical Analysis of Alveolar Beams (NAAB). O DPAB é responsável por declarar a classe, definir seus atributos e especificar seus métodos. Enquanto isso, o NAAB é encarregado de criar o objeto (vigas alveolares) e realizar a análise numérica. Isso envolve o incremento de cargas e a utilização dos atributos e métodos relevantes associados a esse objeto específico.

Para mais informações ver em [Aguiar (2023)](https://lume.ufrgs.br/handle/10183/259192#).

# Exemplos

Exemplos | $d$ [mm] | $t_w$ [mm] | $t_f$ [mm] | $b_f$ [mm] | $b_w$ [mm] | $d_g$ [mm] | $L$ [m] | $D_0$ [mm] | $f_y$ [MPa] | $E$ [GPa] | $\nu$
--- | --- | --- | --- |--- |--- |--- |--- |--- |--- |--- |--- 
01 | 30.30 | 0.51 | 0.57 | 10.10 | 8.18 | 45.45 | 4.55 | 27.27 | 345 | 200 | 0.30 
02 | 40.30 | 0.70 | 1.12 | 14.00 | 8.06 | 60.45 | 12.09 | 40.30 | 345 | 200 | 0.30 
03 | 52.90 | 0.97 | 1.36 | 16.60 | 31.74 | 79.35 | 3.97 | 63.48 | 345 | 200 | 0.30 

# Referências 

AGUIAR, L. A. Análise paramétrica sobre procedimentos para dimensionamento de vigas alveolares de aço. 180 p. Thesis (master) — Federal University of Rio Grande do Sul, 2023.

ASSOCIAÇÃO BRASILEIRA DE NORMAS TÉCNICAS (ABNT). NBR 8800: Projeto de estruturas de aço e de estruturas mistas de aço e concreto de edifícios. Rio de Janeiro, 2008. 247 p.

CEN - EUROPEAN COMMITTEE FOR STANDARDIZATION. ENV 1993-1-1. Eurocode 3: Design of Steel Structures - Part 1-1 - General rules and rules for building. Amendment A2: N: Openings in webs (draft). Europe, 1998.

FARES, S.; COULSON, J. DINEHART, D. AISC design guide 31: Castellated and cellular beam design, American Institute of Steel Construction, USA, 2016.

GRILO, L. F.; FAKURY, R. H.; VERÍSSIMO, G. S. Design procedure for the web-post buckling of steel cellular beams. Journal of Constructional Steel Research, v. 148, p. 525–541, 2018.

VERÍSSIMO, G. S.; VIEIRA, W. B.; SILVEIRA, E. G.; RIBEIRO, J. C. L.; PAES, J. L. R.; BEZERRA, E. M.; Castro e Silva, A. L. R. Estados limites aplicáveis às vigas alveolares de aço. Revista de Estrutura de Aço, v. 2, n. 2, p. 126-144, 2013.

WARD, J. K. Design of Composite and Non-Composite Cellular Beams: Publication nº100. Ascot, UK, 1990.




